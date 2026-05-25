#!/usr/bin/env python3
"""
Generate CSV and interactive HTML grid heatmap for binned ATAC fragments.
Color scale uses log10(nFrags+1), hover shows raw integer nFrags.
Toggles: All/Tissue barcodes × Filtered/chrM/Blacklist metrics.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from collections import defaultdict
from typing import Dict, Tuple, Optional, List, Any

import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", required=True, help="binned fragment TSV or TSV.GZ")
    parser.add_argument("--mask", required=True, help="tissue mask TIF file")
    parser.add_argument("--mask-bin", type=int, default=20, help="mask bin size (pixels)")
    parser.add_argument("--mask-summary", default="", help="summary.json (auto if empty)")
    parser.add_argument("--bin-size", type=int, required=True, help="fragment bin size (e.g., 100)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (.csv and .html)")
    parser.add_argument("--exclude-chr", default="chrM", help="Chromosome to exclude (default: chrM, empty string = none)")
    parser.add_argument("--blacklist", default="", help="Optional ENCODE blacklist BED")
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_mask(path: str) -> np.ndarray:
    try:
        import tifffile
        arr = tifffile.imread(path)
    except Exception:
        from PIL import Image
        arr = np.array(Image.open(path))
    if arr.ndim != 2:
        raise ValueError(f"mask must be 2D, got shape={arr.shape}")
    return (arr > 0).astype(np.uint8)


def load_origin(summary_path: str, mask_bin: int) -> Tuple[int, int]:
    with open(summary_path, "r", encoding="utf-8") as f:
        summary = json.load(f)
    bin_size = int(summary["bin_size"])
    if bin_size != mask_bin:
        raise ValueError(f"mask-bin mismatch: summary={bin_size}, requested={mask_bin}")
    return int(summary["grid_origin_minx"]), int(summary["grid_origin_miny"])


def load_blacklist(path: str) -> Dict[str, np.ndarray]:
    intervals: Dict[str, list] = defaultdict(list)
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            intervals[parts[0]].append((start, end))

    result = {}
    for chrom, items in intervals.items():
        items.sort()
        arr = np.array(items, dtype=np.int64)
        cummax = np.maximum.accumulate(arr[:, 1])
        result[chrom] = np.column_stack([arr, cummax])
    return result


def in_blacklist(chrom: str, start: int, end: int, blacklist: Dict[str, np.ndarray]) -> bool:
    arr = blacklist.get(chrom)
    if arr is None:
        return False
    starts = arr[:, 0]
    cummax_ends = arr[:, 2]
    idx = np.searchsorted(starts, end, side="left")
    if idx <= 0:
        return False
    return bool(cummax_ends[idx - 1] > start)


def parse_cb(cb: str) -> Optional[Tuple[int, int]]:
    try:
        x_str, y_str = cb.split("_", 1)
        return int(x_str), int(y_str)
    except ValueError:
        return None


def cb_in_mask(x: int, y: int, mask: np.ndarray, mask_bin: int, minx: int, miny: int) -> bool:
    ax = (x // mask_bin) * mask_bin
    ay = (y // mask_bin) * mask_bin
    gx = (ax - minx) // mask_bin
    gy = (ay - miny) // mask_bin
    if gx < 0 or gy < 0 or gx >= mask.shape[1] or gy >= mask.shape[0]:
        return False
    return bool(mask[gy, gx] > 0)


def custom_colorscale() -> List[List]:
    colors = [
        "#0E458F", "#0F5298", "#0E6BA8", "#0C86B8", "#3399A1", "#3B9C9C",
        "#B2C061", "#F2CE38", "#F2AB38", "#F2AB38", "#EB7232", "#E65B2E",
        "#E14428", "#DC2E22", "#DB2921", "#CC2623"
    ]
    n = len(colors)
    positions = np.linspace(0, 1, n)
    return [[pos, col] for pos, col in zip(positions, colors)]


def main() -> int:
    args = parse_args()

    if not args.mask_summary:
        args.mask_summary = os.path.join(os.path.dirname(args.mask), "summary.json")

    print("[1/5] Loading tissue mask and origin...", file=sys.stderr)
    mask = load_mask(args.mask)
    minx, miny = load_origin(args.mask_summary, args.mask_bin)

    blacklist = load_blacklist(args.blacklist) if args.blacklist else {}
    if blacklist:
        n_intervals = sum(arr.shape[0] for arr in blacklist.values())
        print(f"[INFO] Loaded {n_intervals} blacklist intervals", file=sys.stderr)

    exclude_chr = args.exclude_chr or ""

    print("[2/5] Counting fragments by category...", file=sys.stderr)
    filtered_counts: Dict[str, int] = defaultdict(int)
    chrm_counts: Dict[str, int] = defaultdict(int)
    black_counts: Dict[str, int] = defaultdict(int)
    raw_counts: Dict[str, int] = defaultdict(int)
    coords: Dict[str, Tuple[int, int]] = {}

    total_lines = 0
    with open_text(args.fragments) as f:
        for line in f:
            total_lines += 1
            if total_lines % 5_000_000 == 0:
                print(f"  ... processed {total_lines:,} lines", file=sys.stderr)
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            chrom = parts[0]
            try:
                f_start = int(parts[1])
                f_end = int(parts[2])
            except ValueError:
                continue

            cb = parts[3]
            raw_counts[cb] += 1

            is_exclude = (exclude_chr and chrom == exclude_chr)
            is_black = False
            if blacklist:
                is_black = in_blacklist(chrom, f_start, f_end, blacklist)

            if is_black:
                black_counts[cb] += 1
            if is_exclude:
                chrm_counts[cb] += 1
            if not is_exclude and not is_black:
                filtered_counts[cb] += 1

            if cb not in coords:
                xy = parse_cb(cb)
                coords[cb] = xy if xy is not None else (-1, -1)

    print(f"[3/5] Unique barcodes: {len(raw_counts):,}", file=sys.stderr)

    # Build records and tissue membership
    records = []
    x_vals, y_vals = [], []
    for cb in raw_counts.keys():
        x, y = coords[cb]
        if x == -1:
            continue
        in_tissue = cb_in_mask(x, y, mask, args.mask_bin, minx, miny)
        records.append({
            "cb": cb,
            "x": x,
            "y": y,
            "in_tissue": in_tissue,
            "filtered": filtered_counts.get(cb, 0),
            "chrm": chrm_counts.get(cb, 0),
            "blacklist": black_counts.get(cb, 0),
            "raw": raw_counts[cb],
        })
        x_vals.append(x)
        y_vals.append(y)

    if not records:
        print("Error: No valid barcodes found!", file=sys.stderr)
        return 1

    # Determine grid dimensions
    x_unique = sorted(set(x_vals))
    y_unique = sorted(set(y_vals))
    x_to_idx = {x: i for i, x in enumerate(x_unique)}
    y_to_idx = {y: i for i, y in enumerate(y_unique)}
    shape = (len(y_unique), len(x_unique))

    # Matrices for raw values (0 or NaN for tissue view out-of-tissue)
    matrices_raw = {
        "all_filtered": np.zeros(shape, dtype=np.float32),
        "all_chrm": np.zeros(shape, dtype=np.float32),
        "all_blacklist": np.zeros(shape, dtype=np.float32),
        "tissue_filtered": np.full(shape, np.nan, dtype=np.float32),
        "tissue_chrm": np.full(shape, np.nan, dtype=np.float32),
        "tissue_blacklist": np.full(shape, np.nan, dtype=np.float32),
    }

    for rec in records:
        xi = x_to_idx[rec["x"]]
        yi = y_to_idx[rec["y"]]
        matrices_raw["all_filtered"][yi, xi] = rec["filtered"]
        matrices_raw["all_chrm"][yi, xi] = rec["chrm"]
        matrices_raw["all_blacklist"][yi, xi] = rec["blacklist"]
        if rec["in_tissue"]:
            matrices_raw["tissue_filtered"][yi, xi] = rec["filtered"]
            matrices_raw["tissue_chrm"][yi, xi] = rec["chrm"]
            matrices_raw["tissue_blacklist"][yi, xi] = rec["blacklist"]

    # Write CSV
    csv_path = args.output_prefix + ".csv"
    print(f"[4/5] Writing CSV to {csv_path}...", file=sys.stderr)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["CB", "x", "y", "in_tissue", "nFrags_filtered", "nFrags_chrM", "nFrags_blacklist", "nFrags_raw"])
        for rec in records:
            writer.writerow([
                rec["cb"], rec["x"], rec["y"], rec["in_tissue"],
                rec["filtered"], rec["chrm"], rec["blacklist"], rec["raw"]
            ])

    # Compute log10 matrices and max values
    log_matrices = {}
    for key, mat in matrices_raw.items():
        with np.errstate(divide='ignore', invalid='ignore'):
            log_mat = np.log10(mat + 1)
            log_mat[~np.isfinite(log_mat)] = np.nan
        log_matrices[key] = log_mat

    # Precompute max for each log matrix (for colorbar range)
    log_max_vals = {}
    for key, mat in log_matrices.items():
        finite_max = np.nanmax(mat)
        log_max_vals[key] = finite_max if np.isfinite(finite_max) else 0

    # Create initial trace (All / Filtered)
    initial_key = "all_filtered"
    initial_z = log_matrices[initial_key]
    initial_zmax = log_max_vals[initial_key]
    initial_customdata = matrices_raw[initial_key]
    initial_title = (f"All barcodes – Filtered nFrags (log10 scale)<br>"
                     f"excl. {args.exclude_chr or 'none'} & blacklist | bin size = {args.bin_size} bp")

    colorscale = custom_colorscale()
    trace = go.Heatmap(
        z=initial_z,
        x=x_unique,
        y=y_unique,
        colorscale=colorscale,
        zmin=0,
        zmax=initial_zmax,
        colorbar=dict(title="log10(nFrags+1)", len=0.6),
        hoverongaps=False,
        customdata=initial_customdata,
        hovertemplate="Barcode: %{x}_%{y}<br>nFrags: %{customdata:,.0f}<extra></extra>",
    )

    fig = go.Figure(data=[trace])

    # Prepare 6 buttons: each updates z, zmax, customdata, title
    buttons = []

    # Define all combinations: (region, metric) -> key_suffix
    combos = [
        ("All", "Filtered", "all_filtered"),
        ("All", "chrM", "all_chrm"),
        ("All", "Blacklist", "all_blacklist"),
        ("Tissue", "Filtered", "tissue_filtered"),
        ("Tissue", "chrM", "tissue_chrm"),
        ("Tissue", "Blacklist", "tissue_blacklist"),
    ]

    for region_label, metric_label, key in combos:
        button = dict(
            label=f"{region_label} / {metric_label}",
            method="update",
            args=[
                {
                    "z": [log_matrices[key]],
                    "zmax": log_max_vals[key],
                    "customdata": [matrices_raw[key]],
                    "colorbar.title": "log10(nFrags+1)",
                },
                {
                    "title": (f"{region_label} barcodes – {metric_label} nFrags (log10 scale)<br>"
                              f"bin size = {args.bin_size} bp")
                }
            ]
        )
        buttons.append(button)

    updatemenu = dict(
        buttons=buttons,
        direction="down",
        showactive=True,
        x=0.02,
        y=0.98,
        xanchor="left",
        yanchor="top",
        bgcolor="white",
        bordercolor="black",
        font=dict(color="black", size=12),
    )

    fig.update_layout(
        title=initial_title,
        xaxis_title="x coordinate (bin centre)",
        yaxis_title="y coordinate (bin centre)",
        xaxis=dict(
            tickvals=[],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            tickvals=[],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            autorange="reversed",
        ),
        width=900,
        height=800,
        hovermode="closest",
        updatemenus=[updatemenu],
        plot_bgcolor="black",
        paper_bgcolor="black",
        font=dict(color="white"),
    )

    html_path = args.output_prefix + ".html"
    plot(fig, filename=html_path, auto_open=False)
    print(f"\nOutput files saved:", file=sys.stderr)
    print(f"  CSV: {csv_path}", file=sys.stderr)
    print(f"  HTML: {html_path}", file=sys.stderr)

    # Summary
    total_barcodes = len(records)
    tissue_barcodes = sum(1 for r in records if r["in_tissue"])
    median_filtered_all = np.median([r["filtered"] for r in records if r["filtered"] > 0]) if any(r["filtered"]>0 for r in records) else 0
    print(f"\nSummary:", file=sys.stderr)
    print(f"  total barcodes            : {total_barcodes:,}", file=sys.stderr)
    print(f"  tissue barcodes (mask)    : {tissue_barcodes:,}", file=sys.stderr)
    print(f"  median filtered nFrags (all barcodes with >0): {median_filtered_all:.0f}", file=sys.stderr)
    print(f"  excluded chromosome       : {exclude_chr or 'none'}", file=sys.stderr)
    print(f"  blacklist                 : {args.blacklist or 'none'}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())