#!/usr/bin/env python3
"""Render the StereoATAC HTML report via Jinja2 (`templates/report.html.j2`).

Replaces the previous f-string-heavy renderer. Parses per-lane QC CSVs +
per-bin downstream/nucleosome/saturation/fragment outputs, base64-embeds all
referenced SVG/PNG assets, and hands the assembled context dict to Jinja2.
"""

from __future__ import annotations

import argparse
import base64
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

from jinja2 import Environment, FileSystemLoader, select_autoescape


# ---------- formatting ----------

def parse_path_list(value: str) -> list[str]:
    return [s.strip() for s in value.split(",") if s.strip()]


def to_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def fmt_int(value, fallback: str = "0") -> str:
    f = to_float(value)
    return fallback if f is None else f"{int(round(f)):,}"


def fmt_pct(value, digits: int = 2, fallback: str = "0.00%") -> str:
    f = to_float(value)
    return fallback if f is None else f"{f:.{digits}f}%"


def fmt_num(value, digits: int = 2, fallback: str = "0") -> str:
    f = to_float(value)
    return fallback if f is None else f"{f:.{digits}f}"


# ---------- image embedding ----------

def mime_for(path: str) -> str:
    suf = Path(path).suffix.lower()
    if suf == ".svg":
        return "image/svg+xml"
    if suf in {".jpg", ".jpeg"}:
        return "image/jpeg"
    if suf == ".webp":
        return "image/webp"
    return "image/png"


def placeholder_uri(title: str) -> str:
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 400">'
        '<rect width="800" height="400" fill="#f7fafc" stroke="#e2e8f0"/>'
        f'<text x="400" y="200" text-anchor="middle" font-family="Arial" '
        f'font-size="18" fill="#718096">{title}</text></svg>'
    )
    return "data:image/svg+xml;charset=utf-8," + quote(svg)


def embed_image(path: str, title: str) -> str:
    if not path or not os.path.exists(path):
        return placeholder_uri(title)
    with open(path, "rb") as fh:
        return f"data:{mime_for(path)};base64,{base64.b64encode(fh.read()).decode('ascii')}"


# ---------- parsers ----------

def read_csv_rows(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for p in paths:
        if not p or not os.path.exists(p):
            continue
        with open(p, encoding="utf-8", errors="replace") as fh:
            rows.extend(csv.DictReader(fh))
    return rows


def read_metric_map(path: str, delimiter: str = ",") -> dict[str, str]:
    out: dict[str, str] = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter):
            if "metric" in row and "value" in row:
                out[row["metric"]] = row["value"]
    return out


def read_saturation(path: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, encoding="utf-8", errors="replace") as fh:
        fh.readline()  # header
        for line in fh:
            fields = line.strip().split()
            if len(fields) < 5:
                continue
            try:
                rows.append({
                    "percent": float(fields[0]),
                    "unique": float(fields[1]),
                    "total": float(fields[2]),
                    "dup_pct": float(fields[3]),
                    "median_nfrags": float(fields[4]),
                })
            except ValueError:
                continue
    return rows


def read_fragment_size(path: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            try:
                rows.append({
                    "size": int(float(row["size"])),
                    "percentage": round(float(row["percentage"]), 4),
                })
            except (KeyError, TypeError, ValueError):
                continue
    return rows


def read_cluster_first(value: str) -> str:
    if value and os.path.exists(value):
        with open(value, encoding="utf-8", errors="replace") as fh:
            paths = parse_path_list(fh.read())
            if paths:
                return paths[0]
    return ""


# ---------- aggregators ----------

def build_global_metrics(qc_rows: list[dict[str, str]]) -> dict[str, str]:
    sums: defaultdict[str, float] = defaultdict(float)
    latest: dict[str, str] = {}
    for row in qc_rows:
        m, v = row.get("metric", ""), row.get("value", "")
        latest[m] = v
        f = to_float(v)
        if f is not None:
            sums[m] += f

    total_reads = sums.get("total_reads", 0)
    mapped = sums.get("mapped_reads", 0)
    exact = sums.get("barcode_exactly_overlap_reads", 0)
    mis = sums.get("barcode_mis_overlap_reads", 0)
    chromap_total = sums.get("chromap_total_reads", 0)
    chromap_hq = sums.get("chromap_high_quality_reads", 0)

    def pct(num: float, den: float) -> str:
        return f"{(num * 100 / den):.2f}%" if den else "0.00%"

    return {
        "total_reads": fmt_int(total_reads),
        "barcode_mapping_rate": pct(mapped, total_reads),
        "exact_overlap_rate": pct(exact, total_reads),
        "mismatch_overlap_rate": pct(mis, total_reads),
        "genome_mapping_rate": pct(chromap_total, mapped),
        "high_quality_mapping": pct(chromap_hq, chromap_total),
        "low_quality_mapping": pct(chromap_total - chromap_hq, chromap_total),
        "sequencing_saturation": fmt_pct(latest.get("saturation", "0")),
        "mitochondria_ratio": fmt_pct(latest.get("mitochondria_ratio", "0")),
    }


def pad_list(items: list[str], n: int, label: str) -> list[str]:
    items = list(items)
    if len(items) < n:
        print(f"[WARN] {label}: expected {n}, got {len(items)}; padding", file=sys.stderr)
        items += [""] * (n - len(items))
    elif len(items) > n:
        print(f"[WARN] {label}: expected {n}, got {len(items)}; truncating", file=sys.stderr)
        items = items[:n]
    return items


def resolve_template_dir(explicit: str) -> Path:
    """Locate templates/report.html.j2: --template-dir wins, then alongside
    the script, then sibling to the scripts/ folder. The fallback chain lets
    the same renderer work locally and inside docker without recompiling the
    image when the layout differs."""
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    here = Path(__file__).resolve().parent
    candidates.append(here / "templates")
    candidates.append(here.parent / "templates")
    for c in candidates:
        if (c / "report.html.j2").exists():
            return c
    raise FileNotFoundError(
        f"report.html.j2 not found in any of: {[str(c) for c in candidates]}"
    )


# ---------- main ----------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render StereoATAC final HTML report via Jinja2.")
    p.add_argument("sample_name")
    p.add_argument("output_html")
    p.add_argument("barcode_metrics")
    p.add_argument("mapping_metrics")
    p.add_argument("sam2frag_metrics")
    p.add_argument("merge_metrics")
    p.add_argument("--bin-sizes", required=True)
    p.add_argument("--downstream-metrics", required=True)
    p.add_argument("--nucleosome-stats", required=True)
    p.add_argument("--saturation-tables", required=True)
    p.add_argument("--fragment-stats-plots", default="",
                   help="Accepted for WDL compat; not consumed by the Jinja2 renderer")
    p.add_argument("--qc-violins", required=True)
    p.add_argument("--tss-scatters", required=True)
    p.add_argument("--fragment-sizes", required=True)
    p.add_argument("--spatial-qcs", default="",
                   help="Accepted for WDL compat; superseded by --spatial-tiles-bin100")
    p.add_argument("--spatial-qcs-pre", default="",
                   help="Accepted for WDL compat; superseded by --spatial-tiles-bin100")
    p.add_argument("--cluster-plots-list", required=True)
    p.add_argument("--spatial-tiles-bin100", default="",
                   help="Comma-separated 5 SVGs from plot_spatial_tiles.R "
                        "(order: filtered, raw, chrM, blacklist, tissue_filtered)")
    p.add_argument("--spatial-tss-bin100", default="",
                   help="bin100 _bins_under_tissue_TSS.svg from signac_analysis.R")
    p.add_argument("--template-dir", default="",
                   help="Override path to dir containing report.html.j2")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    bin_sizes = [s.strip() for s in args.bin_sizes.split(",") if s.strip()]
    n = len(bin_sizes)

    qc_rows = read_csv_rows([
        *parse_path_list(args.barcode_metrics),
        *parse_path_list(args.mapping_metrics),
        *parse_path_list(args.sam2frag_metrics),
        args.merge_metrics,
    ])
    metrics = build_global_metrics(qc_rows)

    downstream_list = pad_list(parse_path_list(args.downstream_metrics), n, "downstream-metrics")
    nucleosome_list = pad_list(parse_path_list(args.nucleosome_stats), n, "nucleosome-stats")
    saturation_list = pad_list(parse_path_list(args.saturation_tables), n, "saturation-tables")
    qc_violins = pad_list(parse_path_list(args.qc_violins), n, "qc-violins")
    tss_scatters = pad_list(parse_path_list(args.tss_scatters), n, "tss-scatters")
    fragment_sizes = pad_list(parse_path_list(args.fragment_sizes), n, "fragment-sizes")
    cluster_lists = pad_list(
        [s.strip() for s in args.cluster_plots_list.split(";")], n, "cluster-plots-list"
    )

    bins_ctx: list[dict] = []
    primary_ns: dict[str, str] = {}
    primary_blacklist = "0"
    saturation_payload = {"duplicate_rate": [], "median_nfrags": []}

    for i, b in enumerate(bin_sizes):
        ds = read_metric_map(downstream_list[i], delimiter=",")
        ns = read_metric_map(nucleosome_list[i], delimiter="\t")
        if i == 0:
            primary_ns = ns
            primary_blacklist = ds.get("median_blacklist", "0")

        sat_rows = read_saturation(saturation_list[i])
        saturation_payload["duplicate_rate"].append({
            "bin": b,
            "points": [[round(r["total"] / 1_000_000.0, 3), round(r["dup_pct"], 2)] for r in sat_rows],
        })
        saturation_payload["median_nfrags"].append({
            "bin": b,
            "points": [[round(r["total"] / 1_000_000.0, 3), round(r["median_nfrags"], 1)] for r in sat_rows],
        })

        bins_ctx.append({
            "size": b,
            "total_spots": fmt_int(ds.get("total_barcodes", "0")),
            "valid_spots": fmt_int(ds.get("valid_barcodes", "0")),
            "median_nFrags": fmt_int(ds.get("median_nFrags", "0")),
            "median_TSS": fmt_num(ds.get("median_TSS", "0")),
            "median_FRiP": fmt_num(ds.get("median_FRiP", "0")),
            "qc_violin": embed_image(qc_violins[i], f"bin{b} QC violin"),
            "tss_scatter": embed_image(tss_scatters[i], f"bin{b} TSS scatter"),
            "cluster": embed_image(read_cluster_first(cluster_lists[i]),
                                   f"bin{b} cluster res=0.8"),
        })

    metrics["blacklist_ratio"] = fmt_pct(primary_blacklist)
    metrics["fraction_nfr"] = fmt_pct(primary_ns.get("fraction_nfr", "0"))
    metrics["fraction_mono"] = fmt_pct(primary_ns.get("fraction_mono_nucleosome", "0"))
    metrics["fraction_multi"] = fmt_pct(primary_ns.get("fraction_multi_nucleosome", "0"))

    tile_paths = pad_list(parse_path_list(args.spatial_tiles_bin100), 5, "spatial-tiles-bin100")
    # Order MUST match plot_spatial_tiles.R output:
    #   _all_bins_nFrags_filtered, _all_bins_nFrags_raw, _all_bins_nFrags_chrM,
    #   _all_bins_nFrags_blacklist, _bins_under_tissue_nFrags_filtered
    tile_keys = ["all_filtered", "all_raw", "all_chrm", "all_blacklist", "tissue_filtered"]
    tile_titles = {
        "all_filtered": "All bins - nFrags_filtered",
        "all_raw": "All bins - nFrags_raw",
        "all_chrm": "All bins - nFrags_chrM",
        "all_blacklist": "All bins - nFrags_blacklist",
        "tissue_filtered": "Bins under tissue - nFrags_filtered",
    }
    tiles_ctx = {k: embed_image(p, tile_titles[k]) for k, p in zip(tile_keys, tile_paths)}
    tiles_ctx["tissue_tss"] = embed_image(args.spatial_tss_bin100, "Bins under tissue - TSS")

    primary_fragment = read_fragment_size(fragment_sizes[0]) if fragment_sizes else []

    context = {
        "sample_name": args.sample_name,
        "bins_label": ", ".join(f"bin{b}" for b in bin_sizes),
        "metrics": metrics,
        "tiles": tiles_ctx,
        "bins": bins_ctx,
        "saturation_data_json": json.dumps(saturation_payload),
        "fragment_data_json": json.dumps(primary_fragment),
    }

    template_dir = resolve_template_dir(args.template_dir)
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    rendered = env.get_template("report.html.j2").render(**context)

    out_path = Path(args.output_html)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
