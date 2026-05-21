#!/usr/bin/env python3
"""Create tissue barcode tables for binned StereoATAC fragments.

The input fragment file is kept intact for downstream analysis.  This script
only uses the tissue mask to decide which bin barcodes are inside tissue, then
reports each barcode with its fragment count.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from collections import defaultdict
from statistics import median

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute bin tissue barcodes and nFrags from a tissue mask."
    )
    parser.add_argument("--fragments", required=True, help="Binned fragment TSV or TSV.GZ")
    parser.add_argument("--mask", required=True, help="Mask tif from stereoatac_tissuecut_v1.1.py")
    parser.add_argument("--mask-bin", type=int, default=20, help="Mask bin size")
    parser.add_argument(
        "--mask-summary",
        default="",
        help="summary.json next to mask.tif. Defaults to <mask_dir>/summary.json",
    )
    parser.add_argument("--bin-size", type=int, required=True, help="Bin size of fragment CB")
    parser.add_argument("--output", required=True, help="Output CSV with CB,nFrags columns")
    parser.add_argument("--summary", required=True, help="Output metric/value TSV summary")
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_mask(path: str) -> np.ndarray:
    try:
        import tifffile  # type: ignore

        arr = tifffile.imread(path)
    except Exception:
        from PIL import Image  # type: ignore

        arr = np.array(Image.open(path))
    if arr.ndim != 2:
        raise ValueError(f"mask must be 2D, got shape={arr.shape}")
    return (arr > 0).astype(np.uint8)


def load_origin(summary_path: str, mask_bin: int) -> tuple[int, int]:
    with open(summary_path, "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    bin_size = int(summary["bin_size"])
    if bin_size != int(mask_bin):
        raise ValueError(f"mask-bin mismatch: summary={bin_size}, requested={mask_bin}")
    return int(summary["grid_origin_minx"]), int(summary["grid_origin_miny"])


def parse_cb(cb: str) -> tuple[int, int] | None:
    try:
        x_str, y_str = cb.split("_", 1)
        return int(x_str), int(y_str)
    except ValueError:
        return None


def cb_in_mask(cb: str, mask: np.ndarray, mask_bin: int, minx: int, miny: int) -> bool:
    parsed = parse_cb(cb)
    if parsed is None:
        return False
    x, y = parsed
    ax = (x // mask_bin) * mask_bin
    ay = (y // mask_bin) * mask_bin
    gx = (ax - minx) // mask_bin
    gy = (ay - miny) // mask_bin
    if gx < 0 or gy < 0 or gx >= mask.shape[1] or gy >= mask.shape[0]:
        return False
    return bool(mask[gy, gx] > 0)


def main() -> int:
    args = parse_args()
    if not args.mask_summary:
        args.mask_summary = os.path.join(os.path.dirname(args.mask), "summary.json")

    mask = load_mask(args.mask)
    minx, miny = load_origin(args.mask_summary, args.mask_bin)

    all_counts: dict[str, int] = defaultdict(int)
    tissue_cache: dict[str, bool] = {}
    tissue_counts: dict[str, int] = defaultdict(int)

    with open_text(args.fragments) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            cb = fields[3]
            try:
                count = int(float(fields[4]))
            except ValueError:
                count = 1
            all_counts[cb] += count
            keep = tissue_cache.get(cb)
            if keep is None:
                keep = cb_in_mask(cb, mask, args.mask_bin, minx, miny)
                tissue_cache[cb] = keep
            if keep:
                tissue_counts[cb] += count

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["CB", "nFrags"])
        for cb, count in sorted(tissue_counts.items()):
            writer.writerow([cb, count])

    values = list(tissue_counts.values())
    med = median(values) if values else 0
    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"bin_size\t{args.bin_size}\n")
        handle.write(f"total_barcodes\t{len(all_counts)}\n")
        handle.write(f"valid_barcodes\t{len(tissue_counts)}\n")
        handle.write(f"median_nFrags\t{med:.0f}\n")

    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
