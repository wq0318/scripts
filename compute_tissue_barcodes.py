#!/usr/bin/env python3
"""Create tissue barcode tables for binned StereoATAC fragments.

The input fragment file is kept intact for downstream analysis.  This script
only uses the tissue mask to decide which bin barcodes are inside tissue, then
reports each barcode with its unique-fragment count.

nFrags semantics
----------------
``nFrags`` here is the *unique-fragment* count per barcode after dropping
``--exclude-chr`` rows and ``--blacklist`` overlaps. Each row of the input
fragment file represents one unique fragment; column 5 (the duplicate-read
count from PISA) is **not** used as a multiplier. This matches the convention
used by ``atac_saturation.py`` and Cell Ranger ATAC's per-cell
``unique_fragments``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
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
    parser.add_argument(
        "--exclude-chr",
        default="chrM",
        help=(
            "Chromosome name to exclude from nFrags accumulation (default: chrM). "
            "Mitochondrial reads are PCR-prone and inflate duplicate rate; every other "
            "chromosome (including chrY and unplaced scaffolds) is real signal and kept. "
            "Pass an empty string to disable the filter entirely."
        ),
    )
    parser.add_argument(
        "--blacklist",
        default="",
        help=(
            "Optional ENCODE-style blacklist BED (chr/start/end at columns 1-3). "
            "Fragments whose [start, end) overlaps any blacklist interval are dropped "
            "before nFrags accumulation. Pass an empty string (default) to disable."
        ),
    )
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


def load_blacklist(path: str) -> dict[str, np.ndarray]:
    """Parse a blacklist BED into ``{chrom: ndarray[(N, 3)]}``.

    Columns 0/1 of each row are ``(start, end)``; column 2 is the running
    ``max(ends[0..i])`` (a prefix-max). With that, ``in_blacklist`` can answer
    "is there any interval i with starts[i] < frag.end AND ends[i] > frag.start"
    via one binary search plus one comparison.
    """
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("track") or stripped.startswith("browser"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if end <= start:
                continue
            intervals[fields[0]].append((start, end))
    result: dict[str, np.ndarray] = {}
    for chrom, items in intervals.items():
        items.sort()
        arr = np.array(items, dtype=np.int64)
        cummax = np.maximum.accumulate(arr[:, 1])
        result[chrom] = np.column_stack([arr, cummax])
    return result


def in_blacklist(chrom: str, start: int, end: int, blacklist: dict[str, np.ndarray]) -> bool:
    arr = blacklist.get(chrom)
    if arr is None:
        return False
    starts = arr[:, 0]
    cummax_ends = arr[:, 2]
    idx = int(np.searchsorted(starts, end, side="left"))
    if idx <= 0:
        return False
    return bool(cummax_ends[idx - 1] > start)


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

    exclude_chr = args.exclude_chr or ""
    mask = load_mask(args.mask)
    minx, miny = load_origin(args.mask_summary, args.mask_bin)
    blacklist = load_blacklist(args.blacklist) if args.blacklist else {}
    if blacklist:
        n_intervals = sum(arr.shape[0] for arr in blacklist.values())
        print(
            f"[INFO] Loaded {n_intervals} blacklist intervals across "
            f"{len(blacklist)} chromosomes from {args.blacklist}",
            file=sys.stderr,
        )

    all_counts: dict[str, int] = defaultdict(int)
    tissue_cache: dict[str, bool] = {}
    tissue_counts: dict[str, int] = defaultdict(int)
    excluded_rows = 0
    excluded_blacklist_rows = 0

    with open_text(args.fragments) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            if exclude_chr and fields[0] == exclude_chr:
                excluded_rows += 1
                continue
            if blacklist:
                try:
                    f_start = int(fields[1])
                    f_end = int(fields[2])
                except ValueError:
                    continue
                if in_blacklist(fields[0], f_start, f_end, blacklist):
                    excluded_blacklist_rows += 1
                    continue
            cb = fields[3]
            # Each fragment row counts as 1 unique fragment. Column 5 (PISA's
            # duplicate-read count) is intentionally ignored so nFrags matches
            # the unique-fragment semantics used by atac_saturation.py and
            # Cell Ranger ATAC.
            all_counts[cb] += 1
            keep = tissue_cache.get(cb)
            if keep is None:
                keep = cb_in_mask(cb, mask, args.mask_bin, minx, miny)
                tissue_cache[cb] = keep
            if keep:
                tissue_counts[cb] += 1

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
        handle.write(f"excluded_chr\t{exclude_chr}\n")
        handle.write(f"excluded_rows\t{excluded_rows}\n")
        handle.write(f"blacklist\t{args.blacklist or 'none'}\n")
        handle.write(f"excluded_blacklist_rows\t{excluded_blacklist_rows}\n")

    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
