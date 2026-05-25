#!/usr/bin/env python3
"""Adapt the 8-column grid CSV from generate_nfrags_grid.py into the legacy
2-column tissue_barcode.csv (`CB,nFrags`) + summary TSV consumed by
atac_saturation.py / signac_analysis.R / parse_downstream_analysis_outputs.py.

nFrags semantics match compute_tissue_barcodes.py: filtered count (excludes
exclude_chr and blacklist), tissue mask membership taken from in_tissue.
"""

from __future__ import annotations

import argparse
import csv
import sys
from statistics import median


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="8-column grid CSV from generate_nfrags_grid.py")
    p.add_argument("--output", required=True, help="compatibility tissue_barcode.csv (CB,nFrags)")
    p.add_argument("--summary", required=True, help="compatibility tissue_barcode_summary.tsv")
    p.add_argument("--bin-size", type=int, required=True)
    p.add_argument("--exclude-chr", default="chrM")
    p.add_argument("--blacklist", default="", help="blacklist path to record in summary; empty = 'none'")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    total_barcodes = 0
    valid_barcodes = 0
    excluded_rows = 0
    excluded_blacklist_rows = 0
    tissue_nfrags: list[int] = []

    with open(args.input, "r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin)
        required = {"CB", "in_tissue", "nFrags_filtered", "nFrags_chrM", "nFrags_blacklist"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            print(f"[ERROR] Missing columns in {args.input}: {sorted(missing)}", file=sys.stderr)
            return 2

        with open(args.output, "w", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout)
            writer.writerow(["CB", "nFrags"])
            rows: list[tuple[str, int]] = []
            for row in reader:
                total_barcodes += 1
                try:
                    nf = int(row["nFrags_filtered"])
                    nc = int(row["nFrags_chrM"])
                    nb = int(row["nFrags_blacklist"])
                except (TypeError, ValueError):
                    continue
                excluded_rows += nc
                excluded_blacklist_rows += nb
                if row["in_tissue"].strip().lower() == "true":
                    valid_barcodes += 1
                    tissue_nfrags.append(nf)
                    rows.append((row["CB"], nf))
            rows.sort(key=lambda r: r[0])
            for cb, nf in rows:
                writer.writerow([cb, nf])

    med = median(tissue_nfrags) if tissue_nfrags else 0
    blacklist_label = args.blacklist if args.blacklist else "none"
    with open(args.summary, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"bin_size\t{args.bin_size}\n")
        f.write(f"total_barcodes\t{total_barcodes}\n")
        f.write(f"valid_barcodes\t{valid_barcodes}\n")
        f.write(f"median_nFrags\t{med:.0f}\n")
        f.write(f"excluded_chr\t{args.exclude_chr}\n")
        f.write(f"excluded_rows\t{excluded_rows}\n")
        f.write(f"blacklist\t{blacklist_label}\n")
        f.write(f"excluded_blacklist_rows\t{excluded_blacklist_rows}\n")

    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
