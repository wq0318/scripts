#!/usr/bin/env python3
"""Compute fragment-size distribution + nucleosome fractions from a sorted
merged fragment.tsv.gz (output of merge_fragment2.0.py).

Bin-independent: fragment sizes are intrinsic to read pairs, so the computation
runs once on the merged bin1 fragment file (previously was duplicated inside
signac_analysis.R per bin).

Outputs:
  <prefix>_fragment_size.csv     -- size,count,percentage rows for sizes 1..max_size
  <prefix>_fragment_data.js      -- ECharts payload (var fragmentData = [...])
  <prefix>_nucleosome_stats.tsv  -- fraction_nfr / fraction_mono / fraction_multi (%)

Excludes the chrM (or any user-specified) chromosome to match the previous
signac behavior. Insert sizes are clamped to (0, max_size].
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fragments", required=True, help="merged fragment tsv(.gz)")
    p.add_argument("--prefix", required=True,
                   help="output path prefix, e.g. .../03.fragment/${SN}")
    p.add_argument("--exclude-chr", default="chrM",
                   help="chromosome to exclude (default chrM, empty = none)")
    p.add_argument("--max-size", type=int, default=750,
                   help="upper insert-size bound (inclusive, default 750)")
    p.add_argument("--nfr-cut", type=int, default=147,
                   help="size < this = NFR (default 147)")
    p.add_argument("--mono-cut", type=int, default=294,
                   help="NFR <= size < this = mono-nucleosome (default 294)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not os.path.exists(args.fragments):
        print(f"[ERROR] fragments not found: {args.fragments}", file=sys.stderr)
        return 2

    counts = [0] * (args.max_size + 1)  # index 0 unused
    exclude = args.exclude_chr or ""

    with open_text(args.fragments) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if exclude and parts[0] == exclude:
                continue
            try:
                size = int(parts[2]) - int(parts[1])
            except ValueError:
                continue
            if 0 < size <= args.max_size:
                counts[size] += 1

    total = sum(counts)
    csv_path = f"{args.prefix}_fragment_size.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("size,count,percentage\n")
        for size in range(1, args.max_size + 1):
            pct = (counts[size] / total * 100) if total else 0.0
            f.write(f"{size},{counts[size]},{pct:.4f}\n")

    js_payload = ",".join(
        f"{{size:{size}, percentage:{(counts[size]/total*100 if total else 0):.2f}}}"
        for size in range(1, args.max_size + 1)
    )
    js_path = f"{args.prefix}_fragment_data.js"
    with open(js_path, "w", encoding="utf-8") as f:
        f.write(f"var fragmentData = [{js_payload}];\n")

    if total > 0:
        nfr = sum(counts[1:args.nfr_cut])
        mono = sum(counts[args.nfr_cut:args.mono_cut])
        multi = sum(counts[args.mono_cut:args.max_size + 1])
        fraction_nfr = nfr / total * 100
        fraction_mono = mono / total * 100
        fraction_multi = multi / total * 100
    else:
        fraction_nfr = fraction_mono = fraction_multi = 0.0

    stats_path = f"{args.prefix}_nucleosome_stats.tsv"
    with open(stats_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"fraction_nfr\t{fraction_nfr:.2f}\n")
        f.write(f"fraction_mono_nucleosome\t{fraction_mono:.2f}\n")
        f.write(f"fraction_multi_nucleosome\t{fraction_multi:.2f}\n")
        f.write(f"total_insert_pairs\t{total}\n")

    print(f"[fragment_size] total={total:,} NFR={fraction_nfr:.2f}% mono={fraction_mono:.2f}% multi={fraction_multi:.2f}%",
          file=sys.stderr)
    print(f"[fragment_size] wrote {csv_path}", file=sys.stderr)
    print(f"[fragment_size] wrote {js_path}", file=sys.stderr)
    print(f"[fragment_size] wrote {stats_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
