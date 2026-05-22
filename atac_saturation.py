#!/usr/bin/env python3
"""Analytical ATAC-seq saturation curve (方案 D / closed-form expectation).

Replaces ``atac_saturation.pl``. For each fragment row with duplicate read
count ``num``, the probability that the fragment survives a sub-sampling at
percent ``p`` equals ``1 - (1 - p)^num``. Summing this over rows gives the
expected number of unique fragments at ``p`` — no Monte-Carlo random
sampling, no in-memory expansion of every read.

Input  : fragment.tsv(.gz) with 5 columns (chr, start, end, barcode, num).
Output : ``{output_prefix}_result.txt`` with header

    Percent  Unique_Fragments  Total_Fragments  Percent_Duplicates  Median_nFrags

at 8 sampling points (12.5, 25, ..., 100). Format identical to
``atac_saturation.pl`` so downstream parsers don't need to change.

Tissue-aware mode
-----------------
When ``--barcode-whitelist`` is given (e.g. the ``CB`` column of
``compute_tissue_barcodes.py``'s ``_tissue_barcode.csv``), fragments whose
barcode is **not** in the whitelist are dropped before any statistic is
computed. This matches the Cell Ranger ATAC Library Complexity convention,
where the per-cell saturation curve is computed over called cells only —
the spatial analogue being tissue-internal bins only. Both Y axes of the
HTML report's Section 2 charts (Duplicate Rate and Median nFrags / Bin)
then reflect the in-tissue signal.

Notes on semantics vs the legacy Perl implementation
----------------------------------------------------
* ``Unique_Fragments`` here is the *expected* count under the closed-form
  model. The Perl script reports a single random realisation (Fisher-Yates
  shuffle). Both are statistically valid; the analytical form has zero
  Monte-Carlo variance and is what Cell Ranger ATAC / ArchR use.
* ``Median_nFrags`` is computed over the (whitelisted, if any) barcodes
  with at least one fragment in the file, using each barcode's expected
  unique-fragment count at percent ``p``.
* Only the row's ``chr`` column matching ``--exclude-chr`` (default
  ``chrM``) is skipped. Mitochondrial reads are PCR-prone and would
  inflate the duplicate rate; every other chromosome is real signal and
  counted as-is. The script exits with code 2 if no rows survive (e.g.
  the fragment file is empty after the filter).
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime

import numpy as np
import pandas as pd

PERCENTS = np.array([12.5, 25.0, 37.5, 50.0, 62.5, 75.0, 87.5, 100.0])
P_FRAC = PERCENTS / 100.0


def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analytical ATAC saturation curve (drop-in for atac_saturation.pl)",
    )
    parser.add_argument("-i", "--input", required=True, help="fragment.tsv(.gz)")
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="output prefix; writes {prefix}_result.txt",
    )
    parser.add_argument(
        "-b",
        "--barcode-whitelist",
        default=None,
        help=(
            "optional CSV whose first column is a barcode (e.g. compute_tissue_barcodes.py's "
            "_tissue_barcode.csv with header CB,nFrags). When supplied, only rows whose barcode "
            "is in this set contribute to the saturation curve — the tissue-aware analogue of "
            "Cell Ranger ATAC's per-called-cell library complexity."
        ),
    )
    parser.add_argument(
        "--exclude-chr",
        default="chrM",
        help=(
            "Chromosome name to drop before counting (default: chrM). Mitochondrial "
            "reads are PCR-prone and would inflate the duplicate rate; every other "
            "chromosome (chrY, unplaced scaffolds, etc.) is real signal and kept. "
            "Pass an empty string to disable."
        ),
    )
    parser.add_argument(
        "--blacklist",
        default="",
        help=(
            "Optional ENCODE-style blacklist BED (chr/start/end at columns 1-3). "
            "Fragments whose [start, end) overlaps any blacklist interval are dropped "
            "before saturation accumulation. Pass an empty string (default) to disable."
        ),
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=2_000_000,
        help="rows per pandas chunk (default 2M, ~230MB peak per chunk)",
    )
    return parser.parse_args()


def load_whitelist(path: str) -> set[str]:
    """Read the first column of a CSV (with or without header) as a barcode set."""
    df = pd.read_csv(path, usecols=[0], dtype=str)
    col = df.columns[0]
    barcodes = df[col].dropna().astype(str).str.strip()
    barcodes = barcodes[barcodes != ""]
    barcodes = barcodes[barcodes.str.lower() != "cb"]  # drop header if read positionally
    return set(barcodes)


def load_blacklist(path: str) -> dict[str, np.ndarray]:
    """Parse a blacklist BED into ``{chrom: ndarray[N, 3]}``.

    Columns 0/1 are sorted ``(start, end)``; column 2 is the prefix-max of the
    ``end`` column. The prefix-max lets the chunk-level filter answer "does any
    interval i with starts[i] < frag.end have ends[i] > frag.start?" with a
    single np.searchsorted plus one comparison — fully vectorised per chrom.
    """
    intervals: dict[str, list[tuple[int, int]]] = {}
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
            intervals.setdefault(fields[0], []).append((start, end))
    result: dict[str, np.ndarray] = {}
    for chrom, items in intervals.items():
        items.sort()
        arr = np.array(items, dtype=np.int64)
        cummax = np.maximum.accumulate(arr[:, 1])
        result[chrom] = np.column_stack([arr, cummax])
    return result


def blacklist_drop_mask(
    chr_arr: np.ndarray,
    start_arr: np.ndarray,
    end_arr: np.ndarray,
    blacklist: dict[str, np.ndarray],
) -> np.ndarray:
    """Return a boolean ndarray: ``True`` where the row overlaps blacklist."""
    drop = np.zeros(len(chr_arr), dtype=bool)
    if not blacklist:
        return drop
    for chrom, arr in blacklist.items():
        sel = chr_arr == chrom
        if not sel.any():
            continue
        starts = arr[:, 0]
        cummax_ends = arr[:, 2]
        sub_starts = start_arr[sel]
        sub_ends = end_arr[sel]
        idx = np.searchsorted(starts, sub_ends, side="left")
        sub_drop = np.zeros(sub_starts.shape[0], dtype=bool)
        valid = idx > 0
        if valid.any():
            sub_drop[valid] = cummax_ends[idx[valid] - 1] > sub_starts[valid]
        drop[sel] = sub_drop
    return drop


def main() -> int:
    args = parse_args()
    exclude_chr = args.exclude_chr or ""

    whitelist: set[str] | None = None
    if args.barcode_whitelist:
        whitelist = load_whitelist(args.barcode_whitelist)
        if not whitelist:
            log(f"[ERROR] barcode whitelist {args.barcode_whitelist} is empty after parsing.")
            return 2
        log(f"Loaded {len(whitelist):,} whitelisted barcodes from {args.barcode_whitelist}")

    blacklist: dict[str, np.ndarray] = {}
    if args.blacklist:
        blacklist = load_blacklist(args.blacklist)
        n_intervals = sum(arr.shape[0] for arr in blacklist.values())
        log(
            f"Loaded {n_intervals:,} blacklist intervals across {len(blacklist)} "
            f"chromosomes from {args.blacklist}"
        )

    out_path = f"{args.output_prefix}_result.txt"

    expected_unique = np.zeros_like(P_FRAC)
    total_reads = 0
    total_rows = 0
    blacklist_dropped = 0
    bc_acc: dict[str, np.ndarray] = {}

    log(f"Reading {args.input} (chunksize={args.chunksize:,})")

    reader = pd.read_csv(
        args.input,
        sep="\t",
        header=None,
        comment="#",
        names=["chr", "start", "end", "bc", "num"],
        dtype={"chr": str, "start": np.int64, "end": np.int64, "bc": str, "num": np.int64},
        chunksize=args.chunksize,
        engine="c",
        compression="infer",
    )

    for chunk_idx, chunk in enumerate(reader):
        if exclude_chr:
            chunk = chunk.loc[chunk["chr"] != exclude_chr]
        if blacklist and not chunk.empty:
            drop_mask = blacklist_drop_mask(
                chunk["chr"].to_numpy(),
                chunk["start"].to_numpy(),
                chunk["end"].to_numpy(),
                blacklist,
            )
            blacklist_dropped += int(drop_mask.sum())
            chunk = chunk.loc[~drop_mask]
        if chunk.empty:
            continue
        chunk = chunk[["bc", "num"]]
        if whitelist is not None and not chunk.empty:
            chunk = chunk.loc[chunk["bc"].isin(whitelist)]
        if chunk.empty:
            continue

        num = chunk["num"].to_numpy(dtype=np.float64)
        survival = 1.0 - np.power(1.0 - P_FRAC[None, :], num[:, None])  # (N, 8)

        total_reads += int(num.sum())
        total_rows += len(chunk)
        expected_unique += survival.sum(axis=0)

        # Per-barcode aggregation: groupby once per chunk and merge into global dict.
        df = pd.DataFrame(survival, columns=[f"p{i}" for i in range(len(P_FRAC))])
        df["bc"] = chunk["bc"].to_numpy()
        grouped = df.groupby("bc", sort=False).sum()
        bcs = grouped.index.to_numpy()
        sums = grouped.to_numpy(dtype=np.float64)
        for i, bc in enumerate(bcs):
            if bc in bc_acc:
                bc_acc[bc] += sums[i]
            else:
                bc_acc[bc] = sums[i].copy()

        if (chunk_idx + 1) % 10 == 0:
            log(
                f"  processed {total_rows:,} rows | total_reads={total_reads:,} | "
                f"barcodes={len(bc_acc):,}"
            )

    log(
        f"Scan complete: {total_rows:,} unique fragments, {total_reads:,} reads, "
        f"{len(bc_acc):,} barcodes (blacklist dropped: {blacklist_dropped:,})"
    )

    if total_rows == 0:
        if whitelist is not None:
            log(
                "[ERROR] No fragments left after applying --exclude-chr and --barcode-whitelist; "
                "nothing to write. Check whether the tissue mask is empty or the fragment "
                "file contains only the excluded chromosome."
            )
        else:
            log(
                "[ERROR] No fragments left after applying --exclude-chr; nothing to write. "
                "Check input or pass --exclude-chr '' to disable the filter."
            )
        return 2

    if bc_acc:
        stacked = np.vstack(list(bc_acc.values()))  # (B, 8)
        medians = np.median(stacked, axis=0)
    else:
        medians = np.zeros_like(P_FRAC)

    log(f"Writing {out_path}")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Percent\tUnique_Fragments\tTotal_Fragments\tPercent_Duplicates\tMedian_nFrags\n")
        for i, p in enumerate(PERCENTS):
            uniq = expected_unique[i]
            total_at_p = total_reads * P_FRAC[i]
            dup_pct = (1.0 - uniq / total_at_p) * 100.0 if total_at_p > 0 else 0.0
            fh.write(
                f"{p:.1f}\t{int(round(uniq))}\t{int(round(total_at_p))}\t"
                f"{dup_pct:.2f}\t{medians[i]:.2f}\n"
            )

    log("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
