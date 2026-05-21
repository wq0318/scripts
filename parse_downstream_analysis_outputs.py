#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import sys
from datetime import datetime


def read_key_value_tsv(path: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            parsed[row["metric"]] = row["value"]
    return parsed


def read_signac_stats(path: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            parsed[row["metric"]] = row["value"]
    return parsed


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "Usage: parse_downstream_analysis_outputs.py <summary_tsv> <signac_stats> <sample_name> <bin_size> <output_csv>"
        )
        return 1

    summary_file, signac_stats_file, sample_name, bin_size, output_path = sys.argv[1:6]
    summary = read_key_value_tsv(summary_file)
    signac = read_signac_stats(signac_stats_file)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "metric", "value", "unit", "source", "timestamp"])
        writer.writerow(["downstream", "sample_name", sample_name, "-", f"bin{bin_size}", timestamp])
        writer.writerow(["downstream", "bin_size", bin_size, "-", f"bin{bin_size}", timestamp])
        for metric, value in summary.items():
            writer.writerow(["downstream", metric, value, "-", os.path.basename(summary_file), timestamp])
        for metric, value in signac.items():
            writer.writerow(["downstream", metric, value, "-", os.path.basename(signac_stats_file), timestamp])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
