#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import re
import sys
from datetime import datetime


PATTERNS = [
    ("raw_reads", r"^Raw reads,(\d+)", "reads"),
    ("mapped_reads", r"^Mapped reads,(\d+)", "reads"),
    ("mapped_reads_rate", r"^Mapped reads,\d+ \(([\d.]+)%\)", "%"),
    ("properly_paired_reads", r"^Properly paired reads,(\d+)", "reads"),
    ("singleton_reads", r"^Singleton reads,(\d+)", "reads"),
    ("mitochondria_ratio", r"^Mitochondria ratio,([\d.]+)%", "%"),
]


def main() -> int:
    if len(sys.argv) != 5:
        print("Usage: parse_sam2frag_report.py <report_file> <sample_name> <lane_label> <output_csv>")
        return 1

    report_file, sample_name, lane_label, output_path = sys.argv[1:5]
    found: dict[str, str] = {}
    with open(report_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            for metric, pattern, _unit in PATTERNS:
                match = re.search(pattern, line)
                if match:
                    found[metric] = match.group(1)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "metric", "value", "unit", "source", "timestamp"])
        writer.writerow(["sam2frag", "sample_name", sample_name, "-", lane_label, timestamp])
        writer.writerow(["sam2frag", "lane_label", lane_label, "-", lane_label, timestamp])
        for metric, _pattern, unit in PATTERNS:
            writer.writerow(["sam2frag", metric, found.get(metric, "0"), unit, os.path.basename(report_file), timestamp])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
