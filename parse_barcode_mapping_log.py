#!/usr/bin/env python3
"""
Normalize ST_BarcodeMap log output into a simple metrics table.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from datetime import datetime


METRIC_SPECS = [
    ("total_reads", r"^total_reads:\s*(\d+)", "reads"),
    ("pass_filter_reads", r"^pass_filter_reads:\s*(\d+)", "reads"),
    ("mapped_reads", r"^mapped_reads:\s*(\d+)", "reads"),
    ("mapped_reads_rate", r"^mapped_reads:\s*\d+\s+([\d.]+)%", "%"),
    (
        "barcode_exactly_overlap_reads",
        r"^barcode_exactlyOverlap_reads:\s*(\d+)",
        "reads",
    ),
    (
        "barcode_exactly_overlap_rate",
        r"^barcode_exactlyOverlap_reads:\s*\d+\s+([\d.]+)%",
        "%",
    ),
    ("barcode_mis_overlap_reads", r"^barcode_misOverlap_reads:\s*(\d+)", "reads"),
    ("barcode_mis_overlap_rate", r"^barcode_misOverlap_reads:\s*\d+\s+([\d.]+)%", "%"),
    ("q10_bases_in_barcode", r"^Q10_bases_in_barcode:\s*([\d.]+)%", "%"),
    ("q20_bases_in_barcode", r"^Q20_bases_in_barcode:\s*([\d.]+)%", "%"),
    ("q30_bases_in_barcode", r"^Q30_bases_in_barcode:\s*([\d.]+)%", "%"),
]


def extract_metrics(log_path: str, sample_name: str, lane_label: str) -> list[dict[str, str]]:
    if not os.path.exists(log_path):
        raise FileNotFoundError(f"log file not found: {log_path}")

    found: dict[str, str] = {}
    with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            for metric, pattern, _unit in METRIC_SPECS:
                match = re.search(pattern, line)
                if match:
                    found[metric] = match.group(1)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    rows = [
        {
            "section": "barcode_mapping",
            "metric": "sample_name",
            "value": sample_name,
            "unit": "-",
            "source": os.path.basename(log_path),
            "timestamp": timestamp,
        }
    ]
    rows.append(
        {
            "section": "barcode_mapping",
            "metric": "lane_label",
            "value": lane_label,
            "unit": "-",
            "source": os.path.basename(log_path),
            "timestamp": timestamp,
        }
    )

    for metric, _pattern, unit in METRIC_SPECS:
        rows.append(
            {
                "section": "barcode_mapping",
                "metric": metric,
                "value": found.get(metric, "0"),
                "unit": unit,
                "source": os.path.basename(log_path),
                "timestamp": timestamp,
            }
        )

    return rows


def write_metrics(rows: list[dict[str, str]], output_path: str) -> None:
    fieldnames = ["section", "metric", "value", "unit", "source", "timestamp"]
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, str]]) -> None:
    print("=" * 60)
    print("BARCODE MAPPING METRICS")
    print("=" * 60)
    sample_name = next((row["value"] for row in rows if row["metric"] == "sample_name"), "N/A")
    print(f"sample_name                 : {sample_name}")
    for row in rows:
        if row["metric"] == "sample_name":
            continue
        print(f"{row['metric']:28s}: {row['value']:>12} {row['unit']}")
    print("=" * 60)


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "Usage: python parse_barcode_mapping_log.py <log_file> <sample_name> <lane_label> <output_csv>"
        )
        return 1

    log_path = sys.argv[1]
    sample_name = sys.argv[2]
    lane_label = sys.argv[3]
    output_path = sys.argv[4]

    rows = extract_metrics(log_path, sample_name, lane_label)
    write_metrics(rows, output_path)
    print_summary(rows)
    print(f"metrics_file                : {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
