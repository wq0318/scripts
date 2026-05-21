#!/usr/bin/env python3

from __future__ import annotations

import csv
import io
import os
import re
import sys
from datetime import datetime


def read_summary(path: str) -> tuple[list[str], list[str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    if "Number of" in text:
        return read_chromap_log(text)

    first_line = text.splitlines()[0] if text.splitlines() else ""
    if "," in first_line:
        return read_csv_summary(text, path)

    tokens = [token.strip() for token in text.replace("\n", ",").split(",") if token.strip()]
    if not tokens:
        raise ValueError(f"empty mapping summary: {path}")
    if any(not token.replace(".", "", 1).isdigit() for token in tokens[:5]):
        header = tokens[:5]
        values = tokens[5:10]
    else:
        header = ["barcode", "total", "duplicate", "unmapped", "lowmapq"]
        values = tokens[:5]
    if len(values) < len(header):
        values.extend(["0"] * (len(header) - len(values)))
    return header, values


def read_csv_summary(text: str, path: str) -> tuple[list[str], list[str]]:
    reader = csv.DictReader(io.StringIO(text))
    row = next(reader, None)
    if row is None:
        raise ValueError(f"empty mapping summary: {path}")
    header = reader.fieldnames or []
    values = [row.get(field, "0") or "0" for field in header]
    return header, values


def read_chromap_log(text: str) -> tuple[list[str], list[str]]:
    patterns = {
        "total": r"Number of reads:\s*([\d,]+)",
        "mapped": r"Number of mapped reads:\s*([\d,]+)",
        "unique": r"Number of uniquely mapped reads:\s*([\d,]+)",
    }
    parsed: dict[str, str] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        parsed[key] = match.group(1).replace(",", "") if match else "0"

    total = float(parsed["total"] or 0)
    mapped = float(parsed["mapped"] or 0)
    unmapped = max(total - mapped, 0)
    return (
        ["barcode", "total", "duplicate", "unmapped", "lowmapq", "unique"],
        ["0", str(int(total)), "0", str(int(unmapped)), "0", str(int(float(parsed["unique"] or 0)))],
    )


def main() -> int:
    if len(sys.argv) != 5:
        print("Usage: parse_mapping_summary.py <summary_csv> <sample_name> <lane_label> <output_csv>")
        return 1

    summary_path, sample_name, lane_label, output_path = sys.argv[1:5]
    header, values = read_summary(summary_path)
    parsed = dict(zip(header, values))
    total = float(parsed.get("total", "0") or 0)
    unmapped = float(parsed.get("unmapped", "0") or 0)
    lowmapq = float(parsed.get("lowmapq", "0") or 0)
    duplicate = float(parsed.get("duplicate", "0") or 0)
    mapped = max(total - unmapped, 0)
    high_quality = float(parsed.get("unique", "0") or 0) if "unique" in parsed else max(total - unmapped - lowmapq, 0)

    rows = [
        ("mapping", "sample_name", sample_name, "-", lane_label),
        ("mapping", "lane_label", lane_label, "-", lane_label),
        ("mapping", "chromap_total_reads", f"{int(total)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_unmapped_reads", f"{int(unmapped)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_duplicate_reads", f"{int(duplicate)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_lowmapq_reads", f"{int(lowmapq)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_mapped_reads", f"{int(mapped)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_high_quality_reads", f"{int(high_quality)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_mapping_rate", f"{(mapped * 100 / total) if total else 0:.2f}", "%", os.path.basename(summary_path)),
        ("mapping", "chromap_high_quality_rate", f"{(high_quality * 100 / total) if total else 0:.2f}", "%", os.path.basename(summary_path)),
    ]

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "metric", "value", "unit", "source", "timestamp"])
        for row in rows:
            writer.writerow([*row, timestamp])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
