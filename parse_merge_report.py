#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import os
import sys
from datetime import datetime


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: parse_merge_report.py <report_file> <sample_name> <output_csv>")
        return 1

    report_file, sample_name, output_path = sys.argv[1:4]
    parsed: dict[str, str] = {}
    with open(report_file, "r", encoding="utf-8", errors="replace") as handle:
        content = handle.read().strip()

    if content.startswith("{"):
        loaded = json.loads(content)
        parsed = {key: str(value) for key, value in loaded.items()}
    else:
        for line in content.splitlines():
            line = line.strip()
            if ":" in line:
                key, value = line.split(":", 1)
                parsed[key.strip()] = value.strip()

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    mapping = {
        "total_frag_count": "fragments",
        "total_unique_frag_count": "fragments",
        "saturation": "ratio",
    }
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "metric", "value", "unit", "source", "timestamp"])
        writer.writerow(["merge", "sample_name", sample_name, "-", os.path.basename(report_file), timestamp])
        for metric, unit in mapping.items():
            writer.writerow(["merge", metric, parsed.get(metric, "0"), unit, os.path.basename(report_file), timestamp])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
