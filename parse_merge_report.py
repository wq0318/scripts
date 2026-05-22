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

    def parse_kv(text: str) -> dict[str, str]:
        kv: dict[str, str] = {}
        for line in text.splitlines():
            line = line.strip()
            if ":" in line:
                key, value = line.split(":", 1)
                kv[key.strip()] = value.strip()
        return kv

    if content.startswith("{"):
        try:
            loaded = json.loads(content)
            parsed = {key: str(value) for key, value in loaded.items()}
        except json.JSONDecodeError as exc:
            print(f"[WARN] JSON parse failed ({exc}); falling back to KV parser: {report_file}", file=sys.stderr)
            parsed = parse_kv(content)
    else:
        parsed = parse_kv(content)

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
