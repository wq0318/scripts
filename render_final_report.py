#!/usr/bin/env python3

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


def parse_path_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def read_metrics(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def read_cluster_plot_set(value: str) -> list[str]:
    if value and os.path.exists(value):
        with open(value, "r", encoding="utf-8", errors="replace") as handle:
            value = handle.read().strip()
    return parse_path_list(value)


def read_saturation_table(path: str) -> list[dict[str, float]]:
    if not path or not os.path.exists(path):
        return []
    rows: list[dict[str, float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        for line in handle:
            fields = line.strip().split()
            if len(fields) < 5:
                continue
            try:
                rows.append({
                    "percent": float(fields[0]),
                    "unique": float(fields[1]),
                    "total": float(fields[2]),
                    "dup_pct": float(fields[3]),
                    "median_nfrags": float(fields[4]),
                })
            except ValueError:
                continue
    return rows


def mime_type_for_suffix(path: str) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".svg":
        return "image/svg+xml"
    if suffix in {".jpg", ".jpeg"}:
        return "image/jpeg"
    if suffix == ".webp":
        return "image/webp"
    return "image/png"


def placeholder_data_uri(title: str) -> str:
    safe_title = html.escape(title)
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 400">'
        '<rect width="800" height="400" fill="#f7fafc" stroke="#e2e8f0"/>'
        f'<text x="400" y="200" text-anchor="middle" font-family="Arial" font-size="18" fill="#718096">{safe_title}</text>'
        '</svg>'
    )
    return "data:image/svg+xml;charset=utf-8," + quote(svg)


def embed_image(path: str, title: str) -> str:
    if not path or not os.path.exists(path):
        return placeholder_data_uri(title)
    with open(path, "rb") as handle:
        encoded = base64.b64encode(handle.read()).decode("ascii")
    return f"data:{mime_type_for_suffix(path)};base64,{encoded}"


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def format_int(value: str, fallback: str = "0") -> str:
    parsed = to_float(value)
    return fallback if parsed is None else f"{int(round(parsed)):,}"


def format_num(value: str, digits: int = 2, suffix: str = "", fallback: str = "0") -> str:
    parsed = to_float(value)
    return fallback if parsed is None else f"{parsed:.{digits}f}{suffix}"


def aggregate_global_metrics(rows: list[dict[str, str]]) -> dict[str, str]:
    numeric_sums: defaultdict[str, float] = defaultdict(float)
    latest_values: dict[str, str] = {}

    for row in rows:
        metric = row["metric"]
        value = row["value"]
        latest_values[metric] = value
        parsed = to_float(value)
        if parsed is None:
            continue
        numeric_sums[metric] += parsed

    total_reads = numeric_sums.get("total_reads", 0)
    mapped_reads = numeric_sums.get("mapped_reads", 0)
    exact_reads = numeric_sums.get("barcode_exactly_overlap_reads", 0)
    mis_reads = numeric_sums.get("barcode_mis_overlap_reads", 0)
    mapping_total = numeric_sums.get("chromap_total_reads", 0)
    mapping_hq = numeric_sums.get("chromap_high_quality_reads", 0)

    latest_values["qc_total_reads"] = str(int(total_reads))
    latest_values["qc_valid_reads"] = str(int(mapped_reads))
    latest_values["qc_barcode_mapping_rate"] = f"{(mapped_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["qc_exact_overlap_rate"] = f"{(exact_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["qc_mismatch_overlap_rate"] = f"{(mis_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["qc_mapping_total_reads"] = str(int(mapping_total))
    latest_values["qc_genome_mapping_rate"] = f"{(mapping_total * 100 / mapped_reads) if mapped_reads else 0:.2f}"
    latest_values["qc_high_quality_mapping_rate"] = f"{(mapping_hq * 100 / mapping_total) if mapping_total else 0:.2f}"
    latest_values["qc_low_quality_mapping_rate"] = f"{((mapping_total - mapping_hq) * 100 / mapping_total) if mapping_total else 0:.2f}"
    latest_values["qc_mitochondria_ratio"] = latest_values.get("mitochondria_ratio", "0")
    latest_values["qc_saturation"] = latest_values.get("saturation", "0")

    return latest_values


def aggregate_bin_metrics(rows: list[dict[str, str]]) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for row in rows:
        metrics[row["metric"]] = row["value"]
    return metrics


def read_nucleosome_stats(path: str) -> dict[str, str]:
    stats: dict[str, str] = {}
    if not path or not os.path.exists(path):
        return stats
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if "metric" in row and "value" in row:
                stats[row["metric"]] = row["value"]
    return stats


def render_section1_overview(qc: dict[str, str]) -> str:
    total_reads = format_int(qc.get("qc_total_reads", "0"))
    barcode_rate = format_num(qc.get("qc_barcode_mapping_rate", "0"), suffix="%")
    exact_rate = format_num(qc.get("qc_exact_overlap_rate", "0"), suffix="%")
    mismatch_rate = format_num(qc.get("qc_mismatch_overlap_rate", "0"), suffix="%")
    genome_rate = format_num(qc.get("qc_genome_mapping_rate", "0"), suffix="%")
    hq_rate = format_num(qc.get("qc_high_quality_mapping_rate", "0"), suffix="%")
    lq_rate = format_num(qc.get("qc_low_quality_mapping_rate", "0"), suffix="%")
    saturation_pct = format_num(qc.get("qc_saturation", "0"), suffix="%")
    mito_rate = format_num(qc.get("qc_mitochondria_ratio", "0"), suffix="%")

    return f"""
    <h2>1. Alignment &amp; Quality Overview</h2>
    <div class="qc-overview-container">
        <div class="qc-table-box">
            <div class="summary-row">
                <span class="summary-label">Total Reads</span><span class="summary-val">{total_reads}</span>
            </div>
            <details>
                <summary><span class="summary-label">Barcode Mapping Rate</span><span class="summary-val">{barcode_rate}</span></summary>
                <div class="nested-content">
                    <div class="nested-item"><span class="summary-label">Exact Overlap Rate</span><span>{exact_rate}</span></div>
                    <div class="nested-item"><span class="summary-label">Mismatch Overlap Rate</span><span>{mismatch_rate}</span></div>
                </div>
            </details>
            <details>
                <summary><span class="summary-label">Genome Mapping Rate</span><span class="summary-val">{genome_rate}</span></summary>
                <div class="nested-content">
                    <div class="nested-item"><span class="summary-label">High-quality Mapping Rate</span><span>{hq_rate}</span></div>
                    <div class="nested-item"><span class="summary-label">Low-quality Mapping Rate</span><span>{lq_rate}</span></div>
                </div>
            </details>
        </div>
        <div class="qc-table-box">
            <div class="summary-row">
                <span class="summary-label">Sequencing Saturation</span><span class="summary-val">{saturation_pct}</span>
            </div>
            <div class="summary-row">
                <span class="summary-label">Mitochondria Ratio</span><span class="summary-val">{mito_rate}</span>
            </div>
        </div>
    </div>
    """


def render_section2_saturation(saturation_per_bin: dict[str, list[dict[str, float]]]) -> str:
    js_payload = {}
    for bin_label, rows in saturation_per_bin.items():
        js_payload[bin_label] = [
            {
                "total_m": round(row["total"] / 1_000_000.0, 3),
                "dup_pct": round(row["dup_pct"], 2),
                "median_nfrags": round(row["median_nfrags"], 1),
            }
            for row in rows
        ]
    data_json = json.dumps(js_payload)

    return f"""
    <h2>2. Sequencing Saturation</h2>
    <div class="chart-grid">
        <div class="chart-card">
            <span class="chart-title">Sequencing Saturation Plot</span>
            <div id="saturationChart" class="chart-container"></div>
        </div>
        <div class="chart-card">
            <span class="chart-title">Library Sensitivity Plot</span>
            <div id="sensitivityChart" class="chart-container"></div>
        </div>
    </div>
    <script>
      window.__saturationData = {data_json};
    </script>
    """


def render_section3_fragment(bin_labels: list[str], bin_metrics_map: dict[str, tuple[dict[str, str], dict[str, str]]],
                              fragment_size_images: dict[str, str]) -> str:
    primary = bin_labels[0]
    _, ns = bin_metrics_map[primary]
    nfr = format_num(ns.get("fraction_nfr", "0"), suffix="%")
    mono = format_num(ns.get("fraction_mono_nucleosome", "0"), suffix="%")
    multi = format_num(ns.get("fraction_multi_nucleosome", "0"), suffix="%")

    return f"""
    <h2>3. Fragment Size Distribution <span class="bin-context">({primary})</span></h2>
    <div class="qc-overview-container">
        <div style="flex: 0 0 350px;">
            <span class="chart-title" style="text-align: left;">Fragment Statistics</span>
            <table class="metrics-table">
                <tr><td align="left">Fraction of NFR</td><td>{nfr}</td></tr>
                <tr><td align="left">Fraction of Mono-nucleosome</td><td>{mono}</td></tr>
                <tr><td align="left">Fraction of Multi-nucleosome</td><td>{multi}</td></tr>
            </table>
        </div>
        <div class="chart-card" style="flex: 1;">
            <span class="chart-title">Insert Size Distribution</span>
            <img class="full-img" src="{fragment_size_images[primary]}" alt="{primary} fragment size" />
        </div>
    </div>
    """


def render_section4_downstream(bin_labels: list[str],
                                bin_metrics_map: dict[str, tuple[dict[str, str], dict[str, str]]],
                                bin_images: dict[str, dict[str, str]],
                                cluster_map: dict[str, dict[str, str]]) -> str:
    header = "".join([
        "<th>Resolution</th>",
        "<th>Total spots</th>",
        "<th>Valid spots</th>",
        "<th>Median nFrags</th>",
        "<th>Median TSS</th>",
        "<th>Median FRiP</th>",
    ])
    body_rows = []
    for bl in bin_labels:
        dm, ns = bin_metrics_map[bl]
        body_rows.append(
            "<tr>"
            f"<td><strong>Bin {bl.replace('bin', '')}</strong></td>"
            f"<td>{format_int(dm.get('total_barcodes', '0'))}</td>"
            f"<td>{format_int(dm.get('valid_barcodes', '0'))}</td>"
            f"<td>{format_int(dm.get('median_nFrags', '0'))}</td>"
            f"<td>{format_num(dm.get('median_TSS', '0'), digits=2)}</td>"
            f"<td>{format_num(dm.get('median_FRiP', '0'), digits=2)}</td>"
            "</tr>"
        )
    table_html = (
        '<table class="metrics-table">'
        f'<thead><tr>{header}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody>'
        '</table>'
    )

    tab_btns = "".join(
        f'<button class="tab-btn{" active" if i == 0 else ""}" data-tab="{bl}">{bl.replace("bin", "Bin ")} View</button>'
        for i, bl in enumerate(bin_labels)
    )

    tab_contents = []
    resolutions = ["0.4", "0.6", "0.8", "1.0", "1.2"]
    for i, bl in enumerate(bin_labels):
        active_cls = " active" if i == 0 else ""
        images = bin_images[bl]
        cluster_entries = ", ".join(
            f'"{res}": "{cluster_map[bl].get(res, "")}"' for res in resolutions
        )
        tab_contents.append(f"""
        <div class="tab-content{active_cls}" id="{bl}-content">
            <div class="bin-toolbar">
                <button class="layout-btn active" data-layout="post" data-bin="{bl}">Post-cut View</button>
                <button class="layout-btn" data-layout="pre" data-bin="{bl}">Pre-cut View</button>
            </div>

            <div class="chart-grid layout-stage" data-bin="{bl}">
                <div class="chart-card">
                    <span class="chart-title">Spatial nFrags</span>
                    <img class="full-img stage-img" data-img-post="{images['spatial']}" data-img-pre="{images['spatial_pre']}" src="{images['spatial']}" alt="{bl} spatial nFrags" />
                </div>
                <div class="chart-card">
                    <span class="chart-title">TSS Enrichment Scatter</span>
                    <img class="full-img" src="{images['tss']}" alt="{bl} TSS scatter" />
                </div>
            </div>

            <div class="chart-grid" style="margin-top: 20px;">
                <div class="chart-card"><span class="chart-title">QC Violin Plots</span><img class="full-img" src="{images['qc']}" alt="{bl} qc violin" /></div>
                <div class="chart-card"><span class="chart-title">Fragment Size</span><img class="full-img" src="{images['fragment']}" alt="{bl} fragment size" /></div>
            </div>

            <div class="cluster-control-panel">
                <div class="slider-box">
                    <strong>Cluster Resolution:</strong>
                    <input type="range" class="res-slider" min="0" max="4" step="1" value="2" data-bin="{bl}" aria-label="{bl} cluster resolution" />
                    <span class="res-val" data-bin="{bl}">0.8</span>
                </div>
                <div class="chart-card cluster-canvas">
                    <span class="chart-title">UMAP + Spatial Clusters</span>
                    <img class="full-img cluster-img" data-bin="{bl}" src="{cluster_map[bl].get('0.8', '')}" alt="{bl} cluster" />
                </div>
                <script>
                  window.__clusterMaps = window.__clusterMaps || {{}};
                  window.__clusterMaps["{bl}"] = {{{cluster_entries}}};
                </script>
            </div>
        </div>
        """)

    return f"""
    <h2>4. Downstream Analysis Metrics</h2>
    {table_html}

    <div class="interactive-container">
        <div class="tab-buttons">
            {tab_btns}
        </div>
        {"".join(tab_contents)}
    </div>
    """


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render multi-bin StereoATAC HTML report aligned with preview.html template"
    )
    parser.add_argument("sample_name")
    parser.add_argument("output_html")
    parser.add_argument("barcode_metrics", help="Comma-separated barcode mapping metrics CSV files")
    parser.add_argument("mapping_metrics", help="Comma-separated mapping metrics CSV files")
    parser.add_argument("sam2frag_metrics", help="Comma-separated sam2frag metrics CSV files")
    parser.add_argument("merge_metrics", help="Merge metrics CSV file")
    parser.add_argument("--bin-sizes", required=True)
    parser.add_argument("--downstream-metrics", required=True)
    parser.add_argument("--nucleosome-stats", required=True)
    parser.add_argument("--saturation-plots", required=True, help="Per-bin Saturation.svg (legacy, unused in template)")
    parser.add_argument("--saturation-tables", required=True, help="Per-bin _result.txt from atac_saturation.pl")
    parser.add_argument("--fragment-stats-plots", required=True)
    parser.add_argument("--qc-violins", required=True)
    parser.add_argument("--tss-scatters", required=True)
    parser.add_argument("--fragment-sizes", required=True)
    parser.add_argument("--spatial-qcs", required=True, help="Per-bin post-cut spatial_qc.svg")
    parser.add_argument("--spatial-qcs-pre", required=True, help="Per-bin pre-cut spatial_qc_pre.svg")
    parser.add_argument("--cluster-plots-list", required=True)
    args = parser.parse_args()

    bin_sizes = [s.strip() for s in args.bin_sizes.split(",") if s.strip()]
    bin_labels = [f"bin{s}" for s in bin_sizes]
    num_bins = len(bin_labels)

    qc_rows = read_metrics([
        *parse_path_list(args.barcode_metrics),
        *parse_path_list(args.mapping_metrics),
        *parse_path_list(args.sam2frag_metrics),
        args.merge_metrics,
    ])
    qc_metrics = aggregate_global_metrics(qc_rows)

    def parse_per_bin_list(value: str, expected: int) -> list[str]:
        items = [s.strip() for s in value.split(",") if s.strip()]
        if len(items) != expected:
            raise ValueError(f"Expected {expected} items, got {len(items)}: {value}")
        return items

    downstream_list = parse_per_bin_list(args.downstream_metrics, num_bins)
    nucleosome_list = parse_per_bin_list(args.nucleosome_stats, num_bins)
    saturation_table_list = parse_per_bin_list(args.saturation_tables, num_bins)
    qc_violins = parse_per_bin_list(args.qc_violins, num_bins)
    tss_scatters = parse_per_bin_list(args.tss_scatters, num_bins)
    fragment_sizes = parse_per_bin_list(args.fragment_sizes, num_bins)
    spatial_qcs = parse_per_bin_list(args.spatial_qcs, num_bins)
    spatial_qcs_pre = parse_per_bin_list(args.spatial_qcs_pre, num_bins)

    cluster_plots_per_bin = [s.strip() for s in args.cluster_plots_list.split(";")]
    if len(cluster_plots_per_bin) != num_bins:
        raise ValueError(f"Expected {num_bins} cluster plot sets, got {len(cluster_plots_per_bin)}")

    bin_images: dict[str, dict[str, str]] = {}
    cluster_map: dict[str, dict[str, str]] = {}
    bin_metrics_map: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    saturation_per_bin: dict[str, list[dict[str, float]]] = {}

    resolutions = ["0.4", "0.6", "0.8", "1.0", "1.2"]
    fragment_size_images: dict[str, str] = {}

    for i, bl in enumerate(bin_labels):
        bin_images[bl] = {
            "qc": embed_image(qc_violins[i], f"{bl} QC violin"),
            "tss": embed_image(tss_scatters[i], f"{bl} TSS scatter"),
            "fragment": embed_image(fragment_sizes[i], f"{bl} fragment size"),
            "spatial": embed_image(spatial_qcs[i], f"{bl} spatial QC post-cut"),
            "spatial_pre": embed_image(spatial_qcs_pre[i], f"{bl} spatial QC pre-cut"),
        }
        fragment_size_images[bl] = bin_images[bl]["fragment"]

        cluster_paths = read_cluster_plot_set(cluster_plots_per_bin[i])
        cluster_map[bl] = {}
        for j, res in enumerate(resolutions):
            if j < len(cluster_paths):
                cluster_map[bl][res] = embed_image(cluster_paths[j], f"{bl} cluster res={res}")
            else:
                cluster_map[bl][res] = placeholder_data_uri(f"{bl} cluster res={res}")

        dm = aggregate_bin_metrics(read_metrics([downstream_list[i]]))
        ns = read_nucleosome_stats(nucleosome_list[i])
        bin_metrics_map[bl] = (dm, ns)

        saturation_per_bin[bl] = read_saturation_table(saturation_table_list[i])

    section1 = render_section1_overview(qc_metrics)
    section2 = render_section2_saturation(saturation_per_bin)
    section3 = render_section3_fragment(bin_labels, bin_metrics_map, fragment_size_images)
    section4 = render_section4_downstream(bin_labels, bin_metrics_map, bin_images, cluster_map)

    bin_labels_json = json.dumps(bin_labels)
    default_bin = bin_labels[0]
    title = f"{args.sample_name} StereoATAC-seq Quality Control Report"

    html_text = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{html.escape(title)}</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.4.3/dist/echarts.min.js"></script>
    <style>
        :root {{
            --primary-color: #1a365d;
            --secondary-color: #2d3748;
            --accent-color: #3182ce;
            --bg-color: #f7fafc;
            --border-color: #e2e8f0;
            --highlight-bg: #edf2f7;
            --card-bg: #ffffff;
            --text-primary: #1a202c;
            --text-secondary: #4a5568;
            --text-muted: #718096;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
            font-size: 20px;
            color: var(--text-primary);
            background-color: var(--bg-color);
            line-height: 1.6;
            margin: 0;
            padding: 30px 20px;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background: var(--card-bg);
            padding: 35px 40px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08), 0 4px 16px rgba(0,0,0,0.04);
            border: 1px solid var(--border-color);
        }}
        .main-title {{
            font-size: 26px;
            font-weight: 600;
            margin-bottom: 8px;
            color: var(--primary-color);
            letter-spacing: 0.5px;
        }}
        .sub-title {{
            font-size: 17px;
            color: var(--text-muted);
            margin-bottom: 30px;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 15px;
        }}
        h2 {{
            font-size: 23px;
            font-weight: 600;
            color: var(--primary-color);
            border-left: 4px solid var(--accent-color);
            padding-left: 14px;
            margin-top: 40px;
            margin-bottom: 18px;
            background: linear-gradient(to right, var(--highlight-bg), transparent);
            padding-top: 8px;
            padding-bottom: 8px;
        }}
        .bin-context {{
            font-size: 14px;
            color: var(--text-muted);
            font-weight: 400;
        }}
        .qc-overview-container {{
            display: flex;
            gap: 40px;
            align-items: flex-start;
        }}
        .qc-table-box {{ flex: 1; }}
        .summary-row {{
            padding: 10px 0;
            font-size: 17px;
            display: flex;
            justify-content: space-between;
            border-bottom: 1px solid #f0f0f0;
        }}
        details {{
            margin-bottom: 0;
            cursor: pointer;
            border-bottom: 1px solid var(--border-color);
            background: white;
        }}
        details:hover {{ background: var(--highlight-bg); }}
        summary {{
            padding: 12px 0;
            font-size: 17px;
            list-style: none;
            display: flex;
            align-items: center;
            justify-content: space-between;
            color: var(--text-secondary);
            font-weight: 500;
        }}
        summary::-webkit-details-marker {{ display: none; }}
        summary::before {{
            content: "▶";
            font-size: 10px;
            margin-right: 10px;
            transition: transform 0.2s;
            color: var(--text-muted);
        }}
        details[open] summary::before {{ transform: rotate(90deg); }}
        .summary-label {{ flex: 1; }}
        .summary-val {{
            font-weight: 600;
            color: var(--accent-color);
            font-family: "Arial", "Helvetica", sans-serif;
        }}
        .nested-content {{
            padding-left: 30px;
            padding-bottom: 12px;
            font-size: 16px;
            color: var(--text-secondary);
            background: var(--highlight-bg);
        }}
        .nested-item {{
            display: flex;
            justify-content: space-between;
            padding: 5px 0;
        }}
        .chart-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 25px;
            margin-top: 20px;
        }}
        .chart-card {{
            border: 1px solid var(--border-color);
            padding: 20px;
            background: var(--card-bg);
            border-radius: 4px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.06);
        }}
        .chart-title {{
            font-size: 17px;
            font-weight: 600;
            margin-bottom: 15px;
            display: block;
            text-align: center;
            color: var(--text-secondary);
        }}
        .chart-container {{
            width: 100%;
            height: 250px;
        }}
        .full-img {{
            display: block;
            width: 100%;
            height: auto;
        }}
        .metrics-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 16px;
        }}
        .metrics-table th, .metrics-table td {{
            border: 1px solid var(--border-color);
            padding: 10px 12px;
            text-align: center;
        }}
        .metrics-table th {{
            background: var(--highlight-bg);
            font-weight: 600;
            color: var(--text-primary);
        }}
        .metrics-table tr:hover {{ background: var(--highlight-bg); }}
        .interactive-container {{
            border: 1px solid var(--border-color);
            margin-top: 25px;
            padding: 25px;
            border-radius: 6px;
            background: var(--card-bg);
            box-shadow: 0 1px 3px rgba(0,0,0,0.06);
        }}
        .tab-buttons {{
            display: flex;
            gap: 12px;
            margin-bottom: 25px;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 12px;
        }}
        .tab-btn, .layout-btn {{
            padding: 8px 24px;
            border: 1px solid var(--border-color);
            background: white;
            cursor: pointer;
            transition: all 0.2s;
            font-size: 16px;
            font-weight: 500;
            color: var(--text-secondary);
            border-radius: 3px;
        }}
        .tab-btn:hover, .layout-btn:hover {{
            background: var(--highlight-bg);
            border-color: var(--accent-color);
        }}
        .tab-btn.active, .layout-btn.active {{
            background: var(--accent-color);
            color: white;
            border-color: var(--accent-color);
        }}
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
        .bin-toolbar {{
            display: flex;
            gap: 8px;
            margin-bottom: 20px;
        }}
        .cluster-control-panel {{
            margin-top: 30px;
            background: var(--highlight-bg);
            padding: 20px;
            border-radius: 6px;
            border: 1px solid var(--border-color);
        }}
        .slider-box {{
            display: flex;
            align-items: center;
            gap: 20px;
            background: white;
            padding: 12px 25px;
            margin: 0 0 15px 0;
            border-radius: 4px;
            font-size: 17px;
            border: 1px solid var(--border-color);
        }}
        .slider-box input[type="range"] {{
            flex: 0 0 200px;
            accent-color: var(--accent-color);
        }}
        .res-val {{
            font-family: monospace;
            font-weight: bold;
            color: var(--accent-color);
            min-width: 30px;
        }}
        .cluster-canvas {{
            margin-top: 12px;
        }}
        footer {{
            margin-top: 50px;
            text-align: center;
            color: var(--text-muted);
            font-size: 15px;
            padding-top: 20px;
            border-top: 1px solid var(--border-color);
        }}
    </style>
</head>
<body>
<div class="container">
    <div class="main-title">{html.escape(title)}</div>
    <div class="sub-title">Sample: {html.escape(args.sample_name)} &nbsp;|&nbsp; Bins: {", ".join(bin_labels)}</div>

    {section1}
    {section2}
    {section3}
    {section4}

    <footer>StereoATAC-seq Pipeline | Standardized QC Report</footer>
</div>

<script>
    (function() {{
        const satData = window.__saturationData || {{}};
        const binLabels = {bin_labels_json};
        const defaultBin = "{default_bin}";

        const pal = {{ saturation: "#3182ce", sensitivity: "#38a169" }};

        let saturationChart, sensitivityChart;
        function renderSaturation(bin) {{
            const rows = satData[bin] || [];
            const xAxis = rows.map(r => r.total_m);
            const dup = rows.map(r => r.dup_pct);
            const median = rows.map(r => r.median_nfrags);

            const baseGrid = {{ left: '12%', right: '8%', bottom: '18%', top: '12%' }};
            const xAxisCommon = {{
                type: 'category',
                name: 'Total Fragments (M)',
                nameLocation: 'middle',
                nameGap: 28,
                data: xAxis,
                axisLine: {{ lineStyle: {{ color: '#666' }} }},
                axisLabel: {{ color: '#444', fontSize: 11 }}
            }};

            if (!saturationChart) saturationChart = echarts.init(document.getElementById('saturationChart'));
            saturationChart.setOption({{
                tooltip: {{ trigger: 'axis', axisPointer: {{ type: 'cross' }} }},
                grid: baseGrid,
                xAxis: xAxisCommon,
                yAxis: {{
                    type: 'value',
                    name: 'Duplicate Rate (%)',
                    nameLocation: 'middle',
                    nameGap: 40,
                    axisLine: {{ lineStyle: {{ color: '#666' }} }},
                    axisLabel: {{ formatter: '{{value}}%' }},
                    splitLine: {{ lineStyle: {{ color: '#e2e8f0', type: 'dashed' }} }}
                }},
                series: [{{
                    name: 'Duplicate %',
                    type: 'line', smooth: true, symbol: 'circle', symbolSize: 6,
                    data: dup,
                    lineStyle: {{ color: pal.saturation, width: 2.5 }},
                    itemStyle: {{ color: pal.saturation }},
                    areaStyle: {{
                        color: new echarts.graphic.LinearGradient(0, 0, 0, 1, [
                            {{ offset: 0, color: 'rgba(49, 130, 206, 0.3)' }},
                            {{ offset: 1, color: 'rgba(49, 130, 206, 0.05)' }}
                        ])
                    }}
                }}]
            }}, true);

            if (!sensitivityChart) sensitivityChart = echarts.init(document.getElementById('sensitivityChart'));
            sensitivityChart.setOption({{
                tooltip: {{ trigger: 'axis', axisPointer: {{ type: 'cross' }} }},
                grid: baseGrid,
                xAxis: xAxisCommon,
                yAxis: {{
                    type: 'value',
                    name: 'Median nFrags / Bin',
                    nameLocation: 'middle',
                    nameGap: 48,
                    axisLine: {{ lineStyle: {{ color: '#666' }} }},
                    splitLine: {{ lineStyle: {{ color: '#e2e8f0', type: 'dashed' }} }}
                }},
                series: [{{
                    name: 'Median nFrags',
                    type: 'line', smooth: true, symbol: 'circle', symbolSize: 6,
                    data: median,
                    lineStyle: {{ color: pal.sensitivity, width: 2.5 }},
                    itemStyle: {{ color: pal.sensitivity }},
                    areaStyle: {{
                        color: new echarts.graphic.LinearGradient(0, 0, 0, 1, [
                            {{ offset: 0, color: 'rgba(56, 161, 105, 0.3)' }},
                            {{ offset: 1, color: 'rgba(56, 161, 105, 0.05)' }}
                        ])
                    }}
                }}]
            }}, true);
        }}

        renderSaturation(defaultBin);
        window.addEventListener('resize', function() {{
            if (saturationChart) saturationChart.resize();
            if (sensitivityChart) sensitivityChart.resize();
        }});

        const tabButtons = document.querySelectorAll('.tab-btn');
        const tabContents = document.querySelectorAll('.tab-content');
        function switchBin(tabId) {{
            tabButtons.forEach(b => b.classList.toggle('active', b.dataset.tab === tabId));
            tabContents.forEach(c => c.classList.toggle('active', c.id === tabId + '-content'));
            renderSaturation(tabId);
        }}
        tabButtons.forEach(b => b.addEventListener('click', () => switchBin(b.dataset.tab)));

        document.querySelectorAll('.layout-btn').forEach(btn => {{
            btn.addEventListener('click', function() {{
                const bin = this.dataset.bin;
                const layout = this.dataset.layout;
                document.querySelectorAll(`.layout-btn[data-bin="${{bin}}"]`).forEach(b =>
                    b.classList.toggle('active', b === this));
                document.querySelectorAll(`.stage-img`).forEach(img => {{
                    const parentBin = img.closest('.tab-content').id.replace('-content', '');
                    if (parentBin !== bin) return;
                    const src = layout === 'pre' ? img.dataset.imgPre : img.dataset.imgPost;
                    if (src) img.src = src;
                    img.alt = `${{bin}} spatial nFrags (${{layout}}-cut)`;
                }});
            }});
        }});

        const clusterMaps = window.__clusterMaps || {{}};
        const resolutionLabels = ["0.4", "0.6", "0.8", "1.0", "1.2"];
        document.querySelectorAll('.res-slider').forEach(slider => {{
            const bin = slider.dataset.bin;
            const valEl = document.querySelector(`.res-val[data-bin="${{bin}}"]`);
            const imgEl = document.querySelector(`.cluster-img[data-bin="${{bin}}"]`);
            slider.addEventListener('input', function() {{
                const res = resolutionLabels[Number(this.value)];
                valEl.textContent = res;
                const src = (clusterMaps[bin] || {{}})[res];
                if (src) imgEl.src = src;
                imgEl.alt = `${{bin}} cluster res=${{res}}`;
            }});
        }});

        switchBin(defaultBin);
    }})();
</script>
</body>
</html>
"""

    Path(args.output_html).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_html, "w", encoding="utf-8") as handle:
        handle.write(html_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
