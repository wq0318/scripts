#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026-04-21
@author: wq

Version: v1.2

Changelog:
- v1.0: 初始优化版本，修复弃用API，添加日志，增强健壮性
- v1.1: 增加 bgzip 压缩和 tabix 索引功能
- v1.2: 性能优化 - 保留 bin_x/bin_y 避免字符串拆分，添加压缩选项
"""

import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import shutil
from datetime import datetime

import polars as pl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 配置日志
logger = logging.getLogger('binSize')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# 设置Matplotlib缓存目录
if os.path.exists("/data/temp") and os.access("/data/temp", os.W_OK):
    os.environ["MPLCONFIGDIR"] = "/data/temp"
else:
    mpl_temp = tempfile.mkdtemp(prefix="mpl_")
    os.environ["MPLCONFIGDIR"] = mpl_temp
    logger.info(f"Using temp matplotlib cache dir: {mpl_temp}")

parser = argparse.ArgumentParser(
    description='Create binned fragment files, spatial distribution plot, and optionally compress+index',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python binSize.py -i fragment.tsv -o ./out -s sample1 -b 50
  python binSize.py -i fragment.tsv.gz -o ./out -s sample1 -b 100 -l ./binSize.log
  python binSize.py -i fragment.tsv -o ./out -s sample1 -b 50 --no-compress --keep-intermediate
""")
parser.add_argument('-i', '--input', type=str, required=True,
                    help='Input fragment file (bin1 resolution, .tsv or .tsv.gz)')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Output directory')
parser.add_argument('-b', '--binSize', type=int, required=True,
                    help='Bin size for aggregation (must be > 0)')
parser.add_argument('-s', '--sample', type=str, required=True,
                    help='Sample name prefix for output files')
parser.add_argument('-l', '--log', type=str, default=None,
                    help='Log file path (optional)')
parser.add_argument('--no-compress', action='store_true',
                    help='Skip bgzip compression and tabix indexing')
parser.add_argument('--keep-intermediate', action='store_true',
                    help='Keep uncompressed TSV file after compression (ignored if --no-compress)')
args = parser.parse_args()

if args.log:
    file_handler = logging.FileHandler(args.log, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

logger.info(f"=== BinSize Analysis Started at {datetime.now()} ===")
logger.info(f"Sample: {args.sample}, BinSize: {args.binSize}")
logger.info(f"Input: {args.input}")
logger.info(f"Output: {args.output}")

# 参数验证
if args.binSize <= 0:
    logger.error("--binSize must be > 0")
    sys.exit(1)

if not os.path.exists(args.input):
    logger.error(f"Input file not found: {args.input}")
    sys.exit(1)

def check_dependencies():
    missing = []
    for tool in ['bgzip', 'tabix']:
        if not shutil.which(tool):
            missing.append(tool)
    if missing:
        logger.error(f"Missing tools for compression: {', '.join(missing)}")
        logger.error("Install samtools (for bgzip/tabix)")
        sys.exit(1)

if not args.no_compress:
    check_dependencies()

os.makedirs(args.output, exist_ok=True)

# --------------------- 读取数据 ---------------------
logger.info("Reading input file...")
start_time = time.time()
try:
    # 指定列数据类型，避免自动推断，略微提速
    df = pl.read_csv(
        args.input,
        separator="\t",
        has_header=False,
        dtypes=[pl.Utf8, pl.Int64, pl.Int64, pl.Utf8, pl.Int64]
    )
    logger.info(f"Loaded {df.height} rows in {time.time() - start_time:.2f}s")
except Exception as e:
    logger.error(f"Failed to read input file: {e}")
    sys.exit(1)

# --------------------- 解析坐标 ---------------------
logger.info("Parsing coordinates...")
# 一次拆分，得到 x, y
df = df.with_columns(
    pl.col("column_4")
    .str.split_exact("_", 1)
    .struct.rename_fields(["x", "y"])
    .alias("fields")
).unnest("fields")

df = df.with_columns([
    pl.col("x").cast(pl.Int64),
    pl.col("y").cast(pl.Int64)
])

# 重命名列，便于后续处理
df = df.select([
    pl.col("column_1").alias("chr"),
    pl.col("column_2").alias("start"),
    pl.col("column_3").alias("end"),
    pl.col("x"),
    pl.col("y"),
    pl.col("column_5").alias("counts")
])

# --------------------- 分箱处理 ---------------------
logger.info(f"Binning with binSize={args.binSize}...")
bin_start = time.time()
# 计算 bin 坐标（保留数值，用于后续统计）
df = df.with_columns([
    ((pl.col('x') // args.binSize) * args.binSize).cast(pl.Int64).alias('bin_x'),
    ((pl.col('y') // args.binSize) * args.binSize).cast(pl.Int64).alias('bin_y')
])

# 构造 cell_id（仅用于输出文件）
df = df.with_columns(
    (pl.col('bin_x').cast(pl.Utf8) + "_" + pl.col('bin_y').cast(pl.Utf8)).alias('cell_id')
)

# 准备输出数据（保留原始基因组坐标和新 cell_id）
out_df = df.select(["chr", "start", "end", "cell_id", "counts"])
unique_bins = out_df.select(pl.col("cell_id").n_unique()).item()
logger.info(f"Created {unique_bins} unique bins in {time.time() - bin_start:.2f}s")

# 写入分箱结果
bin_file = os.path.join(args.output, f"{args.sample}_bin{args.binSize}_fragment.tsv")
try:
    out_df.write_csv(bin_file, separator="\t", has_header=False)
    logger.info(f"Binned data written: {bin_file}")
except Exception as e:
    logger.error(f"Failed to write binned data: {e}")
    sys.exit(1)

# --------------------- 统计 bin 内片段数（用于绘图） ---------------------
logger.info("Counting fragments per bin...")
count_start = time.time()
# 直接使用 bin_x, bin_y 进行 groupby，避免后续字符串拆分
counts = df.groupby(["bin_x", "bin_y"]).agg(
    pl.count().alias("n_frags")
)
# 重命名为绘图所需的列名
counts = counts.rename({"bin_x": "x", "bin_y": "y"})
logger.info(f"Counting completed in {time.time() - count_start:.2f}s")

# 释放原始 DataFrame 内存（可选）
del df

# --------------------- 绘图 ---------------------
logger.info("Generating spatial distribution plot...")
plot_start = time.time()
cmap = mcolors.LinearSegmentedColormap.from_list("", [
    "#0E458F", "#0F5298", "#0E6BA8", "#0C86B8", "#3399A1", "#3B9C9C",
    "#B2C061", "#F2CE38", "#F2AB38", "#F2AB38", "#EB7232", "#E65B2E",
    "#E14428", "#DC2E22", "#DB2921", "#CC2623"
])

fig_file = os.path.join(args.output, f"{args.sample}_bin{args.binSize}_nFrag_distribution.pdf")
try:
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(counts['x'], counts['y'], c=counts['n_frags'], cmap=cmap, s=3)
    cbar = plt.colorbar(scatter, shrink=0.4)
    plt.gca().set_facecolor('black')
    plt.gca().invert_yaxis()
    plt.title('')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks([])
    plt.yticks([])
    cbar.ax.set_facecolor('white')
    plt.savefig(fig_file, format='pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"Plot saved: {fig_file} ({time.time() - plot_start:.2f}s)")
except Exception as e:
    logger.error(f"Failed to generate plot: {e}")

# --------------------- 压缩与索引 ---------------------
if not args.no_compress:
    compressed_file = bin_file + ".gz"
    logger.info(f"Compressing with bgzip -> {compressed_file}")
    comp_start = time.time()
    try:
        subprocess.run(['bgzip', '-f', bin_file], check=True)
        logger.info(f"Compression completed in {time.time() - comp_start:.2f}s")
    except subprocess.CalledProcessError as e:
        logger.error(f"bgzip failed: {e}")
        sys.exit(1)

    logger.info("Creating tabix index...")
    idx_start = time.time()
    try:
        subprocess.run(['tabix', '-p', 'bed', compressed_file], check=True)
        logger.info(f"Indexing completed in {time.time() - idx_start:.2f}s -> {compressed_file}.tbi")
    except subprocess.CalledProcessError as e:
        logger.error(f"tabix failed: {e}")
        sys.exit(1)

    if not args.keep_intermediate and os.path.exists(bin_file):
        os.remove(bin_file)
        logger.info(f"Removed uncompressed file: {bin_file}")
else:
    logger.info("Skipping compression/indexing (--no-compress enabled)")

total_elapsed = time.time() - start_time
logger.info(f"=== BinSize Analysis Completed in {total_elapsed:.2f}s ===")