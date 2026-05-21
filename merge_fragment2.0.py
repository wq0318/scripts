#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2026-04-21
@author: wq

Version: v2.1

Changelog:
- v1.0: 初始版本，基础合并功能
- v1.1: 性能优化 - 延迟聚合，修复重复赋值bug
- v1.2: 添加完整日志记录功能，支持 --log 参数
- v2.0: 封装完整流程，集成排序/压缩/索引功能，替代 merge.sh
- v2.1: 智能断点续跑、管道排序+压缩、sort参数优化、饱和度结果固化为JSON
"""

import os
import sys
import time
import json
import logging
import subprocess
import shutil
from datetime import datetime
import polars as pl

import argparse

parser = argparse.ArgumentParser(
    description='Merge, sort, compress and index fragment files (replaces merge.sh)',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  # Basic usage (like merge.sh)
  python merge_fragment.py -i file1.tsv.gz file2.tsv.gz -o ./out -s sample1

  # With custom threads for compression
  python merge_fragment.py -i *.tsv.gz -o ./out -s sample1 --threads 16

  # Skip compression (for testing)
  python merge_fragment.py -i file1.tsv.gz -o ./out -s sample1 --no-compress
""")
parser.add_argument('-i', '--input', dest='input_files', nargs='+', required=True,
                    help='Input fragment files (.tsv or .tsv.gz)')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Output directory')
parser.add_argument('-s', '--sample', type=str, required=True,
                    help='Sample name prefix for output files')
parser.add_argument('-l', '--log', type=str, default=None,
                    help='Log file path (optional)')
parser.add_argument('-t', '--threads', type=int, default=10,
                    help='Threads for sort, bgzip compression (default: 10)')
parser.add_argument('--no-compress', action='store_true',
                    help='Skip compression and indexing (for testing)')
parser.add_argument('--keep-intermediate', action='store_true',
                    help='Keep intermediate unsorted file after completion')
parser.add_argument('--sort-memory', type=str, default='80G',
                    help='Memory limit for sort command (e.g. 80G, 100G)')
parser.add_argument('--sort-tmpdir', type=str, default='/tmp',
                    help='Temporary directory for sort (should be local fast disk)')
args = parser.parse_args()

# 配置日志
logger = logging.getLogger('merge_fragment')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

if args.log:
    log_path = args.log if os.path.isabs(args.log) else os.path.join(args.output, args.log)
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

logger.info(f"=== Merge Fragment Started at {datetime.now()} ===")
logger.info(f"Sample: {args.sample}")
logger.info(f"Output: {args.output}")
logger.info(f"Input files: {len(args.input_files)} file(s)")

# 参数验证
if args.threads < 1:
    logger.warning(f"Invalid threads {args.threads}, using 1")
    args.threads = 1

def check_dependencies():
    """检查必要的外部工具是否可用"""
    missing = []
    for tool in ['sort', 'bgzip', 'tabix']:
        if not shutil.which(tool):
            missing.append(tool)
    if missing:
        logger.error(f"Missing required tools: {', '.join(missing)}")
        logger.error("Please install: samtools (for bgzip/tabix), coreutils (for sort)")
        sys.exit(1)

if not args.no_compress:
    check_dependencies()

if not os.path.exists(args.output):
    os.makedirs(args.output, exist_ok=True)

# 定义路径
unsorted_file = os.path.join(args.output, f"{args.sample}_fragment.tsv")
stat_file = os.path.join(args.output, f"{args.sample}_frag_stat.json")
sorted_file = os.path.join(args.output, f"{args.sample}_sorted_fragment.tsv")
compressed_file = sorted_file + ".gz"

files = args.input_files
for f in files:
    if not os.path.exists(f):
        logger.error(f"输入文件不存在: {f}")
        sys.exit(1)

# ================== 检查是否可以断点续跑 ==================
need_aggregate = True
if os.path.exists(unsorted_file):
    unsorted_mtime = os.path.getmtime(unsorted_file)
    if all(os.path.getmtime(f) < unsorted_mtime for f in files):
        logger.info("Found valid unsorted checkpoint, skipping aggregation...")
        need_aggregate = False

# ================== Step 1: 聚合（如果需要） ==================
total_start = time.time()
if need_aggregate:
    start_time = time.time()
    dfs = []
    total_rows = 0
    for file in files:
        try:
            file_start = time.time()
            df = pl.read_csv(file, separator="\t", has_header=False)
            rows = df.height
            total_rows += rows
            file_elapsed = time.time() - file_start
            logger.info(f"Read {file}: {rows} rows in {file_elapsed:.2f}s")
            dfs.append(df)
        except (pl.exceptions.PolarsError, IOError, OSError) as e:
            logger.error(f"Error reading {file}: {e}")
            sys.exit(1)

    read_elapsed = time.time() - start_time
    logger.info(f"Total read time: {read_elapsed:.2f}s for {total_rows} rows")

    if not dfs:
        logger.error("没有成功读取任何文件")
        sys.exit(1)

    logger.info(f"Merging {len(dfs)} dataframes, total {total_rows} rows...")
    merge_start = time.time()
    all_df = pl.concat(dfs)
    merge_elapsed = time.time() - merge_start
    logger.info(f"Merge completed in {merge_elapsed:.2f}s")

    logger.info("Aggregating data...")
    agg_start = time.time()
    unique_before = all_df.height
    all_df = all_df.groupby(["column_1", "column_2", "column_3", "column_4"]).agg(
        pl.col("column_5").sum()
    )
    unique_after = all_df.height
    agg_elapsed = time.time() - agg_start
    logger.info(f"Aggregation completed in {agg_elapsed:.2f}s: {unique_before} -> {unique_after} unique rows")

    # 写出 unsorted 检查点
    try:
        write_start = time.time()
        all_df.write_csv(unsorted_file, separator="\t", has_header=False)
        write_elapsed = time.time() - write_start
        logger.info(f"Unsorted fragment file written: {unsorted_file} ({all_df.height} rows, {write_elapsed:.2f}s)")
    except (IOError, OSError) as e:
        logger.error(f"Error writing output file: {e}")
        sys.exit(1)

    # 计算并保存饱和度统计（JSON格式，方便断点续跑）
    total_frag_count = all_df.select(pl.col("column_5").sum()).item()
    total_unique_frag_count = all_df.height
    saturation = (1 - total_unique_frag_count / total_frag_count) if total_frag_count > 0 else 0.0
    stats = {
        "total_frag_count": total_frag_count,
        "total_unique_frag_count": total_unique_frag_count,
        "saturation": saturation
    }
    with open(stat_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Statistics saved to: {stat_file}")
else:
    # 读取已存在的统计信息
    if os.path.exists(stat_file):
        with open(stat_file) as f:
            stats = json.load(f)
        logger.info(f"Loaded existing statistics: total_frags={stats['total_frag_count']}, "
                     f"unique_frags={stats['total_unique_frag_count']}, saturation={stats['saturation']:.4f}")
    else:
        logger.warning("Stat file not found, will compute from unsorted file (slow)...")
        # 备用：从 unsorted 重新计算
        all_df = pl.read_csv(unsorted_file, separator="\t", has_header=False)
        total_frag_count = all_df.select(pl.col("column_5").sum()).item()
        total_unique_frag_count = all_df.height
        saturation = (1 - total_unique_frag_count / total_frag_count) if total_frag_count > 0 else 0.0
        stats = {
            "total_frag_count": total_frag_count,
            "total_unique_frag_count": total_unique_frag_count,
            "saturation": saturation
        }
        with open(stat_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Statistics computed and saved to: {stat_file}")

# ================== Step 2: 排序（+ 可选压缩） ==================
if not args.no_compress:
    # 使用管道：sort -> bgzip，直接写出 .gz，不产生中间的 sorted.tsv
    logger.info("Step 2: Sorting and compressing via pipe...")
    sort_start = time.time()
    sort_cmd = [
        'sort',
        '-t', '\t',
        '--parallel', str(args.threads),
        '-S', args.sort_memory,
        '-T', args.sort_tmpdir,
        '-k', '1,1V', '-k', '2,2n', '-k', '3,3n',
        unsorted_file
    ]
    bgzip_cmd = ['bgzip', '-@', str(args.threads), '-c']

    try:
        with open(compressed_file, 'wb') as out_f:
            sort_proc = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bgzip_proc = subprocess.Popen(bgzip_cmd, stdin=sort_proc.stdout, stdout=out_f, stderr=subprocess.PIPE)
            sort_proc.stdout.close()  # 允许 SIGPIPE 传递

            sort_ret = sort_proc.wait()
            if sort_ret != 0:
                stderr_msg = sort_proc.stderr.read().decode()
                raise subprocess.CalledProcessError(sort_ret, sort_cmd, stderr=stderr_msg)

            bgzip_ret = bgzip_proc.wait()
            if bgzip_ret != 0:
                stderr_msg = bgzip_proc.stderr.read().decode()
                raise subprocess.CalledProcessError(bgzip_ret, bgzip_cmd, stderr=stderr_msg)

        sort_elapsed = time.time() - sort_start
        logger.info(f"Sort + compress completed in {sort_elapsed:.2f}s -> {compressed_file}")

        # Step 3: 创建 tabix 索引
        logger.info("Step 3: Creating tabix index...")
        idx_start = time.time()
        subprocess.run(['tabix', '-p', 'bed', compressed_file], check=True)
        logger.info(f"Indexing completed in {time.time() - idx_start:.2f}s -> {compressed_file}.tbi")

        # 清理 unsorted 检查点（除非用户要求保留）
        if not args.keep_intermediate and os.path.exists(unsorted_file):
            os.remove(unsorted_file)
            logger.info(f"Removed checkpoint: {unsorted_file}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Pipeline error: {e}")
        if os.path.exists(compressed_file):
            os.remove(compressed_file)
        logger.info("Checkpoint file still available, you can retry the sort+compress step.")
        sys.exit(1)

else:
    # --no-compress：仍然排序，但不压缩
    logger.info("Step 2: Sorting without compression...")
    sort_start = time.time()
    sort_cmd = [
        'sort',
        '-t', '\t',
        '--parallel', str(args.threads),
        '-S', args.sort_memory,
        '-T', args.sort_tmpdir,
        '-k', '1,1V', '-k', '2,2n', '-k', '3,3n',
        unsorted_file
    ]
    try:
        with open(sorted_file, 'w') as out_f:
            subprocess.run(sort_cmd, stdout=out_f, check=True)
        sort_elapsed = time.time() - sort_start
        logger.info(f"Sorting completed in {sort_elapsed:.2f}s -> {sorted_file}")
        if not args.keep_intermediate and os.path.exists(unsorted_file):
            os.remove(unsorted_file)
            logger.info(f"Removed checkpoint: {unsorted_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Sorting failed: {e}")
        sys.exit(1)

# ================== 输出总结 ==================
logger.info("Output files:")
final_files = [stat_file]
if args.no_compress:
    final_files.append(sorted_file)
else:
    final_files.extend([compressed_file, f"{compressed_file}.tbi"])

for f in final_files:
    if os.path.exists(f):
        size = os.path.getsize(f)
        logger.info(f"  - {f} ({size / (1024**2):.1f} MB)")

total_elapsed = time.time() - total_start
logger.info(f"=== Merge Fragment Pipeline Completed in {total_elapsed:.2f}s ===")