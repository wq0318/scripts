#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import polars as pl
import argparse

parser = argparse.ArgumentParser(description='Read multiple fragment files into a Polars DataFrame')
parser.add_argument('-i', dest='input_files', nargs='*', required=True, help='List of input files')
parser.add_argument('--output', '-o', type=str, default=None, help='output path')
parser.add_argument('--sample', '-s', type=str, default=None)
args = parser.parse_args()

files = args.input_files
all_df = pl.read_csv(files[0], separator="\t", has_header=False)

for file in files[1:]:
    try:
        df = pl.read_csv(file, separator="\t", has_header=False)
        print(f"Successfully read {file} with {df.shape[0]} rows and {df.shape[1]} columns")
        all_df = pl.concat([all_df, df])
        # 改动：group_by → groupby
        all_df = (
            all_df.groupby(["column_1","column_2","column_3","column_4"])
            .agg([
                pl.col("column_5").sum()
            ])
        )
    except Exception as e:
        print(f"Error reading {file}: {e}")

File = args.output + "/" + str(args.sample) + "_fragment.tsv"
# 改动：include_header=False → has_header=False
all_df.write_csv(File, separator="\t", has_header=False)

# Calculate saturation
total_frag_count = all_df.select(pl.col("column_5").sum()).item()
total_unique_frag_count = all_df.height
saturation = (1 - total_unique_frag_count / total_frag_count)
data = {
    "total_frag_count": total_frag_count,
    "total_unique_frag_count": total_unique_frag_count,
    "saturation": saturation
}
fn = args.output + "/" + str(args.sample) + "_frag_stat.tsv"
with open(fn, 'w') as f:
    for k, v in data.items():
        f.write(k + ':' + str(v) + '\n')