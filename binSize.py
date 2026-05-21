#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 2024
@author: Maven
"""

import polars as pl
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import argparse

parser = argparse.ArgumentParser(description='create bin files')
parser.add_argument('--input','-i', type=str, default = None, help = 'input bin1 fragment file')
parser.add_argument('--output', '-o', type=str, default=None,help='output path')
parser.add_argument('--binSize', '-b', type=int, default=None)
parser.add_argument('--sample', '-s', type=str, default=None)

args = parser.parse_args()

df = pl.read_csv(args.input,separator="\t",has_header = False)

df = df.with_columns(
    pl.col("column_4")
    .str.split_exact("_", 1)
    .struct.rename_fields(["x", "y"])
    .alias("fields")
).unnest("fields")

df = df.with_columns(
    pl.col("x").cast(pl.Int64),
    pl.col("y").cast(pl.Int64)
)

df = df.select(["column_1", "column_2","column_3","x","y","column_5"])

df = df.rename({
    "column_1": "chr",
    "column_2": "start",
    "column_3": "end",
    "column_5":"counts"
})

def parse_bin_coor_no_offset(df: pl.DataFrame, bin_size: int) -> pl.DataFrame:
    df = df.with_columns([
        ((pl.col('x') // bin_size) * bin_size).cast(pl.Int64).alias('bin_x'),
        ((pl.col('y') // bin_size) * bin_size).cast(pl.Int64).alias('bin_y')
    ])
    df = df.with_columns(
        (pl.col('bin_x').cast(pl.Utf8) + "_" + pl.col('bin_y').cast(pl.Utf8)).alias('cell_id')
    )
    return df

binSize = int(args.binSize)
bin_df = parse_bin_coor_no_offset(df,binSize)
bin_df = bin_df.select(["chr", "start","end","cell_id","counts"])
binFile = args.output + "/" + str(args.sample) + "_bin" + str(args.binSize) + "_fragment.tsv"

# 修改这一行：include_header → has_header
bin_df.write_csv(binFile, separator="\t", has_header=False)

## plot the number of 
#counts = bin_df.group_by("cell_id").len()
#counts = bin_df.groupby("cell_id").len()
counts = bin_df.groupby("cell_id").count().rename({"count": "len"})
counts = counts.with_columns(
    pl.col("cell_id")
    .str.split_exact("_", 1)
    .struct.rename_fields(["x", "y"])
    .alias("fields")
).unnest("fields")

counts = counts.with_columns(
    pl.col("x").cast(pl.Int64),
    pl.col("y").cast(pl.Int64)
)

cmap = mcolors.LinearSegmentedColormap.from_list("",[
    "#0E458F","#0F5298","#0E6BA8","#0C86B8","#3399A1","#3B9C9C",
    "#B2C061","#F2CE38","#F2AB38","#F2AB38","#EB7232","#E65B2E",
    "#E14428","#DC2E22","#DB2921","#CC2623"
])

figFile = args.output +  "/" + str(args.sample) + "_bin" + str(args.binSize) + "_nFrag_distribution.pdf"
plt.figure(figsize=(8, 6))
scatter = plt.scatter(counts['x'], counts['y'], c=counts['len'], cmap=cmap, s=3)

cbar = plt.colorbar(scatter, shrink=0.4)
plt.gca().set_facecolor('black') 
plt.gca().invert_yaxis()

plt.title('')
plt.xlabel('')
plt.ylabel('')
plt.xticks([]) 
plt.yticks([]) 
cbar.ax.set_facecolor('white')

plt.savefig(figFile, format='pdf', bbox_inches='tight')
plt.show()