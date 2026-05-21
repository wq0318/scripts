#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
import os
from collections import defaultdict, deque

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(
        description="StereoATAC tissuecut (data-driven, no histology image)."
    )
    p.add_argument("--fragments", required=True, help="Fragment file: chr start end x_y pcr_dup")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--target-bin", type=int, default=1, help="Final output bin level (default: 1)")
    p.add_argument(
        "--guide-bin",
        type=int,
        default=100,
        help="Optional coarse bin for guiding target-bin mask. <=0 disables guide.",
    )
    p.add_argument("--pixel-um", type=float, default=0.5, help="Physical size (um) of bin1 pixel")
    p.add_argument("--smooth-um", type=float, default=50.0, help="Smoothing radius in um")
    p.add_argument(
        "--min-component-um2",
        type=float,
        default=40000.0,
        help="Minimum component area in um^2",
    )
    p.add_argument(
        "--use-dedup-weight",
        action="store_true",
        help="Use 1/max(pcr_dup,1) as weight to suppress PCR-heavy noise",
    )
    p.add_argument(
        "--iterative-otsu",
        action="store_true",
        help="Use T1 then Otsu again on values>=T1 to get T2 (recommended for 3-peak distributions)",
    )
    p.add_argument(
        "--tss-file",
        default="",
        help="Optional: x_y<TAB>tss or x<TAB>y<TAB>tss. Used only as negative penalty.",
    )
    p.add_argument(
        "--tss-high-z",
        type=float,
        default=1.5,
        help="Drop pixel if TSS z > this and FRAGS z < 0 inside current mask",
    )
    p.add_argument(
        "--write-mask-tsv",
        action="store_true",
        help="Write dense x,y,keep table (very large for big chips; default off)",
    )
    return p.parse_args()


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_xy(xy):
    parts = xy.split("_")
    if len(parts) != 2:
        raise ValueError
    return int(parts[0]), int(parts[1])


def stream_fragments_to_bins(path, bins, use_dedup_weight):
    """Stream fragments and aggregate into aligned bin coordinates.

    Key change vs old version:
      - Use aligned coordinates ((x//b)*b, (y//b)*b) instead of index (x//b, y//b).
      - Use pandas vectorized parsing + groupby for speed.
      - Read in chunks to avoid loading the full fragment file into memory.
    """
    agg = {b: defaultdict(float) for b in bins}
    raw_rows = 0
    bad_rows = 0

    try:
        import pandas as pd
    except Exception as e:
        raise RuntimeError(f"pandas is required for this script version. Import failed: {e}")

    compression = "gzip" if path.endswith(".gz") else None

    # Use chunks to reduce peak memory. Tune if needed.
    chunksize = 1_000_000
    col_names = ["chrom", "start", "end", "xy", "pcr_dup"]

    for chunk in pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=col_names,
        usecols=[3, 4],
        dtype={"xy": "string", "pcr_dup": "float32"},
        engine="c",
        compression=compression,
        chunksize=chunksize,
    ):
        raw_rows += int(chunk.shape[0])

        coords = chunk["xy"].str.extract(r"(\d+)[_](\d+)")
        bad = coords.isna().any(axis=1)
        if bool(bad.any()):
            bad_rows += int(bad.sum())
        coords = coords[~bad]
        if coords.shape[0] == 0:
            continue

        x = coords[0].astype("int64")
        y = coords[1].astype("int64")
        dup = chunk.loc[~bad, "pcr_dup"].astype("float64")

        if use_dedup_weight:
            w = 1.0 / np.maximum(dup.to_numpy(), 1.0)
        else:
            w = np.ones(dup.shape[0], dtype=np.float64)

        # Aggregate for each requested bin.
        for b in bins:
            bx = ((x // b) * b).to_numpy(dtype=np.int64)
            by = ((y // b) * b).to_numpy(dtype=np.int64)

            dfb = pd.DataFrame({"bx": bx, "by": by, "w": w})
            g = dfb.groupby(["bx", "by"], sort=False)["w"].sum()
            for (abx, aby), val in g.items():
                agg[b][(int(abx), int(aby))] += float(val)

    return agg, raw_rows, bad_rows


def load_tss(path):
    if not path:
        return {}
    tss = {}
    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            try:
                if len(cols) >= 3 and cols[0].isdigit() and cols[1].isdigit():
                    x, y, v = int(cols[0]), int(cols[1]), float(cols[2])
                elif len(cols) >= 2:
                    x, y = parse_xy(cols[0])
                    v = float(cols[1])
                else:
                    continue
                tss[(x, y)] = v
            except Exception:
                continue
    return tss


def aggregate_tss_to_bin(tss_map, bin_size):
    if not tss_map:
        return {}
    tmp = defaultdict(list)
    for (x, y), v in tss_map.items():
        tmp[(((x // bin_size) * bin_size), ((y // bin_size) * bin_size))].append(v)
    return {k: float(np.median(vs)) for k, vs in tmp.items()}


def build_grid(bcount, bin_size):
    keys = list(bcount.keys())
    xs = np.array([k[0] for k in keys], dtype=np.int32)
    ys = np.array([k[1] for k in keys], dtype=np.int32)
    vals = np.array([bcount[k] for k in keys], dtype=np.float32)
    minx, maxx = int(xs.min()), int(xs.max())
    miny, maxy = int(ys.min()), int(ys.max())

    # xs/ys are aligned coordinates (multiples of bin_size). Build a compact grid.
    w = (maxx - minx) // bin_size + 1
    h = (maxy - miny) // bin_size + 1
    grid = np.zeros((h, w), dtype=np.float32)

    gx = (xs - minx) // bin_size
    gy = (ys - miny) // bin_size
    grid[gy, gx] = vals
    return grid, minx, miny


def box_smooth(arr, radius):
    if radius <= 0:
        return arr.copy()
    r = int(radius)
    pad = np.pad(arr, ((r, r), (r, r)), mode="edge")
    ii = np.pad(pad, ((1, 0), (1, 0)), mode="constant", constant_values=0)
    ii = ii.cumsum(axis=0).cumsum(axis=1)
    k = float((2 * r + 1) ** 2)
    s = ii[2 * r + 1 :, 2 * r + 1 :] - ii[: -2 * r - 1, 2 * r + 1 :] - ii[2 * r + 1 :, : -2 * r - 1] + ii[: -2 * r - 1, : -2 * r - 1]
    return s / k


def otsu_threshold(values):
    v = values[np.isfinite(values)]
    v = v[v > 0]
    if v.size == 0:
        return 0.0
    hist, edges = np.histogram(v, bins=256)
    centers = (edges[:-1] + edges[1:]) / 2.0
    w_total = hist.sum()
    sum_total = (hist * centers).sum()
    sum_b = 0.0
    w_b = 0.0
    best = float(np.median(v))
    best_var = -1.0
    for i in range(hist.size):
        w_b += hist[i]
        if w_b == 0:
            continue
        w_f = w_total - w_b
        if w_f == 0:
            break
        sum_b += hist[i] * centers[i]
        m_b = sum_b / w_b
        m_f = (sum_total - sum_b) / w_f
        var_between = w_b * w_f * (m_b - m_f) ** 2
        if var_between > best_var:
            best_var = var_between
            best = centers[i]
    return float(best)


def robust_z(v):
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    if mad < 1e-8:
        return np.zeros_like(v)
    return 0.6745 * (v - med) / mad


def connected_components(mask):
    h, w = mask.shape
    visited = np.zeros_like(mask, dtype=np.uint8)
    comps = []
    for y in range(h):
        for x in range(w):
            if mask[y, x] == 0 or visited[y, x]:
                continue
            q = deque([(y, x)])
            visited[y, x] = 1
            pts = []
            while q:
                cy, cx = q.popleft()
                pts.append((cy, cx))
                if cy > 0 and mask[cy - 1, cx] and not visited[cy - 1, cx]:
                    visited[cy - 1, cx] = 1
                    q.append((cy - 1, cx))
                if cy < h - 1 and mask[cy + 1, cx] and not visited[cy + 1, cx]:
                    visited[cy + 1, cx] = 1
                    q.append((cy + 1, cx))
                if cx > 0 and mask[cy, cx - 1] and not visited[cy, cx - 1]:
                    visited[cy, cx - 1] = 1
                    q.append((cy, cx - 1))
                if cx < w - 1 and mask[cy, cx + 1] and not visited[cy, cx + 1]:
                    visited[cy, cx + 1] = 1
                    q.append((cy, cx + 1))
            comps.append(pts)
    return comps


def remove_small_components(mask, min_size):
    out = mask.copy()
    for comp in connected_components(mask):
        if len(comp) < min_size:
            for y, x in comp:
                out[y, x] = 0
    return out


def fill_holes(mask):
    h, w = mask.shape
    inv = (mask == 0).astype(np.uint8)
    visited = np.zeros_like(inv, dtype=np.uint8)
    q = deque()
    for x in range(w):
        if inv[0, x]:
            q.append((0, x))
        if inv[h - 1, x]:
            q.append((h - 1, x))
    for y in range(h):
        if inv[y, 0]:
            q.append((y, 0))
        if inv[y, w - 1]:
            q.append((y, w - 1))
    while q:
        y, x = q.popleft()
        if visited[y, x]:
            continue
        visited[y, x] = 1
        if y > 0 and inv[y - 1, x] and not visited[y - 1, x]:
            q.append((y - 1, x))
        if y < h - 1 and inv[y + 1, x] and not visited[y + 1, x]:
            q.append((y + 1, x))
        if x > 0 and inv[y, x - 1] and not visited[y, x - 1]:
            q.append((y, x - 1))
        if x < w - 1 and inv[y, x + 1] and not visited[y, x + 1]:
            q.append((y, x + 1))
    holes = (inv == 1) & (visited == 0)
    out = mask.copy()
    out[holes] = 1
    return out


def compute_threshold(smooth, iterative_otsu):
    t1 = otsu_threshold(smooth)
    t2 = t1
    if iterative_otsu:
        sub = smooth[smooth >= t1]
        if sub.size > 0:
            t2 = otsu_threshold(sub)
    return t1, t2


def build_tss_grid(btss, h, w, minx, miny, bin_size):
    if not btss:
        return None
    out = np.full((h, w), np.nan, dtype=np.float32)
    for (bx, by), v in btss.items():
        # bx/by are aligned coordinates (multiples of bin_size). Map to compact indices.
        gx = (int(bx) - int(minx)) // int(bin_size)
        gy = (int(by) - int(miny)) // int(bin_size)
        if 0 <= gy < h and 0 <= gx < w:
            out[gy, gx] = float(v)
    return out


def apply_tss_penalty(mask, log_frags, tss_grid, tss_high_z):
    if tss_grid is None:
        return mask
    idx = np.where(mask > 0)
    if idx[0].size < 10:
        return mask
    fr = log_frags[idx]
    ts = tss_grid[idx]
    good = np.isfinite(ts)
    if good.sum() < 10:
        return mask
    fz = robust_z(fr[good])
    tz = robust_z(ts[good])
    drop = (tz > tss_high_z) & (fz < 0)
    out = mask.copy()
    oy = idx[0][good][drop]
    ox = idx[1][good][drop]
    out[oy, ox] = 0
    return out


def build_gate_from_guide(
    target_shape,
    target_minx,
    target_miny,
    target_bin,
    guide_mask,
    guide_minx,
    guide_miny,
    guide_bin,
):
    h, w = target_shape
    yy, xx = np.indices((h, w), dtype=np.int64)

    # target_minx/miny are aligned coordinates; each step is target_bin.
    bx = int(target_minx) + xx * int(target_bin)
    by = int(target_miny) + yy * int(target_bin)

    # Map each target aligned coord to guide aligned coord by flooring to guide_bin.
    gbx = (bx // int(guide_bin)) * int(guide_bin)
    gby = (by // int(guide_bin)) * int(guide_bin)

    # Convert aligned coords to guide grid indices.
    gx = (gbx - int(guide_minx)) // int(guide_bin)
    gy = (gby - int(guide_miny)) // int(guide_bin)
    gate = np.zeros((h, w), dtype=np.uint8)
    valid = (
        (gy >= 0)
        & (gy < guide_mask.shape[0])
        & (gx >= 0)
        & (gx < guide_mask.shape[1])
    )
    gate[valid] = guide_mask[gy[valid], gx[valid]]
    return gate


def save_tif(mask_uint8, out_tif):
    # mask_uint8 should be 0/1. Write as 0/255.
    img = (mask_uint8 * 255).astype(np.uint8)
    try:
        import tifffile  # type: ignore

        tifffile.imwrite(out_tif, img, photometric="minisblack")
        return "tifffile"
    except Exception:
        pass
    try:
        from PIL import Image  # type: ignore

        Image.fromarray(img, mode="L").save(out_tif, format="TIFF")
        return "PIL"
    except Exception as e:
        raise RuntimeError(
            "Failed to write TIFF. Please install tifffile or pillow. "
            f"Last error: {e}"
        )


def run_one_bin(
    bcount, btss, bin_size, pixel_um, smooth_um, min_component_um2, iterative_otsu, tss_high_z, guide_gate=None
):
    grid, minx, miny = build_grid(bcount, bin_size)
    log_grid = np.log1p(grid)
    pixel_um_now = pixel_um * bin_size
    smooth_r = max(0, int(round(smooth_um / pixel_um_now)))
    smooth = box_smooth(log_grid, smooth_r)

    t1, t2 = compute_threshold(smooth, iterative_otsu)
    mask = (smooth >= t2).astype(np.uint8)
    if guide_gate is not None:
        mask = (mask & guide_gate).astype(np.uint8)

    min_component_px = max(1, int(round(min_component_um2 / (pixel_um_now * pixel_um_now))))
    mask = remove_small_components(mask, min_component_px)
    mask = fill_holes(mask)

    tss_grid = build_tss_grid(btss, mask.shape[0], mask.shape[1], minx, miny, bin_size)
    mask = apply_tss_penalty(mask, log_grid, tss_grid, tss_high_z)
    mask = remove_small_components(mask, min_component_px)

    keep = int(mask.sum())
    total = int(mask.size)
    return {
        "bin_size": bin_size,
        "minx": minx,
        "miny": miny,
        "grid": grid,
        "log_grid": log_grid,
        "smooth": smooth,
        "threshold_t1": t1,
        "threshold_t2": t2,
        "smooth_radius_px": smooth_r,
        "min_component_px": min_component_px,
        "mask": mask,
        "tss_grid": tss_grid,
        "keep_pixels": keep,
        "total_pixels": total,
        "keep_ratio": (keep / total if total else 0.0),
    }


def save_outputs(res, outdir, write_mask_tsv):
    b = res["bin_size"]
    prefix = os.path.join(outdir, f"bin{b}")
    os.makedirs(prefix, exist_ok=True)
    np.save(os.path.join(prefix, "frags_grid.npy"), res["grid"])
    np.save(os.path.join(prefix, "logfrags_grid.npy"), res["log_grid"])
    np.save(os.path.join(prefix, "mask.npy"), res["mask"].astype(np.uint8))
    if res["tss_grid"] is not None:
        np.save(os.path.join(prefix, "tss_grid.npy"), res["tss_grid"])
    if write_mask_tsv:
        with open(os.path.join(prefix, "mask.tsv"), "w", encoding="utf-8") as fo:
            fo.write("x\ty\tkeep\tfrags\n")
            h, w = res["mask"].shape
            for gy in range(h):
                for gx in range(w):
                    # x/y are aligned coordinates in the same unit as fragments.
                    x = int(res["minx"]) + gx * int(res["bin_size"])
                    y = int(res["miny"]) + gy * int(res["bin_size"])
                    fo.write(f"{x}\t{y}\t{int(res['mask'][gy, gx])}\t{float(res['grid'][gy, gx]):.6f}\n")
    summ = {
        "bin_size": b,
        "threshold_t1": res["threshold_t1"],
        "threshold_t2": res["threshold_t2"],
        "smooth_radius_px": res["smooth_radius_px"],
        "min_component_px": res["min_component_px"],
        "keep_pixels": res["keep_pixels"],
        "total_pixels": res["total_pixels"],
        "keep_ratio": res["keep_ratio"],
        "grid_origin_minx": res["minx"],
        "grid_origin_miny": res["miny"],
    }
    with open(os.path.join(prefix, "summary.json"), "w", encoding="utf-8") as fo:
        json.dump(summ, fo, ensure_ascii=False, indent=2)
    return summ


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    bins = [args.target_bin]
    if args.guide_bin > 0 and args.guide_bin != args.target_bin:
        bins.append(args.guide_bin)
        bins = sorted(set(bins))

    binned_counts, raw_rows, bad_rows = stream_fragments_to_bins(
        args.fragments, bins, args.use_dedup_weight
    )
    tss_map = load_tss(args.tss_file)
    btss = {b: aggregate_tss_to_bin(tss_map, b) for b in bins}

    guide_res = None
    if args.guide_bin > 0 and args.guide_bin in bins and args.guide_bin != args.target_bin:
        guide_res = run_one_bin(
            bcount=binned_counts[args.guide_bin],
            btss=btss[args.guide_bin],
            bin_size=args.guide_bin,
            pixel_um=args.pixel_um,
            smooth_um=args.smooth_um,
            min_component_um2=args.min_component_um2,
            iterative_otsu=args.iterative_otsu,
            tss_high_z=args.tss_high_z,
            guide_gate=None,
        )
        save_outputs(guide_res, args.outdir, write_mask_tsv=False)

    target_res = run_one_bin(
        bcount=binned_counts[args.target_bin],
        btss=btss[args.target_bin],
        bin_size=args.target_bin,
        pixel_um=args.pixel_um,
        smooth_um=args.smooth_um,
        min_component_um2=args.min_component_um2,
        iterative_otsu=args.iterative_otsu,
        tss_high_z=args.tss_high_z,
        guide_gate=None,
    )

    # Apply guide gate after target mask generation to ensure coordinate consistency.
    if guide_res is not None:
        gate = build_gate_from_guide(
            target_shape=target_res["mask"].shape,
            target_minx=target_res["minx"],
            target_miny=target_res["miny"],
            target_bin=args.target_bin,
            guide_mask=guide_res["mask"],
            guide_minx=guide_res["minx"],
            guide_miny=guide_res["miny"],
            guide_bin=args.guide_bin,
        )
        target_res["mask"] = (target_res["mask"] & gate).astype(np.uint8)
        target_res["keep_pixels"] = int(target_res["mask"].sum())
        target_res["keep_ratio"] = (
            target_res["keep_pixels"] / target_res["total_pixels"]
            if target_res["total_pixels"]
            else 0.0
        )

    target_summ = save_outputs(target_res, args.outdir, args.write_mask_tsv)

    tif_path = os.path.join(args.outdir, f"bin{args.target_bin}", "mask.tif")
    tif_writer = save_tif(target_res["mask"], tif_path)

    all_summ = {
        "fragments_file": args.fragments,
        "raw_fragment_rows": raw_rows,
        "bad_rows_skipped": bad_rows,
        "target_bin": args.target_bin,
        "guide_bin": args.guide_bin if args.guide_bin > 0 else None,
        "use_dedup_weight": bool(args.use_dedup_weight),
        "iterative_otsu": bool(args.iterative_otsu),
        "tss_file": args.tss_file if args.tss_file else None,
        "tif_writer": tif_writer,
        "target_summary": target_summ,
    }
    if guide_res is not None:
        all_summ["guide_summary"] = {
            "bin_size": guide_res["bin_size"],
            "threshold_t1": guide_res["threshold_t1"],
            "threshold_t2": guide_res["threshold_t2"],
            "keep_ratio": guide_res["keep_ratio"],
        }

    with open(os.path.join(args.outdir, "summary_all.json"), "w", encoding="utf-8") as fo:
        json.dump(all_summ, fo, ensure_ascii=False, indent=2)
    print(json.dumps(all_summ, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
