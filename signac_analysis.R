#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(Signac)
    library(Seurat)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
    library(Matrix)
    library(tidyr)
    library(dplyr)
    library(SeuratObject)
    library(future)
    library(data.table)
    library(ArchR)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
    stop("Usage: Rscript signac.R <fragments.tsv.gz> <sample_name> <output_dir> <bin_tissue_barcode.csv> [threads] [memory_gb] [genome]")
}

fragpath <- args[1]
sample_name <- args[2]
outdir <- args[3]
barcode_csv <- args[4]

parse_positive_int <- function(value, default_value) {
    parsed <- suppressWarnings(as.integer(value))
    if (is.na(parsed) || parsed < 1) {
        return(default_value)
    }
    parsed
}

WORKERS <- if(length(args) >= 5) parse_positive_int(args[5], 8) else 8
MEM_GB <- if(length(args) >= 6) parse_positive_int(args[6], 64) else 64
MAX_MEM <- MEM_GB * 1024^3
GENOME <- if(length(args) >= 7 && nzchar(args[7])) args[7] else "mm10"
if(!(GENOME %in% c("mm10", "hg38"))) {
    stop(sprintf("Unsupported genome '%s'; expected one of: mm10, hg38", GENOME))
}
plan("multisession", workers = WORKERS)
options(future.globals.maxSize = MAX_MEM)

pal <- c('#D51F26', '#272E6A', '#208A42', '#89288F', '#F47D2B', 
            '#FEE500', '#8A9FD1', '#C06CAB', '#E6C2DC', '#90D5E4', 
            '#89C75F', '#F37B7D', '#9983BD', '#D24B27', '#3BBCA8', 
            '#6E4B9E', '#0C727C', '#7E1416', '#D8A767', '#3D3D3D', 
            '#371377', '#7700FF', '#9E0142', '#FF0080', '#DC494C', 
            '#F88D51', '#FAD510', '#FFFF5F', '#88CFA4', '#238B45', 
            '#02401B', '#0AD7D3', '#046C9A', '#A2A475', 'grey35')

message(sprintf("[Signac] Processing: %s", sample_name))
message(sprintf("[Signac] Fragment file: %s", fragpath))
message(sprintf("[Signac] Tissue barcode file: %s", barcode_csv))
message(sprintf("[Signac] Threads: %d", WORKERS))
message(sprintf("[Signac] Future memory: %d GB", MEM_GB))

if(!file.exists(fragpath)) {
    stop(sprintf("Fragment file not found: %s", fragpath))
}
if(!file.exists(barcode_csv)) {
    stop(sprintf("Barcode CSV not found: %s", barcode_csv))
}

barcode_table <- fread(barcode_csv)
if(!all(c("CB", "nFrags") %in% colnames(barcode_table))) {
    stop("Barcode CSV must contain columns: CB,nFrags")
}
barcode_table <- barcode_table[!is.na(CB) & CB != ""]
barcode_table[, nFrags := as.numeric(nFrags)]
barcode_table <- barcode_table[!is.na(nFrags)]
barcodes <- unique(barcode_table$CB)


message(sprintf("[Signac] Barcodes: %d", length(barcodes)))

index_file <- paste0(fragpath, ".tbi")
if (!file.exists(index_file)) {
    message("[Signac] Fragment index not found, creating...")
    IndexFragments(fragment.path = fragpath, verbose = TRUE)
    message("[Signac] Index created")
}

frags <- CreateFragmentObject(path = fragpath, cells = barcodes, validate.fragments = FALSE)
peaks <- CallPeaks(frags)
saveRDS(peaks, file.path(outdir, paste0(sample_name, "_peaks.rds")))

counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)

chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragpath,
    min.cells = 0,
    min.features = 0
)
pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sample_name
)

if(GENOME == "mm10") {
    library(EnsDb.Mmusculus.v79)
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
} else if(GENOME == "hg38") {
    library(EnsDb.Hsapiens.v86)
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
} else {
    stop(sprintf("Unsupported genome: %s", GENOME))
}
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- GENOME
Annotation(pbmc) <- annotations

cell_barcodes <- rownames(pbmc@meta.data)
pbmc$coor_x <- as.numeric(sapply(cell_barcodes, function(x) strsplit(x, "_")[[1]][1]))
pbmc$coor_y <- as.numeric(sapply(cell_barcodes, function(x) strsplit(x, "_")[[1]][2]))


# nFrags already computed upstream (compute_tissue_barcodes.py, chrM excluded).
# We pass it through metadata rather than recomputing it from the fragment file.
pbmc$nFrags <- barcode_table$nFrags[match(Cells(pbmc), barcode_table$CB)]
pbmc$nFrags[is.na(pbmc$nFrags)] <- 0

pbmc <- subset(pbmc, subset = nCount_ATAC > 10)
pbmc <- FRiP(object = pbmc, assay = 'ATAC', total.fragments = 'nFrags')

# The chromatin assay holds references to the (heavy) counts matrix + fragment
# object; once they're inside pbmc we no longer need the bare versions.
rm(counts, chrom_assay, frags); gc(verbose = FALSE)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 20)

n_cells <- ncol(pbmc)
n_features <- nrow(pbmc)
max_svd_dims <- min(n_cells, n_features) - 2
if (max_svd_dims < 8) {
    message(sprintf("[Signac] Warning: Low dimensionality (cells=%d, features=%d), adjusting LSI dims", n_cells, n_features))
    n_svd <- max(5, max_svd_dims)
    dims_to_use <- 2:max(2, max_svd_dims)
} else {
    n_svd <- 30
    dims_to_use <- 2:10
}
message(sprintf("[Signac] SVD dims: %d, Using dims: %s", n_svd, paste(dims_to_use, collapse=":")))

pbmc <- RunSVD(pbmc, n = n_svd)
pbmc <- RunUMAP(pbmc, dims = dims_to_use, reduction = 'lsi')
pbmc <- FindNeighbors(object = pbmc, dims = dims_to_use, reduction = 'lsi')

message("[Signac] Clustering at resolution 0.8")

default_assay <- DefaultAssay(pbmc)
message(sprintf("[Signac] Default assay: %s", default_assay))

res <- 0.8
res_str <- gsub("\\.", "_", as.character(res))
cluster_col <- paste0(default_assay, "_snn_res.", res)

tryCatch({
    pbmc <- FindClusters(object = pbmc, algorithm = 3, resolution = res, verbose = FALSE)
    pbmc[[paste0("clusters_res_", res_str)]] <- Idents(pbmc)

    if (!cluster_col %in% colnames(pbmc@meta.data)) {
        alt_col <- "seurat_clusters"
        if (alt_col %in% colnames(pbmc@meta.data)) {
            cluster_col <- alt_col
            message(sprintf("[Signac] Warning: Using fallback column '%s' for res=%s", cluster_col, res))
        } else {
            stop(sprintf("Cluster column '%s' not found in metadata", cluster_col))
        }
    }

    Idents(pbmc) <- cluster_col

    n_clusters <- length(unique(Idents(pbmc)))
    # Always use colorRampPalette so the palette gracefully extends past the
    # hand-picked colours when cluster count grows (no implicit upper bound).
    cluster_pal <- colorRampPalette(pal)(n_clusters)

    p_umap <- DimPlot(pbmc, label = TRUE, cols = cluster_pal, pt.size = 0.1, reduction = "umap")
    p_spatial <- ggplot() +
        geom_tile(data = pbmc@meta.data,
                  aes(x = coor_x, y = coor_y, fill = .data[[cluster_col]])) +
        scale_fill_manual(values = cluster_pal) +
        theme_void() +
        coord_fixed() +
        theme(text = element_text(size = 16))

    cluster_plot <- p_umap + p_spatial
    plot_file <- file.path(outdir, paste0(sample_name, "_clusters_res", res_str, ".svg"))
    ggsave(plot_file, cluster_plot, width = 12, height = 6, device = "svg")
    message(sprintf("[Signac] Cluster plot (res=%s): %s", res, plot_file))

}, error = function(e) {
    message(sprintf("[Signac] Warning: Failed to generate cluster plot for res=%s: %s", res, e$message))
    blank_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = sprintf("Plot failed\nres=%s", res)) +
        theme_void()
    plot_file <- file.path(outdir, paste0(sample_name, "_clusters_res", res_str, ".svg"))
    ggsave(plot_file, blank_plot, width = 6, height = 6, device = "svg")
})


ns_result <- tryCatch({
    NucleosomeSignal(object = pbmc)
}, error = function(e) {
    message("[Signac] Warning: NucleosomeSignal failed: ", e$message)
    pbmc$nucleosome_signal <- NA
    pbmc
})
if (!is.null(ns_result)) pbmc <- ns_result

tss_result <- tryCatch({
    TSSEnrichment(object = pbmc)
}, error = function(e) {
    message("[Signac] Warning: TSSEnrichment failed (possibly due to sparse test data): ", e$message)
    pbmc$TSS.enrichment <- 0
    pbmc
})
if (!is.null(tss_result)) pbmc <- tss_result

bl_result <- tryCatch({
    blacklist_regions <- if (GENOME == "hg38") Signac::blacklist_hg38_unified else Signac::blacklist_mm10
    pbmc$blacklist_ratio <- FractionCountsInRegion(object = pbmc, regions = blacklist_regions)
    pbmc
}, error = function(e) {
    message("[Signac] Warning: Blacklist ratio calculation failed: ", e$message)
    pbmc$blacklist_ratio <- 0
    pbmc
})
if (!is.null(bl_result)) pbmc <- bl_result


p1 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=nFrags),color='black',fill='#5A8100',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=nFrags),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(median(pbmc@meta.data$nFrags))

p2 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=TSS.enrichment),color='black',fill='#178CA4',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=TSS.enrichment),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(round(median(pbmc@meta.data$TSS.enrichment), 2))

p3 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=FRiP*100),color='black',fill='#FFB400',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=FRiP*100),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(paste0(round(median(pbmc@meta.data$FRiP)*100, 1), "%"))

p4 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=blacklist_ratio*100),color='black',fill='#FF6C02',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=blacklist_ratio*100),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(paste0(round(median(pbmc@meta.data$blacklist_ratio)*100, 2), "%"))

qc_violin_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
tryCatch({
    ggsave(file.path(outdir, paste0(sample_name, "_qc_violin.svg")), qc_violin_plot, width = 16, height = 4, device = "svg")
    message("[Signac] QC violin plot saved")
}, error = function(e) {
    message("[Signac] Warning: QC violin plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="QC plot failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_qc_violin.svg")), blank, width = 16, height = 4, device = "svg")
})

tryCatch({
    # ggPoint() opens a graphics device as a side effect; pdf(NULL) + dev.off()
    # absorbs the implicit PDF so the working directory stays clean.
    pdf(NULL)
    p5 <- ggPoint(
        x = log10(pbmc@meta.data$nFrags),
        y = pbmc@meta.data$TSS.enrichment,
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(log10(pbmc@meta.data$nFrags), probs = 0.99, na.rm = TRUE)),
        ylim = c(0, quantile(pbmc@meta.data$TSS.enrichment, probs = 0.99, na.rm = TRUE))
    ) +
    theme(
        legend.position = "bottom",
        legend.justification = "center"
    )
    dev.off()
    ggsave(file.path(outdir, paste0(sample_name, "_tss_scatter.svg")), p5, width = 6, height = 6, device = "svg")
    message("[Signac] TSS scatter plot saved")
}, error = function(e) {
    message("[Signac] Warning: TSS scatter plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="TSS scatter failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_tss_scatter.svg")), blank, width = 6, height = 6, device = "svg")
})

# Fragment size distribution: data is consumed by render_final_report.py via
# ECharts (responsive to container size), so we emit CSV + JS data files only —
# no static SVG (the previous SVG distorted at non-default aspect ratios).
threads <- WORKERS
cmd <- paste0("pigz -dc -p ", threads, " ", fragpath,
  " | awk '$1 != \"chrM\" && $3-$2 <= 750 && $3-$2 > 0 { counts[$3-$2]++ } END { for (i=1; i<=750; i++) print i, counts[i]+0 }'")
size_df <- fread(cmd = cmd, col.names = c("size", "count"))
size_df[, percentage := (count / sum(count)) * 100]
write.csv(size_df, file.path(outdir, paste0(sample_name, "_fragment_size.csv")), row.names = FALSE)
nfrags_total <- sum(size_df$count)
nfrags_nfr <- sum(size_df[size < 147]$count)
nfrags_mono <- sum(size_df[size >= 147 & size < 294]$count)
nfrags_multi <- sum(size_df[size >= 294]$count)
fraction_nfr <- ifelse(nfrags_total > 0, nfrags_nfr / nfrags_total * 100, 0)
fraction_mono_nucleosome <- ifelse(nfrags_total > 0, nfrags_mono / nfrags_total * 100, 0)
fraction_multi_nucleosome <- ifelse(nfrags_total > 0, nfrags_multi / nfrags_total * 100, 0)

js_data <- paste0(
  "var fragmentData = [",
  paste0(sprintf("{size:%d, percentage:%.2f}", size_df$size, size_df$percentage), collapse = ","),
  "];"
)

writeLines(js_data, file.path(outdir, paste0(sample_name, "_fragment_data.js")))

tryCatch({
    p7 <- ggplot()+
    geom_tile(data=pbmc@meta.data,aes(x=coor_x,y=coor_y,fill=nFrags))+
    scale_fill_gradientn(colours = rainbow(6))+
    theme_void()+
    coord_fixed()+
    theme(text=element_text(size=16))
    p8 <- ggplot()+
    geom_tile(data=pbmc@meta.data,aes(x=coor_x,y=coor_y,fill=TSS.enrichment))+
    scale_fill_gradientn(colours = rainbow(6))+
    theme_void()+
    coord_fixed()+
    theme(text=element_text(size=16))

    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc.svg")), p7 + p8, width = 12, height = 6, device = "svg")
    # Also emit a standalone TSS spatial map for the report's section-4 dropdown
    # (model0525 template: "bins_under_tissue_TSS.svg" option).
    ggsave(file.path(outdir, paste0(sample_name, "_bins_under_tissue_TSS.svg")), p8, width = 8, height = 8, device = "svg")
    message("[Signac] Spatial QC plot saved")
}, error = function(e) {
    message("[Signac] Warning: Spatial QC plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="Spatial QC failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc.svg")), blank, width = 12, height = 6, device = "svg")
    ggsave(file.path(outdir, paste0(sample_name, "_bins_under_tissue_TSS.svg")), blank, width = 8, height = 8, device = "svg")
})

tryCatch({
    pre_cmd <- paste0("pigz -dc -p ", WORKERS, " ", fragpath,
      " | awk 'BEGIN{OFS=\"\\t\"} $1 != \"chrM\" { print $4 }'")
    pre_dt <- fread(cmd = pre_cmd, header = FALSE, col.names = "CB")
    pre_dt <- pre_dt[CB != ""]
    pre_counts <- pre_dt[, .(nFrags = .N), by = CB]
    pre_counts[, coor_x := suppressWarnings(as.numeric(sapply(CB, function(x) strsplit(x, "_")[[1]][1])))]
    pre_counts[, coor_y := suppressWarnings(as.numeric(sapply(CB, function(x) strsplit(x, "_")[[1]][2])))]
    pre_counts <- pre_counts[!is.na(coor_x) & !is.na(coor_y)]

    p_pre <- ggplot(pre_counts) +
        geom_tile(aes(x = coor_x, y = coor_y, fill = nFrags)) +
        scale_fill_gradientn(colours = rainbow(6)) +
        theme_void() + coord_fixed() +
        theme(text = element_text(size = 16))
    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc_pre.svg")), p_pre,
           width = 6, height = 6, device = "svg")
    message("[Signac] Pre-cut spatial nFrags plot saved")
}, error = function(e) {
    message("[Signac] Warning: Pre-cut spatial plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Pre-cut spatial failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc_pre.svg")), blank,
           width = 6, height = 6, device = "svg")
})

saveRDS(pbmc, file.path(outdir, paste0(sample_name, "_filtered_seurat_object.rds")))

stats_file <- file.path(outdir, paste0(sample_name, "_stats.txt"))
stats_df <- data.frame(
    metric = c(
        "total_barcodes",
        "valid_barcodes",
        "median_nFrags",
        "median_TSS",
        "median_FRiP",
        "median_blacklist",
        "fraction_nfr",
        "fraction_mono_nucleosome",
        "fraction_multi_nucleosome",
        "n_cells"
    ),
    value = c(
        length(barcodes),
        ncol(pbmc),
        round(median(pbmc@meta.data$nFrags, na.rm = TRUE), 0),
        round(median(pbmc@meta.data$TSS.enrichment, na.rm = TRUE), 2),
        round(median(pbmc@meta.data$FRiP, na.rm = TRUE) * 100, 1),
        round(median(pbmc@meta.data$blacklist_ratio, na.rm = TRUE) * 100, 2),
        round(fraction_nfr, 2),
        round(fraction_mono_nucleosome, 2),
        round(fraction_multi_nucleosome, 2),
        ncol(pbmc)
    )
)
write.table(stats_df, stats_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message(sprintf("[Signac] Stats file: %s", stats_file))

message(sprintf("[Signac] Analysis complete! Results saved to: %s", outdir))
message(sprintf("[Signac] - Peaks: %s", file.path(outdir, paste0(sample_name, "_peaks.rds"))))
message(sprintf("[Signac] - QC violin: %s", file.path(outdir, paste0(sample_name, "_qc_violin.svg"))))
message(sprintf("[Signac] - TSS scatter: %s", file.path(outdir, paste0(sample_name, "_tss_scatter.svg"))))
message(sprintf("[Signac] - Fragment size data: %s (CSV) + %s (JS)",
        file.path(outdir, paste0(sample_name, "_fragment_size.csv")),
        file.path(outdir, paste0(sample_name, "_fragment_data.js"))))
message(sprintf("[Signac] - Spatial QC: %s", file.path(outdir, paste0(sample_name, "_spatial_qc.svg"))))
message(sprintf("[Signac] - Cluster plot (res=0.8): %s", file.path(outdir, paste0(sample_name, "_clusters_res0_8.svg"))))
message(sprintf("[Signac] - Seurat object: %s", file.path(outdir, paste0(sample_name, "_filtered_seurat_object.rds"))))
