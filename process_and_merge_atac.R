#!/usr/bin/env Rscript
# ==============================================================================
# process_and_merge_atac.R
# ==============================================================================
# Modular, multi-sample pipeline for single-nucleus 10x Multiome ATAC-seq QC,
# preprocessing, barcode harmonization, and recombination with pre-processed
# RNA Seurat objects.
#
# USAGE MODES:
# ------------------------------------------------------------------------------
# 1. Single-sample CLI mode (ideal for SLURM array jobs):
#    Rscript process_and_merge_atac.R \
#      --sample_id GG001AC \
#      --rna_rds /path/to/GG001AC_rna_processed.rds \
#      --peaks /path/to/ganglioglioma_peaks_GG001AC_peaks.narrowPeak \
#      --fragments /path/to/atac_fragments.tsv.gz \
#      --metadata /path/to/per_barcode_metrics.csv \
#      --output_dir /path/to/outputs/GG001AC
#
# 2. Batch mode via sample sheet:
#    Rscript process_and_merge_atac.R --sample_sheet samples.csv
#
# SAMPLE SHEET FORMAT (samples.csv):
# sample_id,rna_rds_path,peaks_path,fragments_path,metadata_path,output_dir
# GG001AC,/data/rna/GG001AC.rds,/data/peaks/GG001AC.narrowPeak,/data/atac/fragments.tsv.gz,/data/meta/per_barcode.csv,/data/out/GG001AC
# ==============================================================================
# Early CLI check for --help / usage without loading heavy libraries
# ==============================================================================
raw_cli_args <- commandArgs(trailingOnly = TRUE)
if (length(raw_cli_args) == 0 || any(c("-h", "--help") %in% raw_cli_args)) {
  cat("\n==============================================================================\n")
  cat("process_and_merge_atac.R: Modular Single-Cell Multiome Processing & Recombination\n")
  cat("==============================================================================\n")
  cat("Usage 1 (Batch Mode via Sample Sheet):\n")
  cat("  Rscript process_and_merge_atac.R --sample_sheet samples.csv [--output_dir <dir>]\n\n")
  cat("Usage 2 (Single-Sample Mode for SLURM Array Jobs):\n")
  cat("  Rscript process_and_merge_atac.R \\\n")
  cat("    --sample_id GG001AC \\\n")
  cat("    --rna_rds /path/to/rna_processed.rds \\\n")
  cat("    --peaks /path/to/peaks.narrowPeak \\\n")
  cat("    --fragments /path/to/atac_fragments.tsv.gz \\\n")
  cat("    [--metadata /path/to/per_barcode_metrics.csv] \\\n")
  cat("    --output_dir /path/to/output/GG001AC\n\n")
  cat("Optional Quality Control & Dimension Reduction Parameters:\n")
  cat("  --min_counts_atac <int>        Minimum ATAC counts per cell [default: 1000]\n")
  cat("  --max_counts_atac <int>        Maximum ATAC counts [default: auto (mean + 2*SD)]\n")
  cat("  --max_nucleosome_signal <num>  Maximum nucleosome signal [default: 2]\n")
  cat("  --min_tss_enrichment <num>     Minimum TSS enrichment score [default: 2]\n")
  cat("  --dims_rna <int>               Number of RNA PCA dimensions for WNN [default: 40]\n")
  cat("  --dims_atac_start <int>        First ATAC LSI dimension to use (skipping depth) [default: 2]\n")
  cat("  --dims_atac_end <int>          Last ATAC LSI dimension to use [default: 30]\n")
  cat("==============================================================================\n\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(future)
})

# Avoid parallel memory crashes on large genomic matrices
options(future.globals.maxSize = 8000 * 1024^2)

# ==============================================================================
# 1. Command-Line Argument Parser
# ==============================================================================
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default parameter values
  params <- list(
    sample_sheet           = NULL,
    sample_id              = NULL,
    rna_rds                = NULL,
    peaks                  = NULL,
    fragments              = NULL,
    metadata               = NULL,
    output_dir             = "./output",
    min_counts_atac        = 1000,
    max_counts_atac        = NULL, # NULL -> auto-calculate (mean + 2*SD)
    max_nucleosome_signal  = 2,
    min_tss_enrichment     = 2,
    dims_rna               = 40,
    dims_atac_start        = 2,
    dims_atac_end          = 30
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      key <- gsub("-", "_", key)
      if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
        val <- args[i + 1]
        i <- i + 2
      } else {
        val <- TRUE
        i <- i + 1
      }
      
      # Type conversions
      if (key %in% c("min_counts_atac", "max_counts_atac", "max_nucleosome_signal", 
                     "min_tss_enrichment", "dims_rna", "dims_atac_start", "dims_atac_end")) {
        val <- as.numeric(val)
      }
      params[[key]] <- val
    } else {
      i <- i + 1
    }
  }
  
  return(params)
}

# ==============================================================================
# 2. Annotation & Reference Setup
# ==============================================================================
setup_annotations <- function() {
  message("[SETUP] Initializing EnsDb.Hsapiens.v86 annotations for hg38...")
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  return(annotations)
}

# ==============================================================================
# 3. Peak File Ingestion
# ==============================================================================
load_and_clean_peaks <- function(peaks_path) {
  message(sprintf("  -> Reading peaks from: %s", peaks_path))
  
  if (grepl("\\.narrowPeak(\\.gz)?$", peaks_path, ignore.case = TRUE)) {
    # MACS2 narrowPeak format (BED6+4)
    peaks_df <- read.table(
      file = peaks_path,
      header = FALSE,
      stringsAsFactors = FALSE,
      col.names = c("chr", "start", "end", "name", "score", "strand", 
                    "fold_change", "neg_log10pvalue", "neg_log10qvalue", "relative_summit_pos")
    )
    
    # Optional q-value filtering (FDR <= 0.05 / -log10(q) >= 1.3)
    if ("neg_log10qvalue" %in% colnames(peaks_df)) {
      n_raw <- nrow(peaks_df)
      peaks_df <- peaks_df[peaks_df$neg_log10qvalue >= 1.3, ]
      message(sprintf("     Retained %d / %d peaks with -log10(qvalue) >= 1.3", nrow(peaks_df), n_raw))
    }
    
    peaks <- makeGRangesFromDataFrame(
      df = peaks_df,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = TRUE
    )
  } else {
    # Standard 3-column or 6-column BED format
    bed_df <- read.table(file = peaks_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(bed_df)[1:3] <- c("chr", "start", "end")
    peaks <- makeGRangesFromDataFrame(
      df = bed_df,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = TRUE
    )
  }
  
  seqlevelsStyle(peaks) <- "UCSC"
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  
  # Remove ENCODE unified blacklisted regions
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  message(sprintf("     Final standard chromosome, blacklist-filtered peaks: %d", length(peaks)))
  
  return(peaks)
}

# ==============================================================================
# 4. Core Single-Sample Processing Function
# ==============================================================================
process_sample <- function(sample_id, rna_rds_path, peaks_path, fragments_path, 
                           metadata_path = NULL, output_dir = "./output", 
                           annotations, params) {
  
  cat("\n==============================================================================\n")
  message(sprintf("[SAMPLE] Starting processing for Sample ID: %s", sample_id))
  cat("==============================================================================\n")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # Step 4.1: Load Pre-Processed RNA Object to Identify Valid Nuclei
  # ----------------------------------------------------------------------------
  message(sprintf("  [1/6] Loading pre-processed RNA Seurat object: %s", rna_rds_path))
  if (!file.exists(rna_rds_path)) {
    stop(sprintf("RNA RDS file not found: %s", rna_rds_path))
  }
  rna_obj <- readRDS(rna_rds_path)
  n_rna_initial <- ncol(rna_obj)
  message(sprintf("     Loaded RNA object containing %d pre-filtered cells.", n_rna_initial))
  
  # Extract raw cell barcodes expected in ATAC fragments (strip sample_id prefix if present)
  # e.g., 'GG001AC_AAACAGCCAAGGAATC-1' -> 'AAACAGCCAAGGAATC-1'
  raw_rna_barcodes <- sub(paste0("^", sample_id, "(_|-)"), "", colnames(rna_obj))
  
  # ----------------------------------------------------------------------------
  # Step 4.2: Ingest ATAC Peaks & Quantify Fragments on Valid Barcodes
  # ----------------------------------------------------------------------------
  message("  [2/6] Building ATAC ChromatinAssay...")
  peaks <- load_and_clean_peaks(peaks_path)
  
  if (!file.exists(fragments_path)) {
    stop(sprintf("ATAC fragments file not found: %s", fragments_path))
  }
  
  # Optimize: constrain Fragment object to RNA barcodes to skip empty droplet quantification
  message(sprintf("     Restricting fragment scanning to %d validated RNA barcodes...", length(raw_rna_barcodes)))
  frag_obj <- CreateFragmentObject(
    path = fragments_path,
    cells = raw_rna_barcodes
  )
  
  message("     Quantifying counts into peak-by-cell FeatureMatrix...")
  atac_counts <- FeatureMatrix(
    fragments = frag_obj,
    features = peaks,
    cells = raw_rna_barcodes
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = frag_obj,
    min.cells = 5,
    annotation = annotations
  )
  
  atac_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sample_id
  )
  
  n_atac_quantified <- ncol(atac_obj)
  
  # Attach per-barcode metadata if provided (e.g. from Cell Ranger ARC)
  if (!is.null(metadata_path) && file.exists(metadata_path)) {
    message(sprintf("     Attaching Cell Ranger ARC per-barcode metrics: %s", metadata_path))
    meta_df <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
    if ("barcode" %in% colnames(meta_df)) {
      rownames(meta_df) <- meta_df$barcode
      matched_meta <- meta_df[intersect(rownames(meta_df), colnames(atac_obj)), ]
      atac_obj <- AddMetaData(atac_obj, metadata = matched_meta)
    }
  }
  
  # ----------------------------------------------------------------------------
  # Step 4.3: Compute Signac ATAC QC Metrics (Explicitly avoiding ATACqc)
  # ----------------------------------------------------------------------------
  message("  [3/6] Computing Signac QC metrics (NucleosomeSignal & TSSEnrichment)...")
  DefaultAssay(atac_obj) <- "ATAC"
  
  # 1. Nucleosome banding pattern (mononucleosomal / nucleosome-free ratio)
  atac_obj <- NucleosomeSignal(object = atac_obj)
  
  # 2. Transcription Start Site (TSS) enrichment score
  atac_obj <- TSSEnrichment(object = atac_obj, fast = FALSE)
  
  # Define TSS enrichment grouping for diagnostic profile visualization
  min_tss <- params$min_tss_enrichment
  atac_obj$high.tss <- ifelse(atac_obj$TSS.enrichment >= min_tss, "High_TSS", "Low_TSS")
  
  # ----------------------------------------------------------------------------
  # Step 4.4: Export Diagnostic QC Plots to PDF
  # ----------------------------------------------------------------------------
  qc_pdf_path <- file.path(output_dir, paste0(sample_id, "_atac_qc_diagnostics.pdf"))
  message(sprintf("     Exporting diagnostic QC plots to: %s", qc_pdf_path))
  
  pdf(qc_pdf_path, width = 12, height = 8)
  
  p_vln <- VlnPlot(
    object = atac_obj,
    features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0.1,
    log = TRUE
  ) + plot_annotation(title = paste0(sample_id, " - Pre-Filter ATAC QC Distributions"))
  print(p_vln)
  
  # TSS Enrichment signal profile across TSS +/- 2000 bp
  tryCatch({
    p_tss <- TSSPlot(atac_obj, group.by = "high.tss") + 
      ggtitle(paste0(sample_id, " - TSS Enrichment Signal Profile (cutoff = ", min_tss, ")"))
    print(p_tss)
  }, error = function(e) {
    message(sprintf("     Notice: Could not generate TSSPlot: %s", e$message))
  })
  
  # Correlation scatter: nCount_ATAC vs TSS.enrichment
  p_scatter <- ggplot(atac_obj@meta.data, aes(x = nCount_ATAC, y = TSS.enrichment, color = nucleosome_signal < params$max_nucleosome_signal)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_x_log10() +
    geom_hline(yintercept = min_tss, linetype = "dashed", color = "firebrick") +
    scale_color_manual(values = c("TRUE" = "#1a73e8", "FALSE" = "#d93025"), name = "Pass Nucleosome") +
    theme_minimal() +
    labs(title = paste0(sample_id, " - TSS vs ATAC Depth"), x = "nCount_ATAC (log10)", y = "TSS Enrichment Score")
  print(p_scatter)
  
  dev.off()
  
  # ----------------------------------------------------------------------------
  # Step 4.5: Filter ATAC Cells using Configurable Thresholds
  # ----------------------------------------------------------------------------
  min_c <- params$min_counts_atac
  max_c <- params$max_counts_atac
  if (is.null(max_c) || is.na(max_c)) {
    # Dynamically compute upper depth cutoff: mean + 2*SD (capped at 100k)
    max_c <- round(mean(atac_obj$nCount_ATAC) + 2 * sd(atac_obj$nCount_ATAC))
    max_c <- min(max(max_c, 50000), 100000)
    message(sprintf("     Auto-calculated upper depth cutoff (mean + 2*SD): %d", max_c))
  }
  max_nuc <- params$max_nucleosome_signal
  
  message(sprintf("     Applying ATAC QC Filters: %d <= nCount_ATAC <= %d, nucleosome_signal <= %.1f, TSS.enrichment >= %.1f",
                  min_c, max_c, max_nuc, min_tss))
  
  atac_filtered <- subset(
    x = atac_obj,
    subset = nCount_ATAC >= min_c &
             nCount_ATAC <= max_c &
             nucleosome_signal <= max_nuc &
             TSS.enrichment >= min_tss
  )
  n_atac_post_qc <- ncol(atac_filtered)
  message(sprintf("     ATAC cells passing QC: %d / %d (%.1f%%)", 
                  n_atac_post_qc, n_atac_quantified, (n_atac_post_qc / n_atac_quantified) * 100))
  
  # ----------------------------------------------------------------------------
  # Step 4.6: Barcode Harmonization & Recombination with RNA Object
  # ----------------------------------------------------------------------------
  message("  [4/6] Harmonizing cell barcodes and intersecting with RNA...")
  
  # Check if RNA barcodes carry a sample_id prefix (e.g. GG001AC_AAACAGCCAAGGAATC-1)
  rna_has_prefix <- any(startsWith(colnames(rna_obj), paste0(sample_id, "_")))
  atac_has_prefix <- any(startsWith(colnames(atac_filtered), paste0(sample_id, "_")))
  
  if (rna_has_prefix && !atac_has_prefix) {
    message(sprintf("     Harmonizing: Prepended '%s_' to ATAC barcodes to match RNA object...", sample_id))
    atac_filtered <- RenameCells(atac_filtered, add.cell.id = sample_id)
  } else if (!rna_has_prefix && atac_has_prefix) {
    message(sprintf("     Harmonizing: Stripping '%s_' prefix from ATAC barcodes to match RNA object...", sample_id))
    new_names <- sub(paste0("^", sample_id, "_"), "", colnames(atac_filtered))
    atac_filtered <- RenameCells(atac_filtered, new.names = new_names)
  }
  
  # Intersect barcodes
  common_cells <- intersect(colnames(rna_obj), colnames(atac_filtered))
  n_final_multiome <- length(common_cells)
  
  message(sprintf("     Common Multiome cells: %d", n_final_multiome))
  message(sprintf("     Retention vs initial RNA:  %.1f%%", (n_final_multiome / n_rna_initial) * 100))
  message(sprintf("     Retention vs post-QC ATAC: %.1f%%", (n_final_multiome / n_atac_post_qc) * 100))
  
  if (n_final_multiome < 50) {
    warning(sprintf("CRITICAL: Extremely low number of overlapping cells (%d) for %s! Check barcode formats.", 
                    n_final_multiome, sample_id))
  }
  
  # Subset both modalities to common cells
  rna_sub <- subset(rna_obj, cells = common_cells)
  atac_sub <- subset(atac_filtered, cells = common_cells)
  
  # Create unified Multiome object
  multiome_obj <- rna_sub
  multiome_obj[["ATAC"]] <- atac_sub[["ATAC"]]
  
  # Merge ATAC-specific metadata into multiome_obj without overwriting existing RNA metadata
  atac_meta_cols <- setdiff(colnames(atac_sub@meta.data), colnames(multiome_obj@meta.data))
  for (col in atac_meta_cols) {
    multiome_obj@meta.data[[col]] <- atac_sub@meta.data[colnames(multiome_obj), col]
  }
  
  # ----------------------------------------------------------------------------
  # Step 4.7: Joint Dimensionality Reduction (RNA + ATAC + WNN)
  # ----------------------------------------------------------------------------
  message("  [5/6] Performing joint dimensionality reduction (LSI, PCA, WNN)...")
  
  # 1. ATAC processing (TF-IDF + SVD / LSI)
  DefaultAssay(multiome_obj) <- "ATAC"
  multiome_obj <- RunTFIDF(multiome_obj, verbose = FALSE)
  multiome_obj <- FindTopFeatures(multiome_obj, min.cutoff = "q0")
  multiome_obj <- RunSVD(multiome_obj, reduction.name = "lsi", reduction.key = "lsi_", verbose = FALSE)
  
  # Skip 1st LSI dimension (technical depth correlation)
  lsi_dims <- params$dims_atac_start:params$dims_atac_end
  multiome_obj <- RunUMAP(
    multiome_obj,
    reduction = "lsi",
    dims = lsi_dims,
    reduction.name = "umap.atac",
    reduction.key = "atacUMAP_",
    verbose = FALSE
  )
  
  # 2. Verify / Refresh RNA PCA on the common cells
  DefaultAssay(multiome_obj) <- "RNA"
  if (!"pca" %in% Reductions(multiome_obj)) {
    message("     RNA PCA not detected; running PCA on normalized RNA assay...")
    multiome_obj <- RunPCA(multiome_obj, verbose = FALSE)
  }
  
  # 3. Weighted Nearest Neighbor (WNN) Multimodal Integration
  rna_dims <- 1:params$dims_rna
  message(sprintf("     Computing WNN neighbors (RNA dims 1:%d, ATAC dims %d:%d)...", 
                  params$dims_rna, params$dims_atac_start, params$dims_atac_end))
  
  multiome_obj <- FindMultiModalNeighbors(
    object = multiome_obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(rna_dims, lsi_dims),
    verbose = FALSE
  )
  
  multiome_obj <- RunUMAP(
    object = multiome_obj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    verbose = FALSE
  )
  
  # Multi-panel diagnostic UMAP plot
  umap_pdf_path <- file.path(output_dir, paste0(sample_id, "_joint_wnn_umaps.pdf"))
  pdf(umap_pdf_path, width = 15, height = 5)
  p1 <- DimPlot(multiome_obj, reduction = if ("umap" %in% Reductions(multiome_obj)) "umap" else "umap.atac", 
                label = FALSE) + ggtitle(paste0(sample_id, " - RNA UMAP"))
  p2 <- DimPlot(multiome_obj, reduction = "umap.atac", label = FALSE) + ggtitle(paste0(sample_id, " - ATAC LSI UMAP"))
  p3 <- DimPlot(multiome_obj, reduction = "wnn.umap", label = FALSE) + ggtitle(paste0(sample_id, " - Joint WNN UMAP"))
  print(p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # ----------------------------------------------------------------------------
  # Step 4.8: Save Final Unified Multiome Object & Summary Stats
  # ----------------------------------------------------------------------------
  final_rds_path <- file.path(output_dir, paste0(sample_id, "_multiome_merged.rds"))
  message(sprintf("  [6/6] Saving unified Multiome Seurat object to: %s", final_rds_path))
  saveRDS(multiome_obj, file = final_rds_path)
  
  summary_stats <- data.frame(
    sample_id            = sample_id,
    initial_rna_cells    = n_rna_initial,
    quantified_atac_cells= n_atac_quantified,
    post_qc_atac_cells   = n_atac_post_qc,
    final_multiome_cells = n_final_multiome,
    pct_rna_retained     = round((n_final_multiome / n_rna_initial) * 100, 2),
    pct_atac_retained    = round((n_final_multiome / n_atac_post_qc) * 100, 2),
    median_nCount_RNA    = median(multiome_obj$nCount_RNA),
    median_nCount_ATAC   = median(multiome_obj$nCount_ATAC),
    median_TSS_enrichment= median(multiome_obj$TSS.enrichment),
    median_nucleosome_sig= median(multiome_obj$nucleosome_signal),
    output_rds           = final_rds_path,
    stringsAsFactors     = FALSE
  )
  
  csv_path <- file.path(output_dir, paste0(sample_id, "_cell_retention_summary.csv"))
  write.csv(summary_stats, file = csv_path, row.names = FALSE)
  message(sprintf("     Saved cell retention summary to: %s", csv_path))
  
  message(sprintf("[DONE] Sample %s processing completed successfully!\n", sample_id))
  return(summary_stats)
}

# ==============================================================================
# 5. Main Execution Driver
# ==============================================================================
main <- function() {
  params <- parse_args()
  
  # Setup genome annotations once
  annotations <- setup_annotations()
  
  all_summaries <- list()
  
  if (!is.null(params$sample_sheet) && file.exists(params$sample_sheet)) {
    # --------------------------------------------------------------------------
    # Batch Mode: Process All Samples in CSV Sample Sheet
    # --------------------------------------------------------------------------
    message(sprintf("[BATCH MODE] Reading sample sheet from: %s", params$sample_sheet))
    samples_df <- read.csv(params$sample_sheet, header = TRUE, stringsAsFactors = FALSE)
    
    required_cols <- c("sample_id", "rna_rds_path", "peaks_path", "fragments_path")
    missing_cols <- setdiff(required_cols, colnames(samples_df))
    if (length(missing_cols) > 0) {
      stop(sprintf("Sample sheet is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    for (i in seq_len(nrow(samples_df))) {
      row <- samples_df[i, ]
      out_d <- if ("output_dir" %in% colnames(row) && !is.na(row$output_dir) && row$output_dir != "") {
        row$output_dir
      } else {
        file.path(params$output_dir, row$sample_id)
      }
      
      meta_p <- if ("metadata_path" %in% colnames(row) && !is.na(row$metadata_path)) row$metadata_path else NULL
      
      stat <- tryCatch({
        process_sample(
          sample_id      = row$sample_id,
          rna_rds_path   = row$rna_rds_path,
          peaks_path     = row$peaks_path,
          fragments_path = row$fragments_path,
          metadata_path  = meta_p,
          output_dir     = out_d,
          annotations    = annotations,
          params         = params
        )
      }, error = function(e) {
        message(sprintf("ERROR processing sample %s: %s", row$sample_id, e$message))
        return(NULL)
      })
      
      if (!is.null(stat)) {
        all_summaries[[row$sample_id]] <- stat
      }
    }
    
    # Save combined multi-sample summary CSV
    if (length(all_summaries) > 0) {
      combined_summary <- do.call(rbind, all_summaries)
      batch_csv_path <- file.path(params$output_dir, "batch_multiome_processing_summary.csv")
      write.csv(combined_summary, file = batch_csv_path, row.names = FALSE)
      message(sprintf("\n[SUMMARY] Batch processing completed. Summary saved to: %s", batch_csv_path))
    }
    
  } else if (!is.null(params$sample_id) && !is.null(params$rna_rds) && 
             !is.null(params$peaks) && !is.null(params$fragments)) {
    # --------------------------------------------------------------------------
    # Single-Sample CLI Mode (SLURM Array Friendly)
    # --------------------------------------------------------------------------
    stat <- process_sample(
      sample_id      = params$sample_id,
      rna_rds_path   = params$rna_rds,
      peaks_path     = params$peaks,
      fragments_path = params$fragments,
      metadata_path  = params$metadata,
      output_dir     = params$output_dir,
      annotations    = annotations,
      params         = params
    )
  } else {
    cat("\nERROR: Insufficient arguments provided.\n\n")
    cat("Usage 1 (Batch Mode):\n")
    cat("  Rscript process_and_merge_atac.R --sample_sheet samples.csv\n\n")
    cat("Usage 2 (Single-Sample Mode):\n")
    cat("  Rscript process_and_merge_atac.R \\\n")
    cat("    --sample_id <ID> \\\n")
    cat("    --rna_rds <path_to_rna.rds> \\\n")
    cat("    --peaks <path_to_peaks.narrowPeak> \\\n")
    cat("    --fragments <path_to_fragments.tsv.gz> \\\n")
    cat("    --output_dir <output_directory>\n\n")
    quit(status = 1)
  }
}

if (!interactive()) {
  main()
}
