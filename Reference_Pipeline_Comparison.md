# Reference Pipeline Comparison: Human Cortical Development Multiome vs. Ganglioglioma Pipeline

**Date:** November 12, 2026  
**Project:** Ganglioglioma 10x Multiome Analysis  
**Reference Study:** Kriegstein Lab Human Cortical Development 10x Multiome Atlas (Wang et al., *Nature*)  
**Reference Codebase:** [`Human_Cortex_Dev_Multiome-main/RScript/`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript)  
**Current Pipeline:** [`process_and_merge_atac.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Multiome_workspace/process_and_merge_atac.R) & [`test_process_and_merge_atac.Rmd`](file:///Users/franz/Documents/Code/Ganglioglioma/Multiome_workspace/test_process_and_merge_atac.Rmd)  

---

## Executive Summary

This report provides a comprehensive architectural breakdown of the 22 non-spatial R scripts from the Kriegstein lab human cortical development 10x Multiome pipeline ([`Human_Cortex_Dev_Multiome-main/RScript/`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript)), compares their methodology against our current single-sample preprocessing pipeline ([`process_and_merge_atac.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Multiome_workspace/process_and_merge_atac.R)), and outlines high-impact modular extensions to integrate into our Ganglioglioma cohort workflow across anatomical zones (`GG001TB` Tumor Bulk, `GG001AMY` Amygdala, `GG001IF` Interface, `GG001AC` Access Cortex).

---

## Part 1: Script-by-Script Walkthrough of the Reference Pipeline

> [!NOTE]
> Scripts `04_MERFISH_data_loading.R` and `05_MERFISH_data_processing.R` pertain exclusively to MERFISH spatial transcriptomics and have been excluded per task specifications. Note that indices `13` and `22` do not exist in the source repository.

```
Human_Cortex_Dev_Multiome-main/RScript/
├── 00_snMultiome_data_integration.R                     # Multi-sample RPCA (RNA) + Reciprocal LSI (ATAC) + WNN
├── 01_snMultiome_cell_fitering.R                        # Ultra-high res over-clustering (wsnn) for doublets & contaminants
├── 02_snMultiome_data_processing.R                      # Final clustering, 3-tier hierarchy annotation & publication QC
├── 03_snMultiome_proportion_analysis.R                  # speckle (logit) + limma mixed model for cell composition
├── 06_AUCell.R                                          # SCENIC+ eRegulon activity scoring (AUCell)
├── 07_slingshot.R                                       # 8D WNN UMAP Gaussian Mixture clustering + Slingshot trajectory
├── 08_tradeSeq.R                                        # GAM modeling along pseudotime & kinetic module discovery
├── 09_L4_IT_subtypes.R                                  # Pseudobulk edgeR + glmmSeq mixed-effects DE modeling
├── 10_macroglial_lineages.R                             # Lineage-specific macroglial sub-clustering & surface markers
├── 11_Sorted_GPC_scRNA-seq_processing.R                 # In vitro glial progenitor scRNA-seq processing
├── 12_Shi2021_processing.R                              # Public ganglionic eminence interneuron reference processing
├── 14_cellchat.R                                        # CellChat ligand-receptor communication networks
├── 15_GPC-derived_IN_type_mapping_to_Shi2021.R          # Interneuron label transfer & SingleCellNet RF classification
├── 16_Tri-IPC-derived_Astrocyte_type_mapping_...        # Orthogene ortholog conversion & DiBella mapping
├── 17_Tri-IPC-derived_Astrocyte_type_mapping_...        # Projection of astrocyte states onto Multiome atlas
├── 18_Motif_analysis.R                                  # JASPAR2020 motif curation & chromVAR deviation scoring
├── 19_GWAS_finemapping_processing.R                     # coloc Bayesian fine-mapping around lead GWAS variants
├── 20_SCAVENGE_TRS_inference.R                          # gchromVAR + SCAVENGE network propagation for cell TRS
├── 21_Correlation_test_of_TRS_and_AUCell_...            # Trajectory correlation: Disease TRS vs. TF eRegulons
├── 23_SST_DEG_analysis.R                                # Pseudobulk limma-voom DE & ranked GSEA
├── 24_GPC_to_snMultiome_mapping_and_Sankey_plot.R       # SingleCellNet RF prediction & D3 Sankey river diagrams
└── 25_GBM_to_snMultiome_mapping.R                       # Glioblastoma-to-normal developmental atlas mapping
```

---

### Phase A: Core Multiome Integration, Filtering & Annotation (`00`–`02`)

#### 1. [`00_snMultiome_data_integration.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R) — *Multi-Sample RNA + ATAC Batch Integration & WNN*
Rather than keeping samples separate, this script pools Cell Ranger ARC outputs across all human fetal donors into a unified, batch-corrected multimodal object:
- **GEX Processing & RPCA Integration ([L15–L120](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R#L15-L120)):**
  - Loads `filtered_feature_bc_matrix.h5` and `per_barcode_metrics.csv` for each sample.
  - Scores cell-cycle phases (`CellCycleScoring`) and defines **Cell Cycle Difference (`CCD = S.Score - G2M.Score`)**. This strategy isolates G1 vs. S/G2M phase variance while preserving biological proliferative progenitor states.
  - Filters low-quality cells (`nCount_RNA > 1000 & < 25000`, `nFeature_RNA > 500`, `percent.mt < 10`) and strips mitochondrial (`^MT-`) and sex-linked (`chrX`/`chrY`) genes prior to normalization.
  - Runs `SCTransform(method = "glmGamPoi", vars.to.regress = "CCD")` per sample, selects 3,000 integration features via `SelectIntegrationFeatures()`, runs `PrepSCTIntegration()`, and performs **Reciprocal PCA (`rpca`) integration** (`FindIntegrationAnchors(reduction = "rpca", k.anchor = 20, reference = 2)` + `IntegrateData()`).
- **Consensus ATAC Peak Set & Reciprocal LSI (`rlsi`) Integration ([L122–L280](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R#L122-L280)):**
  - Reads each sample's `atac_peaks.bed` into `GRanges` and merges all sample peaks into a single master consensus peak set using `GenomicRanges::reduce()`, filtering peak boundaries to `20 < width < 10000`.
  - Quantifies the consensus peak matrix across samples using `Signac::FeatureMatrix()`, wraps each in a `ChromatinAssay` annotated with `EnsDb.Hsapiens.v86` (`hg38`), and computes `NucleosomeSignal()`, `TSSEnrichment(fast = FALSE)`, and Fraction of Reads in Peaks (`pct_reads_in_peaks`).
  - Filters ATAC cells (`atac_peak_region_fragments > 500 & < 50000`, `pct_reads_in_peaks > 15`, `nucleosome_signal < 2`, `TSS.enrichment > 2`), subsets to barcodes passing both GEX and ATAC filters, and merges the objects.
  - Runs `FindTopFeatures(min.cutoff = "q80")`, `RunTFIDF()`, `RunSVD()`, and computes **Reciprocal LSI (`rlsi`) batch integration** via `FindIntegrationAnchors(reduction = "rlsi", dims = 2:30)` followed by `IntegrateEmbeddings(new.reduction.name = "integrated_lsi")`.
- **Joint Weighted Nearest Neighbor (WNN) Graph ([L281–L320](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R#L281-L320)):**
  - Fuses integrated `SCT` PCA (`dims = 1:50`) and integrated ATAC `integrated_lsi` (`dims = 2:40`) via `FindMultiModalNeighbors()` and computes the joint embedding `RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap")`.

#### 2. [`01_snMultiome_cell_fitering.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/01_snMultiome_cell_fitering.R) — *High-Resolution Over-Clustering for Doublet & Contaminant Removal*
- Evaluates ~150 canonical brain marker genes across the multimodal manifold (`wnn.umap`) using `scCustomize::FeaturePlot_scCustom()`.
- Uses the SLM algorithm (`algorithm = 3`) to **intentionally over-cluster** the joint graph (`wsnn`) across fine resolutions (`3, 4, 5, 6, 7`), generating up to **219 micro-clusters** at `wsnn_res.7`.
- Systematically inspects micro-clusters with `Cluster_Highlight_Plot()` to flag and discard:
  1. **Anatomical contaminants:** Off-target subpallial/diencephalic clusters expressing *ISL1*, *SIX3*, *OTX2*, *GBX2* or lacking *FOXG1*/*EMX1*.
  2. **Heterotypic doublet micro-clusters:** Droplets co-expressing mutually exclusive lineage markers (e.g., *NRGN* [excitatory neuron] + *MBP* [oligodendrocyte]) that slipped past global single-cell metric thresholds.
- Re-calculates `FindMultiModalNeighbors()` and `wnn.umap` on the cleaned manifold (`multiome_031123_filtered.rds`).

#### 3. [`02_snMultiome_data_processing.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/02_snMultiome_data_processing.R) — *Final Clustering, Marker Discovery & Publication QC Figures*
- Clusters the filtered graph at `wsnn_res.5` (176 clusters), runs `FindAllMarkers(assay = "SCT", recorrect_umi = FALSE)`, and maps clusters via a curated lookup table (`annotation.csv`) into a 3-tier hierarchy: `class` -> `subclass` -> `type` (34 discrete developmental cell types).
- Computes cell-type marker signatures across all 34 `type` labels.
- Assembles comprehensive QC panels using `cowplot::plot_grid()`: aligned violin/boxplots (`nCount_RNA`, `nFeature_RNA`, `log10(atac_peak_region_fragments)`, `TSS.enrichment`, `nucleosome_signal`) alongside sample-level stacked composition barplots, accompanied by modality-isolated UMAPs (`gex.umap` on PCA 1:50; `atac.umap` on `integrated_lsi` 2:40).

---

### Phase B: Compositional & Lineage-Specific Sub-Analyses (`03`, `09`, `10`, `23`)

#### 4. [`03_snMultiome_proportion_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/03_snMultiome_proportion_analysis.R) — *Statistical Cell-Type Proportion Testing (`speckle` + `limma`)*
- Applies `speckle::getTransformedProps(dat$type, dat$dataset, transform = "logit")` to compute logit-transformed cell-type proportions per sample, overcoming compositional closure constraints.
- Employs `limma::duplicateCorrelation()` (blocking on donor ID to handle paired cortical regions from identical brains), `lmFit()`, and empirical Bayes moderation (`eBayes(robust = TRUE)`) to test for continuous age (`log2_age`) and regional enrichment across cell types.

#### 5. [`09_L4_IT_subtypes.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/09_L4_IT_subtypes.R) — *Pseudobulk Mixed-Effects Differential Expression (`edgeR` + `glmmSeq`)*
- Isolates Layer 4 intratelencephalic (`EN-L4-IT`) neurons, re-clusters them into fine subtypes, and aggregates raw single-cell counts into `(subtype × sample)` pseudobulks via `scuttle::aggregateAcrossCells()` (enforcing $\ge 50$ cells per pseudobulk).
- Filters lowly expressed transcripts (`edgeR::filterByExpr()`), calculates TMM normalization factors (`calcNormFactors`), estimates dispersions (`estimateDisp`), and fits a **Negative Binomial Generalized Linear Mixed Model** using `glmmSeq::glmmSeq(~ subtype + log2_age + (1 | dataset))` against a reduced model (`~ log2_age + (1 | dataset)`). This controls for donor-level random effects and prevents pseudoreplication artifacts.

#### 6. [`10_macroglial_lineages.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/10_macroglial_lineages.R) — *Glial Sub-Atlas & Surface Marker Profiling*
- Subsets the macroglial compartment (`RG-vRG`, `RG-tRG`, `RG-oRG`, `IPC-Glia` [Tri-IPC], `Astrocyte-Immature`, `Astrocyte-Protoplasmic`, `Astrocyte-Fibrous`, `OPC`, `Oligodendrocyte-Immature`, `Oligodendrocyte`).
- Recomputes a dedicated glial multimodal graph (`wsnn.Glia`) and UMAP, and evaluates cell-surface sorting markers (*ITGA2*, *CD38*, *F3*, *EGFR*, *PDGFRA*) using `scCustomize::Stacked_VlnPlot()`.

#### 7. [`23_SST_DEG_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/23_SST_DEG_analysis.R) — *Pseudobulk Condition DE (`limma-voom`) & Ranked GSEA*
- Aggregates single-cell counts by `(cell_type, individual, treatment)` using `scran::aggregateAcrossCells()` ($\ge 10$ cells).
- Executes paired pseudobulk testing via `scran::pseudoBulkDGE(..., design = ~ individual + treatment, method = "voom")` and performs ranked Gene Set Enrichment Analysis (`clusterProfiler::gseGO()`) on the resulting $t$-statistic distributions.

---

### Phase C: Multimodal Trajectory, Pseudotime & Gene Regulatory Networks (`06`–`08`)

#### 8. [`06_AUCell.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/06_AUCell.R) — *SCENIC+ Enhancer-Regulon (eRegulon) Activity Scoring*
- Loads enhancer-driven gene regulatory network gene sets (`582.sceplus.gmt`, inferred in Python SCENIC+ through joint RNA expression and ATAC co-accessibility linking).
- Uses `AUCell::AUCell_run()` on raw RNA count matrices to quantify per-cell TF eRegulon activation scores.
- Instantiates a dedicated `scenicplus_AUC` Seurat assay and maps transcription factor regulatory activities directly onto `wnn.umap`.

#### 9. [`07_slingshot.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/07_slingshot.R) — *8-Dimensional WNN Trajectory & Pseudotime Inference (`mclust` + `slingshot`)*
- Subsets the excitatory neurogenic lineage (`Radial glia`, `IPC-EN`, `Glutamatergic neuron`) and recalculates local multimodal connectivity (`wknn.EN`, `wsnn.EN`).
- **Multidimensional Projection Optimization ([L54–L90](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/07_slingshot.R#L54-L90)):** Computes an **8-dimensional WNN UMAP** (`n.components = 8`) to preserve high-dimensional manifold topology for trajectory calculations, in parallel with a 2D UMAP (`wnn.EN.umap2`) strictly for visualization.
- Clusters the 8D embeddings via Gaussian Mixture Models (`mclust::Mclust(rd, G = 23:28)`), runs `slingshot::getLineages()` (incorporating donor age `times = log2_age` to anchor root clusters in early Radial Glia toward terminal cortical plate layers), fits smooth splines with `getCurves(approx_points = 200)`, and projects the 8D curves back onto the 2D visual UMAP with `embedCurves()`.

#### 10. [`08_tradeSeq.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/08_tradeSeq.R) — *Dynamic Gene & eRegulon Modeling Along Pseudotime (`tradeSeq`)*
- Connects Slingshot pseudotime and lineage weights to `tradeSeq::fitGAM()` with 6 quantile-spaced knots. Employs `family = "betar"` (Beta regression) for bounded `[0, 1]` AUCell scores and negative binomial for raw transcript counts.
- Reconstructs smoothed temporal curves (`predictSmooth()`), clusters kinetic profiles into $k = 6$ temporal wave modules, and annotates functional transitions with `clusterProfiler::compareCluster(..., fun = "enrichGO")`.
- Executes four statistical tests across lineages:
  - `associationTest()`: Tests for significant expression dynamics over pseudotime.
  - `startVsEndTest()`: Identifies genes undergoing net activation/silencing from root to terminus.
  - `diffEndTest()`: Identifies lineage-divergent genes distinguishing terminal fates.
  - `earlyDETest(knots = c(2, 3))`: Pins down transient driver regulators acting at lineage bifurcation points.

---

### Phase D: Cell-Cell Communication (`14`)

#### 11. [`14_cellchat.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/14_cellchat.R) — *Ligand-Receptor Signaling Networks (`CellChat`)*
- Converts normalized `SCT` expression profiles and consensus `type` labels into a `CellChat` object querying `CellChatDB.human`.
- Computes communication probabilities (`computeCommunProb()`), prunes low-abundance interactions (`filterCommunication(min.cells = 10)`), models signaling pathways (`computeCommunProbPathway()`), calculates network centrality metrics (`netAnalysis_computeCentrality()`), and generates interaction circle plots (`netVisual_circle()`).

---

### Phase E: Reference Mapping & Tumor-to-Normal Developmental Projection (`11`–`12`, `15`–`17`, `24`–`25`)

#### 12. [`11_Sorted_GPC_scRNA-seq_processing.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/11_Sorted_GPC_scRNA-seq_processing.R) & [`12_Shi2021_processing.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/12_Shi2021_processing.R)
- Standard scRNA-seq processing of FACS-sorted in vitro Glial Progenitor Cell (GPC) cultures (`11`) and public human ganglionic eminence interneuron reference data (`GSE135827`, Shi et al. 2021) integrated via RPCA (`12`).

#### 13. [`15_GPC-derived_IN_type_mapping_to_Shi2021.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/15_GPC-derived_IN_type_mapping_to_Shi2021.R), [`16_Tri-IPC-derived_Astrocyte_type_mapping_to_DiBella2021.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/16_Tri-IPC-derived_Astrocyte_type_mapping_to_DiBella2021.R), [`17_Tri-IPC-derived_Astrocyte_type_mapping_to_snMultiome_data.R.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/17_Tri-IPC-derived_Astrocyte_type_mapping_to_snMultiome_data.R.R) & [`24_GPC_to_snMultiome_mapping_and_Sankey_plot.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/24_GPC_to_snMultiome_mapping_and_Sankey_plot.R)
- Establishes a **dual-classifier validation framework** for querying external single cells against atlas references:
  1. **Seurat Anchor Transfer:** `FindTransferAnchors(reference.reduction = "pca")` + `MapQuery(reduction.model = "umap")` + `TransferData()` (enforcing prediction confidence threshold `score > 0.5`).
  2. **`singleCellNet` (Random Forest Pair-Based Classifier):** Trains a multi-class Random Forest model on top rank-invariant gene pairs (`scn_train(nTopGenes = 60..200, nTopGenePairs = 150..200, nTrees = 1000)`), validates out-of-bag classification with `assess_comm()`, and scores query cells (`scn_predict()`). Results are rendered as probability heatmaps (`sc_hmClass()`) and interactive river diagrams (`networkD3::sankeyNetwork` in `24`). Cross-species mouse-to-human orthology conversion in `16` is managed via `orthogene::convert_orthologs()`.

#### 14. [`25_GBM_to_snMultiome_mapping.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/25_GBM_to_snMultiome_mapping.R) — *Brain Tumor (`GBmap`) to Normal Developmental Multiome Mapping*
- **Direct relevance to Ganglioglioma:** Trains a `singleCellNet` Random Forest classifier on the normal fetal multiome atlas and projects malignant primary brain tumor cells (`GBmap_extended_malignant_cell.rds`) onto normal developmental stages.
- Deconvolves neoplastic cell states (`AC-like`, `MES-like`, `NPC-like`, `OPC-like`) and patient-specific tumor biopsies against normal progenitor archetypes (`IPC-Glia` [Tri-IPC], `OPC`, `Astrocyte-Immature`, `Radial Glia`, `IN-dLGE-Immature`).

---

### Phase F: Chromatin Motifs, GWAS Fine-Mapping & Trait Relevance Scoring (`18`–`21`)

#### 15. [`18_Motif_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/18_Motif_analysis.R) — *JASPAR Motif Attachment & `chromVAR` TF Activity*
- Subsets ATAC peaks to canonical chromosomes (`BSgenome.Hsapiens.UCSC.hg38`), extracts vertebrate core Position Frequency Matrices from `JASPAR2020` via `TFBSTools::getMatrixSet()`, registers motif matches in the `ChromatinAssay` using `Signac::AddMotifs()`, and computes single-cell motif accessibility deviation z-scores through `Signac::RunChromVAR()`.

#### 16. [`19_GWAS_finemapping_processing.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/19_GWAS_finemapping_processing.R), [`20_SCAVENGE_TRS_inference.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/20_SCAVENGE_TRS_inference.R) & [`21_Correlation_test_of_TRS_and_AUCell_along_the_pseudotime.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/21_Correlation_test_of_TRS_and_AUCell_along_the_pseudotime.R)
- `19`: Executes Bayesian genetic fine-mapping across GWAS risk loci ($\pm 500\text{ kb}$ around lead SNPs) using Approximate Bayes Factors (`coloc::finemap.abf`), extracts 95% credible causal SNP sets, and lifts positions to `hg38` via `rtracklayer::liftOver()`.
- `20`: Quantifies single-cell causal variant enrichment in open chromatin using `gchromVAR::importBedScore()` and `computeWeightedDeviations()`. Performs **`SCAVENGE`** random-walk network propagation (`getmutualknn()` on `integrated_lsi` + `randomWalk_sparse()`) to compute single-cell **Trait Relevance Scores (TRS)**.
- `21`: Calculates Pearson correlation and FDR statistics between single-cell disease TRS scores and SCENIC+ eRegulon `AUCell` trajectories along Slingshot pseudotime, distinguishing developmental driver TFs from disease-associated regulators.

---

## Part 2: Systematic Pipeline Comparison

| Analytical Dimension | Our Ganglioglioma Pipeline ([`process_and_merge_atac.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Multiome_workspace/process_and_merge_atac.R)) | Reference Cortical Dev Pipeline ([`RScript/`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript)) | Architectural Trade-offs & Recommendations |
| :--- | :--- | :--- | :--- |
| **Execution Architecture** | Modular CLI/batch script driven by `optparse` and SLURM array automation; consumes external pre-processed RNA `.rds` objects. | Monolithic serial R scripts (`00`–`25`) operating interactively from raw Cell Ranger ARC `.h5` files. | **Our architecture is far better suited for HPC array execution** and decoupled RNA/ATAC preprocessing; the reference pipeline provides a superior blueprint for multi-sample cohort integration and downstream biological modeling. |
| **Peak Calling Strategy** | **Custom BAM-derived MACS2 `.narrowPeak`** (`--nomodel --shift -100 --extsize 200`), filtered by FDR ($- \log_{10}(q) \ge 1.3$) and unified ENCODE blacklists. | Standard Cell Ranger ARC `atac_peaks.bed` per sample without explicit blacklisting in `00`. | **Our peak calling is biologically superior** (MACS2 sub-nucleosomal shift + ENCODE blacklisting generates cleaner peak calls than 10x defaults). |
| **Multi-Sample Feature Space** | Currently quantifies each sample (`GG001AC`, `GG001IF`, `GG001AMY`, `GG001TB`) against its own individual MACS2 peak set. | Unifies all sample peaks into a single master consensus set via `GenomicRanges::reduce()` (`20 < width < 10000`) before running `FeatureMatrix()`. | **The reference pipeline's consensus approach is mandatory for multi-sample merging.** Distinct MACS2 peak boundaries prevent direct feature concatenation across samples without a shared reference peak set. |
| **ATAC QC Filtration** | `nCount_ATAC`, `nucleosome_signal < 2`, `TSS.enrichment > 2` (via `NucleosomeSignal` + `TSSEnrichment(fast = FALSE)`). | Same (`nucleosome_signal < 2`, `TSS.enrichment > 2`), plus **`pct_reads_in_peaks > 15` (FRiP)** from `per_barcode_metrics.csv`. | Adding the Fraction of Reads in Peaks (FRiP) metric from Cell Ranger ARC metrics provides an extra filter against ambient background DNA. |
| **Doublet & Artifact Cleansing** | Relies on single-cell thresholds (`mean + 2*SD` or configurable limits). | **High-resolution WNN over-clustering (`wsnn_res.7`, 219 clusters)** in `01` to isolate and purge heterotypic doublets (*NRGN*+/*MBP*+). | **The reference approach is highly recommended for Ganglioglioma**, where glioneuronal doublets can confound true neoplastic dysplastic cells. |
| **Cross-Sample Batch Integration** | Computes independent per-sample WNN manifolds (`pca` + `lsi`). | **RPCA** on `SCT` (with `CCD` cell-cycle regression) + **Reciprocal LSI (`rlsi`)** on `ATAC` $\to$ joint multi-sample WNN. | Reference pipeline `00` provides the exact mathematical framework needed for integrating `GG001AC`, `IF`, `AMY`, and `TB`. |
| **Lineage & Trajectory Modeling** | Not currently implemented. | **8D WNN Slingshot trajectory (`07`)** + **`tradeSeq` GAM kinetic modeling (`08`)** + **SCENIC+ eRegulons (`06`)**. | Provides an established framework for dissecting neoplastic glioneuronal differentiation hierarchies. |
| **Tumor-to-Normal Reference Mapping** | Not currently implemented. | **`singleCellNet` Random Forest pair classifier** + **Seurat `MapQuery` (`25`)**. | Enables direct projection of Ganglioglioma neoplastic cells onto normal fetal cortical progenitors and adult cortical reference states. |

---

## Part 3: High-Priority Modules for Integration into the Ganglioglioma Pipeline

Based on our four-region multiome cohort (`GG001TB` Tumor Bulk, `GG001AMY` Amygdala, `GG001IF` Interface, `GG001AC` Access Cortex), the following six modules provide the highest biological and methodological return on investment:

### 1. Unified Consensus MACS2 Peak Set & Reciprocal LSI (`rlsi`) Integration
- **Source:** [`00_snMultiome_data_integration.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R)
- **Rationale:** Currently, `GG001AC`, `GG001IF`, `GG001AMY`, and `GG001TB` possess slightly shifted MACS2 peak boundaries. To merge them into a single cohort Seurat object, they must share an identical row coordinate space.
- **Implementation Strategy:**
  1. Merge individual sample MACS2 `.narrowPeak` files with `GenomicRanges::reduce()`, enforcing `width > 20 & width < 10000` and discarding `blacklist_hg38_unified`.
  2. Quantify each sample's fragments against this shared consensus set using `FeatureMatrix()`.
  3. Integrate the four ATAC spaces using Signac **Reciprocal LSI (`rlsi`)** (`FindIntegrationAnchors(reduction = "rlsi", dims = 2:30)` + `IntegrateEmbeddings()`) and combine with RNA RPCA to generate a unified multi-sample WNN UMAP.

### 2. High-Resolution WNN Over-Clustering for Glioneuronal Doublet Dissection
- **Source:** [`01_snMultiome_cell_fitering.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/01_snMultiome_cell_fitering.R)
- **Rationale:** Ganglioglioma is a biphasic glioneuronal neoplasm composed of dysplastic ganglion cells and neoplastic glial elements. Mechanical tissue dissociation frequently produces artificial neuron-oligodendrocyte doublets (e.g., *NRGN*+ / *MBP*+ or *SYT1*+ / *MOG*+) that pass standard UMI cutoffs.
- **Implementation Strategy:**
  1. Over-cluster the joint WNN graph (`wsnn`) at `resolution = 4` to `7` using the SLM algorithm (`algorithm = 3`).
  2. Flag and remove artificial micro-clusters that co-express mutually exclusive mature lineage markers, ensuring true neoplastic cells are not artifactual doublets.

### 3. Tumor-to-Normal Reference Mapping via `SingleCellNet` & Seurat `MapQuery`
- **Source:** [`25_GBM_to_snMultiome_mapping.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/25_GBM_to_snMultiome_mapping.R) & [`24_GPC_to_snMultiome_mapping_and_Sankey_plot.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/24_GPC_to_snMultiome_mapping_and_Sankey_plot.R)
- **Rationale:** Gangliogliomas are driven by early developmental somatic mutations (e.g., *BRAF* V600E, *CD34* positivity). Projecting neoplastic cells onto normal developmental references determines whether they stall at a radial glia, intermediate progenitor (`IPC-Glia` / Tri-IPC), or immature astrocyte/neuronal state.
- **Implementation Strategy:**
  1. Train a `singleCellNet` Random Forest classifier on the Kriegstein human fetal atlas (and our own `GG001AC` normal cortex).
  2. Classify neoplastic cells from `GG001TB` and `GG001IF`, and render identity flow via `networkD3::sankeyNetwork` river plots.

### 4. TF Motif Accessibility Profiling with `JASPAR2020` + `chromVAR`
- **Source:** [`18_Motif_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/18_Motif_analysis.R)
- **Rationale:** In Ganglioglioma, hyperactivation of the MAPK/ERK pathway alters the binding of AP-1 complex factors (*FOS*, *JUN*), neurogenic factors (*NEUROD1*, *ASCL1*), and glial master regulators (*SOX9*, *OLIG2*).
- **Implementation Strategy:**
  1. Attach `JASPAR2020` position weight matrices to the unified `ChromatinAssay` via `Signac::AddMotifs()`.
  2. Compute single-cell motif accessibility deviation z-scores via `Signac::RunChromVAR()`, creating a `chromvar` assay for differential TF activity across `TB`, `IF`, `AMY`, and `AC`.

### 5. Spatial-Regional Compositional Modeling & Pseudobulk Differential Expression
- **Source:** [`03_snMultiome_proportion_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/03_snMultiome_proportion_analysis.R) & [`09_L4_IT_subtypes.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/09_L4_IT_subtypes.R)
- **Rationale:** Comparing tumor core against the epileptogenic interface and access cortex requires statistical models that account for biological pseudoreplication and compositional constraints.
- **Implementation Strategy:**
  1. Run `speckle::getTransformedProps()` + `limma` to evaluate cell-type proportion shifts across `AC -> IF -> TB`.
  2. Perform pseudobulk aggregation via `scuttle::aggregateAcrossCells()`, and fit mixed-effects models (`edgeR` + `glmmSeq` or `limma-voom`) to extract robust differential expression and accessibility signatures.

### 6. 8D WNN Lineage Trajectory (`Slingshot`) & Dynamic GAMs (`tradeSeq`)
- **Source:** [`07_slingshot.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/07_slingshot.R) & [`08_tradeSeq.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/08_tradeSeq.R)
- **Rationale:** Resolving whether neoplastic astrocytes and dysplastic ganglion cells diverge from a shared neoplastic progenitor or arise from separate clones requires continuous pseudotime trajectory inference.
- **Implementation Strategy:**
  1. Compute an 8-dimensional WNN UMAP embedding on the tumor glioneuronal compartment.
  2. Fit lineage splines with `slingshot`, and model branch-point driver genes using `tradeSeq::fitGAM()` with `earlyDETest()`.


---

## Part 4: Signac & SCENIC+ Usage in the Reference Study

Both **Signac** and **SCENIC+** were used in the reference pipeline, playing complementary roles:

1. **Signac (Core ATAC-seq & Motif Framework in R):**
   - [`00_snMultiome_data_integration.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/00_snMultiome_data_integration.R#L4) loads `library(Signac)` for peak quantification (`FeatureMatrix`), `ChromatinAssay` construction, ATAC QC (`NucleosomeSignal`, `TSSEnrichment`), LSI dimensionality reduction (`RunTFIDF`, `RunSVD`), and Reciprocal LSI (`rlsi`) batch integration.
   - [`18_Motif_analysis.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/18_Motif_analysis.R#L1-L19) uses `Signac::AddMotifs()` (with `JASPAR2020`) and `Signac::RunChromVAR()` to compute per-cell TF motif accessibility deviation z-scores.

2. **SCENIC+ (External eGRN Inference + `AUCell` Scoring in R):**
   - The study used **SCENIC+** (the joint RNA+ATAC multiome successor to SCENIC) rather than RNA-only SCENIC.
   - SCENIC+ network inference itself is **not** executed in the R scripts—it was run externally (in Python) to generate the enhancer-driven regulon (eGRN) gene-set file (`results/AUCell/582.sceplus.gmt`).
   - [`06_AUCell.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/06_AUCell.R#L18-L30) imports `582.sceplus.gmt`, scores every cell's raw RNA counts with `AUCell::AUCell_run()`, and stores the scores as a `"scenicplus_AUC"` assay in the Seurat object (which is later modeled along pseudotime in [`08_tradeSeq.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/08_tradeSeq.R) and correlated with trait relevance scores in [`21_Correlation_test_of_TRS_and_AUCell_along_the_pseudotime.R`](file:///Users/franz/Documents/Code/Ganglioglioma/Human_Cortex_Dev_Multiome-main/RScript/21_Correlation_test_of_TRS_and_AUCell_along_the_pseudotime.R)).
