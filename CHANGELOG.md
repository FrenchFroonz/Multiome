# Multiome Pipeline Change Log

All notable changes, architectural updates, and code modifications for this repository are documented in this file.

---

## [Current - 2026-10-13]
### Added
- **`SKILLS_BREAKDOWN_AND_RECOMMENDATIONS.md`**: Detailed technical synthesis of newly added `bio-atac-seq-*` skills, comparison with existing pipeline, and scientific rationale for optimizations.
- **`TODO.md`**: Project roadmap defining upcoming scripts (03: TF Motifs & chromVAR, 04: Enhancer-Gene Linking, 05: Pseudobulk DARs, 06: Footprinting) and cluster pseudobulk peak calling.
- **`CHANGELOG.md`**: Chronological log of all major project modifications.
- **TSS Enrichment Distribution Analysis in `ATAC_Ganglioglioma.Rmd`**: Added density and distribution analysis chunk with threshold lines at 2, 3, and 4 (saved to `outputs/TSS_enrichment_distribution.pdf`) to empirically assess cutoff stringency for brain tumor/interface samples.
- **Single-Cell ATAC Noise Metrics in `ATAC_Ganglioglioma.Rmd`**: Added calculation of fraction of reads in peaks (`pct_reads_in_peaks` / FRiP) and ENCODE blacklist ratio (`blacklist_ratio`).
- **Gene Activity Matrix in `ATAC_Ganglioglioma.Rmd`**: Derived chromatin gene activity assay (`GeneActivity()`) as an orthogonal confirmation of cell type identity and chromatin-RNA concordance.
- **Corces et al. 2018 501 bp Fixed-Width Peaks in `02_Peak_Calling_MACS2.Rmd`**: Standardized MACS2 narrow peaks around summits ($\pm 250\text{ bp}$) to eliminate peak length confounding in feature counting and downstream statistical tests.
- **Statistical Significance ($q$-value) Filtering in `02_Peak_Calling_MACS2.Rmd`**: Filtered raw MACS2 calls by `neg_log10qvalue >= 2` ($q \le 0.01$).
- **Pseudobulk Peak Calling Documentation in `02_Peak_Calling_MACS2.Rmd`**: Added in-line documentation and TODO specifying future transition to per-cluster pseudobulk peak calling.

### Changed
- **Calibrated Nucleosome Signal Filter in `ATAC_Ganglioglioma.Rmd`**: Relaxed threshold from `< 2` to `< 4` following canonical Signac and ENCODE single-cell benchmarks to avoid discarding nuclei with intact nucleosomal structure.

---

## [2026-10-12]
### Fixed
- **Assay5 Fragments Dispatch in `02_Peak_Calling_MACS2.Rmd` (`59be756`)**: Fixed `UseMethod("Fragments")` failure in Seurat v5 by setting `DefaultAssay(ATAC_gang) <- "ATAC"` and explicitly targeting `Fragments(ATAC_gang[["ATAC"]])` instead of the top-level object containing the default `Assay5` RNA object.

### Added
- **Repository Documentation (`c4d16f0`, `419999f`)**: Created comprehensive `README.md` containing architecture flowcharts, file registry, cluster operating instructions, and git push history.

### Changed
- **Singularity MACS2 Output Integration (`e800be4`)**: Adapted `02_Peak_Calling_MACS2.Rmd` to ingest `.narrowPeak` output from external cluster Singularity run (`clean_fragments.bed` + MACS2), reconstruct chromatin assay via `FeatureMatrix()`, and recompute LSI and WNN embeddings.
- **Seurat Object Renaming (`23b98ea`)**: Replaced generic variable name `pbmc` with sample-specific name `ATAC_gang` across `ATAC_Ganglioglioma.Rmd` and `02_Peak_Calling_MACS2.Rmd`.

---

## [2026-10-08]
### Added
- **Automated Figure PDF Export (`7ea4b91`)**: Added PDF rendering of unannotated and scType-annotated multi-panel UMAPs directly to `outputs/`.
- **RDS Export & Step 2 Draft (`13385e1`)**: Added serialization of `Multiome_WNN_scType_Annotated.rds` and initial draft of `02_Peak_Calling_MACS2.Rmd`.
- **ATAC Quality Metrics (`de26b7a`)**: Integrated `NucleosomeSignal()` and `TSSEnrichment()` calculations into initial QC workflow.
- **Repository Initialization (`8ff65bb`)**: Set up Git tracking and `.gitignore` file rules for multiome workspace.
