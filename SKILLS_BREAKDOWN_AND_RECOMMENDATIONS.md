# ATAC-seq Skills Breakdown & Pipeline Recommendations

This document synthesizes best practices from the newly integrated ATAC-seq agent skills suite (`bio-atac-seq-*`) and compares them against our current Ganglioglioma single-cell multiome implementation (`ATAC_Ganglioglioma.Rmd` and `02_Peak_Calling_MACS2.Rmd`).

---

## 1. Relevant Skills Evaluated

The analysis specifically leveraged the following core skills from the `bio-atac-seq` suite:
* **`bio-atac-seq-single-cell-atac`**: Best practices for 10x scATAC and Multiome processing in Signac and ArchR, per-cell QC metrics, LSI failure modes, and pseudobulk peak calling.
* **`bio-atac-seq-atac-qc`**: Canonical library and per-cell QC thresholds, ENCODE 4 standards, TSS enrichment methodologies, and FRiP calculations.
* **`bio-atac-seq-atac-peak-calling`**: Shift-extend Tn5 mechanics, `-f BED` vs `-f BAMPE`, genome size calibrations, and MACS2/MACS3 algorithmic parameters.
* **`bio-atac-seq-consensus-peakset`**: Rationale and algorithms for building fixed-width consensus peaksets (Corces et al. 2018 standard) vs variable-width peak bias.
* **`bio-atac-seq-motif-deviation`**, **`bio-atac-seq-enhancer-gene-linking`**, and **`bio-atac-seq-differential-accessibility`**: Downstream module roadmaps.

---

## 2. What Was Already Implemented Correctly

Several key design choices in our existing scripts adhere strictly to canonical standards:

1. **Exclusion of LSI Component 1 for ATAC Clustering and UMAP:**
   * *Implementation:* We specified `dims = 2:50` for ATAC UMAP and in multimodal WNN neighbors.
   * *Skill Validation:* `bio-atac-seq-single-cell-atac` identifies component 1 as strongly confounded by sequencing depth (cellular read coverage) rather than chromatin biology. Skipping component 1 avoids density gradients driven purely by coverage differences.
2. **Canonical MACS2 Tn5 Transposition Cut-Site Recipe:**
   * *Implementation:* The cluster command `macs2 callpeak -t clean_fragments.bed -f BED --nomodel --shift -75 --extsize 150 -g hs` was executed.
   * *Skill Validation:* `bio-atac-seq-atac-peak-calling` highlights this exact configuration ("Pattern A" / ENCODE 4 standard). Because Tn5 inserts as a dimer creating a 9 bp stagger, centering a 150 bp window with a -75 bp shift around each insertion point accurately models the open nucleosome-free region (NFR) cut event.
3. **Genomic Range Cleaning:**
   * *Implementation:* Restricting features to standard chromosomes (`keepStandardChromosomes`) and filtering against ENCODE blacklist regions (`blacklist_hg38_unified`).
   * *Skill Validation:* Prevents artefactual, artificially high-signal regions (such as pericentromeric repeats and rDNA clusters) from skewing TF-IDF and SVD.
4. **Multimodal WNN Graph Construction:**
   * *Implementation:* Joint neighbor weighting (`FindMultiModalNeighbors`) using RNA PCA and ATAC LSI, clustering on `wsnn` using SLM algorithm.
   * *Skill Validation:* Properly pairs modalities at the single-nucleus level rather than analyzing them in isolated silos.

---

## 3. Recommended Modifications for Script 1: `ATAC_Ganglioglioma.Rmd`

### A. Calibrate Nucleosome Signal Threshold
* **Current Value:** `nucleosome_signal < 2`
* **Recommended Value:** `nucleosome_signal < 4`
* **Rationale:** The nucleosome signal represents the ratio of mononucleosomal fragments (147–294 bp) to nucleosome-free fragments (< 147 bp). Signac and ENCODE benchmarks establish that `nucleosome_signal <= 4` represents high-quality libraries with clear periodicity. A threshold of `< 2` is excessively stringent and risks prematurely discarding viable nuclei that display healthy, canonical nucleosomal wrapping.

### B. Investigate TSS Enrichment Distribution for Brain Tumor Samples
* **Current Value:** `TSS.enrichment > 2`
* **Context:** In general single-cell ATAC benchmarks, `TSS.enrichment >= 4` is the standard pass threshold, whereas `2–4` is considered "caution" / intermediate.
* **Rationale for Ganglioglioma Data:** Ganglioglioma samples consist of dysplastic neurons, neoplastic glial cells, reactive astrogliosis, and infiltrate from the brain tumor interface. Tumor and clinical resection tissue frequently exhibit higher background noise and lower baseline chromatin accessibility than cell lines or healthy PBMCs.
* **Implementation:** Instead of blindly increasing the cutoff to 4 (which could eliminate vulnerable tumor subpopulations), we add distribution visualization chunks (density curve, ECDF, and correlation with read count) with reference lines at 2, 3, and 4 to empirically calibrate the threshold.

### C. Add Missing Key Single-Cell ATAC QC Metrics
* **Metrics Added:**
  1. `pct_reads_in_peaks` (FRiP per cell): Measures the fraction of fragments in each cell that fall within accessible chromatin peaks versus random background genomic cleavage. Ideal target is `> 15%`.
  2. `blacklist_ratio`: Measures the fraction of counts falling within ENCODE blacklisted regions via `FractionCountsInRegion()`. High values (`> 5%`) flag apoptotic nuclei or non-specific transposed background.
* **Rationale:** Fragment counts alone do not detect dying nuclei with leaky, non-specific transposition. These two ratios clean the feature space prior to dimensionality reduction.

### D. Generate Gene Activity Assay (`GeneActivity`)
* **Rationale:** Although our pipeline utilizes paired snRNA-seq to annotate cell types via `scType`, deriving an ATAC-based `GeneActivity` matrix (aggregating open chromatin across gene bodies and promoters) enables direct validation of chromatin-RNA concordance. Comparing gene activity versus RNA expression for canonical brain/glioma markers (e.g., *GFAP*, *SOX2*, *MBP*, *RBFOX3*, *SLC1A2*) confirms that our epigenomic assay accurately reflects transcriptional identity.

---

## 4. Recommended Modifications for Script 2: `02_Peak_Calling_MACS2.Rmd`

### Proposition 1: Fixed-Width 501 bp Peaks (Corces et al. 2018 Standard)
* **The Problem:** Raw MACS2 narrow peaks vary widely in length (from 150 bp to several kilobases). When quantifying fragment counts with `FeatureMatrix()`, wider peaks capture more reads simply because they cover more genomic territory, not because accessibility is biologically higher. This introduces severe length-dependent bias into TF-IDF, clustering, differential accessibility, and motif scanning.
* **The Solution:**
  1. Calculate the true summit coordinate for each MACS2 peak:
     $$\text{summit} = \text{start} + \text{relative\_summit\_pos}$$
  2. Re-center each peak on its summit and extend $\pm 250\text{ bp}$ to generate uniform 501 bp fixed-width regions:
     $$\text{fixed\_start} = \text{summit} - 250, \quad \text{fixed\_end} = \text{summit} + 250 + 1$$
* **Rationale:** Adopts the Corces et al. (2018) / ArchR standard, standardizing genomic windows and eliminating peak-width confounding.

### Proposition 2: Filter by Statistical Significance ($q$-value)
* **Implementation:** Filter imported MACS2 peaks by `neg_log10qvalue >= 2` ($q \le 0.01$, FDR 1%).
* **Rationale:** Bulk peak calling identifies tens of thousands of low-significance background peaks. Removing marginal peaks before quantifying the feature matrix reduces sparsity and computational overhead while retaining genuine regulatory elements.

### Proposition 3: Per-Cluster Pseudobulk vs Bulk Peak Calling (Future Enhancement)
* **Concept:** Currently, MACS2 peak calling was performed on bulk fragments pooled across all barcodes.
* **Skill Recommendation:** Bulk peak calling is dominated by high-frequency cells (e.g., astrocytic tumor cells) and fails to identify regulatory elements active exclusively in low-abundance cell types (e.g., tumor-associated microglia, endothelial cells, or rare neuronal populations).
* **Future Implementation:** Split fragment files by cell type (`sctype_classification`), run MACS2 independently on each pseudobulk subpopulation, and merge results using iterative overlap removal into a comprehensive consensus peakset. Documented in the script and tracked in `TODO.md`.

---

## 5. Summary Matrix of Pipeline Improvements

| Pipeline Step | Previous Setup | Proposed / Updated Setup | Target Benefit |
| :--- | :--- | :--- | :--- |
| **Nucleosome Signal** | Cutoff `< 2` | Cutoff `< 4` | Retains valid single nuclei with normal nucleosome periodicity. |
| **TSS Enrichment** | Hard cutoff `> 2` | Empirical distribution check (2, 3, 4) | Prevents accidental loss of sensitive brain tumor / interface cells. |
| **Noise Filtering** | Fragment & RNA counts only | Added `blacklist_ratio` & `pct_reads_in_peaks` | Eliminates apoptotic nuclei with non-specific Tn5 transposition. |
| **Modality Validation** | RNA expression only | Added `GeneActivity` matrix | Verifies chromatin-to-expression concordance for key marker genes. |
| **Peak Width** | Variable width (MACS2 default) | Fixed 501 bp centered on summit ($\pm 250$ bp) | Removes peak length bias in quantification, LSI, and differential tests. |
| **Peak Significance** | All called peaks | Filter `neg_log10qvalue >= 2` ($q \le 0.01$) | Eliminates marginal false-positive background peaks. |
| **Peak Scope** | Bulk fragment pool | Bulk fixed-width now; Per-cluster pseudobulk in TODO | Ensures modularity today while providing a roadmap for rare cell types. |
