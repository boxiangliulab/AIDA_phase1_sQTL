# Plagiarism Cross-Check Report

## AIDA_phase1_sQTL vs. sc-splice (https://github.com/liangy79917/sc-splice)

**Date:** 2026-02-17
**Method:** Automated line-by-line code comparison across all scripts in both repositories

---

## Executive Summary

Extensive code-level comparison reveals **strong evidence that sc-splice is derived from AIDA_phase1_sQTL**. The sc-splice repository appears to be a simplified, templatized adaptation of AIDA_phase1_sQTL's analysis pipeline and figure-generation code, applied to a different dataset (OneK1K instead of AIDA). Evidence includes:

- **Near-verbatim code blocks** in at least 6 scripts
- **20+ distinctive variable names** shared identically across repos
- A **shared misspelling** (`palatte` instead of `palette`) preserved in both repos
- **Character-identical parameter strings** in pipeline commands
- **Identical custom R script names** (e.g., `runFDR.mappability.R`)
- **Identical helper functions** copied wholesale
- **Same output filenames** (e.g., `Fig.5a.pdf`)

The sc-splice repository strips away AIDA_phase1_sQTL's production infrastructure (loops, job scheduling, logging) and replaces dataset-specific values (cell type names, file paths, sample sizes) while preserving the core analytical code, variable names, and algorithmic logic.

---

## 1. Repository Overview

| Aspect | AIDA_phase1_sQTL | sc-splice |
|--------|-------------------|-----------|
| **Dataset** | AIDA (Asian Immune Diversity Atlas) | OneK1K |
| **Cell types** | 19 (e.g., `CD14+_Monocyte`, `CD4+_T_naive`) | 13 (e.g., `MonocM`, `CD4Naive`) |
| **Languages** | R (64 files), Shell (19), Python (1) | R + Shell (~40 files) |
| **Publication** | Tian et al., Nat Genet 56, 2739-2752 (2024) | Unknown |
| **Structure** | `Analysis_pipeline/` + `Figures/` | `analysis_code/` + `Figure_code/` |

---

## 2. Evidence of Plagiarism by Category

### 2.1 Near-Verbatim Code Blocks (Highest Severity)

#### (A) UpSet Plot (trans-sQTL Figure 5A)
- **AIDA**: `Figures/Figure5/Figure5A.R`
- **sc-splice**: `Figure_code/Fig 5/fig 5a.R`

The entire UpSet plot construction is character-for-character identical, including:
- The misspelled variable name `palatte` (instead of `palette`) in both repos
- Identical complex `UpSet()` call with `upset_right_annotation()` and `upset_top_annotation()`
- Same variable names: `lt`, `m1`, `m`, `ss`, `cs`, `setOrder.manual`
- Same output filename: `'Fig.5a.pdf'`
- Same dimensions: `width = 6, height = 5`

```r
# Both repos contain this exact code:
palatte <- color$V2
names(palatte) <- color$V1
...
UpSet(m,
      set_order = match(setOrder.manual, set_name(m)),
      right_annotation = upset_right_annotation(m, gp = gpar(col = NA, fill = palatte[names(ss)])),
      comb_order = order(comb_degree(m), -cs),
      top_annotation = upset_top_annotation(m, gp = gpar(col = NA, fill = c(rep('black', length(grep(1, comb_degree(m), invert = T))), palatte[set_name(m)]))),
      comb_col = '#00004d'
)
```

#### (B) PTPRC and CTLA4/CD44 Heatmaps (Figure 2)
- **AIDA**: `Figures/Figure2/Figure2B.R`, `Figure2C.R`
- **sc-splice**: `Figure_code/Fig 2/fig 2f.R`, `fig 2g.R`

Verbatim copy with only data frame name and output filename changed:

```r
# Both repos (only data frame name differs):
f1 <- ggplot(DATA, aes(x=iid_label, y=celltype, fill=average_normalized_splicing)) + geom_tile() +
  scale_fill_viridis_c() +
  theme(legend.position = 'top',
        axis.text.y = element_text(colour = celltype_pal[as.vector(ord_level2)]),
        axis.text.x = element_text(colour = ifelse(ord_level %in% iid_label[cluster.sig], 'red', 'grey'), angle = 45, hjust = 1))
```

Shared distinctive variable names: `iid_label`, `celltype_pal`, `ord_level2`, `ord_level`, `cluster.sig`, `average_normalized_splicing`

#### (C) sQTL Position Distribution Algorithm (Figure 3B)
- **AIDA**: `Figures/Figure3/Figure3BC.R`
- **sc-splice**: `Figure_code/Fig 3/fig 3b.R`

Identical algorithm with identical variable names (`distance`, `frequency`, `pro`), identical matrix dimensions `matrix(0,110,2)`, and a character-for-character identical labeling loop:

```r
# Identical in both repos:
for(i in 1:50){
  frequency[51-i,2]<-paste(i,"*0.5kb upstream",sep="")
  frequency[60+i,2]<-paste(i,"*0.5kb downstream",sep="")
}
for(i in 1:10){
  frequency[50+i,2]<-paste(i," internal",sep="")
}
colnames(frequency)<-c("Number of sQTL","Position")
```

#### (D) Colocalization Core Analysis
- **AIDA**: `Analysis_pipeline/06_colocalization/01.colocalization_process.sh` (lines 139-235)
- **sc-splice**: `analysis_code/9 Colocalization.sh`

The entire `coloc.abf()` block (~12 lines) is near-identical, including:
- Identical guard: `if(nrow(input)>10){`
- Identical trait-type branching on `trait_type=="cc"` and `trait_type=="quant"`
- Identical dataset specifications with same variable names
- Identical result extraction: `t(as.data.frame(result$summary))[1,1:6]`
- Identical variable names throughout: `test_rs`, `snp_1mb`, `input`, `GWAS_pre`, `MAF`, `num_1`, `num_2`, `num_3`, `trait_type`
- Same debug print: `print("start read in")`
- Same variance calculation: `input$varbeta<-(input$se)^2`
- Same merge logic (character-identical): `merge(snp_1mb,MAF,by="variant_id")` and `merge(input,GWAS_pre,by="variant_id",suffixes=c("_sqtl","_gwas"))`

Note: sc-splice fixes a likely bug in AIDA's cc-trait `N=` parameter (`quant_trait_num[num_1]` -> `cc_trait_case[num_1]+cc_trait_control[num_1]`).

#### (E) Locuscompare Helper Functions (Figure 5F)
- **AIDA**: `Figures/Figure5/Figure5F.R`
- **sc-splice**: `Figure_code/Fig 5/fig 5f.R`

Two helper functions copied verbatim:

```r
# Identical in both repos:
sub <- function(row) {
        chr <- as.character(row["VID"])
        chr <- substr(chr, 1, nchar(chr) - 2)
        return(chr)
}
add_label <- function(merged, snp){
    merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, '')
    return(merged)
}
```

#### (F) Library Loading in Colocalization
- **AIDA**: `Analysis_pipeline/06_colocalization/01.colocalization_process.sh`
- **sc-splice**: `analysis_code/9 Colocalization.sh`

Character-identical, including inconsistent quoting:
```r
library("coloc")
library(dplyr)
library(data.table)
```

---

### 2.2 Identical Pipeline Commands (High Severity)

#### (A) STAR Alignment Parameters
- **AIDA**: `Analysis_pipeline/00_preprocessing/02_run_STARsolo.sh`
- **sc-splice**: `analysis_code/1 Demultiplexing sample.sh`

The 17-attribute `--outSAMattributes` string is character-for-character identical:
```
NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW
```
Additional identical parameters: `--twopassMode Basic`, `--soloStrand Forward`, `--soloType CB_UMI_Simple`, `--soloCBmatchWLtype 1MM`, `--soloUMIdedup 1MM_Directional_UMItools`, `--waspOutputMode SAMtag`, whitelist `737K-august-2016.txt`.

#### (B) cis-sQTL Calling (QTLtools)
- **AIDA**: `Analysis_pipeline/03_cis-sQTL/01.phenotype_prepare&&nominal_pass.sh`, `02.permutation_conditional_pass.sh`
- **sc-splice**: `analysis_code/6 Call cis-sQTL.sh`

All three QTLtools cis invocations (permutation, conditional, nominal) use identical flag sets in identical order. The nominal pass line-break formatting is distinctive and identical:
```bash
--nominal 1 --normal \
--grp-best --out ...
```
Same FDR script (`runFDR_cis.R`) with same threshold (`0.05`). Same file naming patterns throughout.

#### (C) trans-sQTL Calling (QTLtools)
- **AIDA**: `Analysis_pipeline/05_trans-sQTL/` (3 shell scripts)
- **sc-splice**: `analysis_code/8 Call trans-sQTL.sh`

Character-identical parameter strings:
```
--threshold 1e-5 --normal --window 5000000 --bin 1000
--threshold 0.1 --normal --window 5000000 --bin 1000
```
Same custom R script name: `runFDR.mappability.R` with same calling convention and FDR threshold.

#### (D) Picard MarkDuplicates
Both repos use identical parameters: `REMOVE_DUPLICATES=true`, `BARCODE_TAG=CB`, `MOLECULAR_IDENTIFIER_TAG=UB`, same metrics file pattern `marked_dup_metrics.SAMPLE.txt`.

#### (E) Pseudobulk BAM Processing
- **AIDA**: `Analysis_pipeline/00_preprocessing/pseudobulk_bam.py`
- **sc-splice**: `analysis_code/4 Pseudobulk.sh` (actually Python, mislabeled)

Identical argument names (`--bam_file`, `--csv_file`, `--out_prefix`), identical filtering logic (WASP vW==1 with NH==1 fallback), identical output naming (`PREFIX.CELLTYPE.bam`). sc-splice is a refactored OOP version of the same algorithm.

---

### 2.3 Shared Distinctive Variable Names (High Severity)

The following project-specific variable names appear identically in both repos:

| Variable Name | Context | Notes |
|---------------|---------|-------|
| `palatte` | UpSet plot | Misspelling preserved in both |
| `setOrder.manual` | UpSet plot | Custom ordering variable |
| `iid_label` | Heatmap plots | Individual ID labels |
| `average_normalized_splicing` | Heatmap plots | Column name in data frame |
| `celltype_pal` | Heatmap plots | Cell type palette |
| `ord_level`, `ord_level2` | Heatmap plots | Ordering variables |
| `cluster.sig` | Heatmap plots | Significant clusters |
| `plot_data_1`, `plot_data_2`, `plot_data_3` | sQTL summary | Sequential processing |
| `cell_type_order` | sQTL summary | Cell type ordering |
| `norm_const` | sQTL summary | Normalization constant |
| `custom_sqtl_color` | sQTL summary | Color palette name |
| `multi_sqtl_gene` | Heterogeneity | Multi-sQTL gene data |
| `total_sgene`, `prop_multi_sqtl_gene` | Heterogeneity | Column names |
| `gene_coord`, `celltype_coord` | Scatterpie plots | Coordinate variables |
| `percent.of.colocalized.GWAS.loci` | Coloc heatmap | Long distinctive column name |
| `number.of.colocalization` | Coloc heatmap | Column name |
| `test_rs`, `snp_1mb` | Colocalization | Analysis variables |
| `num_1`, `num_2`, `num_3`, `trait_type` | Colocalization | Index/type variables |
| `label_dict` | Multiple figures | Cell type renaming |

---

### 2.4 Structural Similarities (Moderate Severity)

1. **Repository organization**: Both split code into analysis pipeline and figure generation directories
2. **Figure numbering**: sc-splice figures correspond to AIDA figures (Fig 2 -> Fig 2, Fig 3 -> Fig 3, Fig 5 -> Fig 5, Fig 6 -> Fig 4)
3. **Analysis workflow**: Identical pipeline steps in same order (STAR -> dedup -> pseudobulk -> leafcutter -> cis-sQTL -> trans-sQTL -> colocalization)
4. **`label_dict` + `for(i in seq_along(label_dict))` pattern**: Same cell-type renaming idiom used throughout both repos

---

## 3. Key Differences

| Aspect | AIDA_phase1_sQTL | sc-splice |
|--------|-------------------|-----------|
| **Dataset** | AIDA (multi-ethnic Asian cohort) | OneK1K |
| **Cell types** | 19 types, long names | 13 types, abbreviated names |
| **Color palette** | Project-specific hex codes | Different project-specific hex codes |
| **Production infrastructure** | Full loops, job scheduling, logging | Template/placeholder style, no loops |
| **Computing environment** | NUS HPC + AWS (`/ebs1/`, `/home/users/nus/`) | Different cluster (`/public/home/liangyan/`) |
| **Unique analyses in AIDA** | Dynamic sQTL, sex-biased sQTL, LDSC, SMR, PIP fine-mapping, spliZ | Not present |
| **Unique analyses in sc-splice** | Age-stratified sQTL, eQTL (TensorQTL), DoubletFinder, Seurat/Harmony | Not present |
| **Some scripts in sc-splice** | - | Appear independently written (Seurat, DoubletFinder, violin plot template) |

---

## 4. Conclusion

The evidence strongly indicates that **sc-splice was derived from AIDA_phase1_sQTL**. The pattern is consistent: sc-splice takes AIDA's analytical code, strips production infrastructure (loops, job management, detailed paths), replaces dataset-specific values (cell types, sample sizes, file paths), and in some cases refactors the code style (e.g., OOP pseudobulk script). However, core algorithmic logic, variable naming, function implementations, parameter choices, and even typos are preserved from the original.

The strongest individual pieces of evidence are:
1. The **shared `palatte` misspelling** in the UpSet plot code
2. The **verbatim PTPRC/CTLA4 heatmap code** with 6+ identical distinctive variable names
3. The **character-identical `coloc.abf()` block** (~12 lines) with a bug fix applied by sc-splice
4. The **identical 17-attribute `--outSAMattributes` STAR string**
5. The **identical `runFDR.mappability.R` custom script name** and calling convention

These cannot be explained by independent implementation following tool documentation or common coding conventions.
