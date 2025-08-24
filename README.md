# Single-Cell: Endothelial and Immune Remodeling in IPF and SSc-ILD Lungs

## Overview

This repository contains the implementation of a single-cell RNA-seq analysis of the endothelial (CD31+) and immune (CD45+) compartments of fibrotic human lung. Starting from the sorted-fraction GEO dataset GSE159354 (control, idiopathic pulmonary fibrosis and systemic-sclerosis-associated interstitial lung disease), it integrates the two compartments across ten donors, resolves the endothelium into arterial, venous, capillary and lymphatic subtypes, subclusters the 3,532 lymphatic endothelial cells (LECs) into eight states, follows them along pseudotime, tests them for differential expression and GO enrichment at the donor level, and asks which LEC states signal to which immune populations with CellChat. The analysis identifies a disease-enriched LEC subpopulation that has lost lymphatic marker expression and gained hypoxia-inducible, venous and proinflammatory programs, and that signals most strongly to M1-like macrophages.

## Research Background

### Lymphatic Endothelium in Fibrotic Lung Disease

Pulmonary fibrosis is a devastating condition that affects more than 5 million people globally; idiopathic pulmonary fibrosis (IPF) and interstitial lung disease associated with autoimmune disorders, including SSc-ILD, are its most common types, and median survival after diagnosis is only 2-5 years. Immune abnormality in pulmonary fibrosis is characterized by aberrant activation of both innate and adaptive immune cells, yet high-resolution studies of immune-stromal crosstalk in human disease are scarce. The lymphatic vasculature maintains local fluid homeostasis and regulates immune trafficking, and emerging data suggest it undergoes pathological remodeling in pulmonary fibrosis; however, it is unclear whether altered lymphatic endothelial cells (LECs) in fibrotic lungs directly modulate immune cell activity and shape the inflammatory microenvironment.

### The Challenge

The goal of this study is to delineate phenotypic alterations in LECs and determine whether they directly modulate immune cell differentiation and activity in pulmonary fibrosis. The widely used bleomycin-induced mouse IPF model is characterized by an acute phase of immune activation followed by tissue fibrosis, which does not fully recapitulate human disease; thus, immune-related studies benefit greatly from analyses of human samples. Since LECs express the canonical EC marker CD31, they represent a subset of CD31+ ECs, and in single-cell datasets without enrichment of CD31+ ECs, LEC populations are generally small. We therefore identified and used a single available dataset which contained enriched CD31+ ECs among 4 cell types derived from 3 control, 4 IPF, and 3 scleroderma-associated (SSc)-ILD patients, with a total of 200,500 cells.

Sorting makes the EC and immune compartments abundant enough to study, but the sort is imperfect (each fraction carries cells of the other type, plus doublets), the four EC subtypes do not separate by clustering alone, and the 3,532 LECs are spread over ten donors, so cell-level tests would be pseudoreplicated.

### The Solution: Staged Cleaning, Module-Score Subtyping, Donor-Level Pseudobulk and Ligand-Receptor Mapping

This project addresses these problems with a staged pipeline. The two fractions are integrated with Harmony on donor and sort fraction, then cleaned in three passes (joint-cluster removal, doublet flagging, cross-compartment module-score filtering), each checked with the same single-cell marker-panel heatmap. ECs are assigned to arterial, venous, capillary or lymphatic identity by the highest of four `AddModuleScore` signatures, re-clustered, and filtered to the clusters that are bona fide ECs. LECs are re-clustered on their own into LEC0–LEC7; subclusters are labelled by their disease-group composition; Slingshot pseudotime is run through them; differential expression is done on donor × subcluster pseudobulks with DESeq2 and followed by GO enrichment; and the LEC states are joined to the 21 annotated immune clusters in CellChat.

## Research Context

We did this work as part of the Kao Lab at Stanford University. It is written up in *Hypoxia-Reprogrammed Lymphatic Endothelial Cells Promote Pathogenic Immunity in Lung Fibrosis* (Xia & Wang), and the numbers quoted under Key Results are taken from that manuscript. We developed the analysis interactively in one R notebook (kept, lightly cleaned, in `notebooks/`) and later re-organised it into the numbered stages in `scripts/` with the shared code in `R/`; `docs/changes.md` lists every edit we made in that re-organisation, including the bugs it fixed. The analysis is exploratory: we made the cluster-removal decisions and the LEC state labels by inspecting marker plots and composition tables, and they are recorded as parameters in the scripts and in `docs/pipeline.md` so they can be revisited.

## Dataset

**GEO GSE159354** – *Single-cell RNA-seq of sorted human lung explants from IPF and SSc-ILD patients* (10x Genomics Chromium). The series contains 48 sorted libraries from 10 lung explant donors: 3 normal controls (N1857, N1865, N1902), 4 IPF (F18536, F18545, F18555, F19609) and 3 SSc-ILD (S19578, S19584, S19600). Each donor was sorted into several fractions (EPCAM+, CD90+, CD31+, CD45+, triple-negative, unsorted, plus senescence- and HT2-280-split epithelial fractions for some donors).

This pipeline uses the two fractions that exist for every donor:

- **CD31+ (EC)** – endothelial cells, 10 samples (GSM4827161, 166, 171, 175, 179, 183, 188, 193, 198, 203)
- **CD45+ (Immune)** – immune cells, 10 samples (GSM4827162, 167, 172, 176, 180, 184, 189, 194, 199, 204)

Per-sample QC keeps genes detected in ≥ 3 cells and cells with 200 < genes < 8,000 and < 20 % mitochondrial reads. The full manifest, the loader and the QC thresholds are in `R/manifest.R`. The raw `.h5` files are not included; `data/README.md` says where to put them.

## Methodology

### Processing Pipeline

1. **Integration** (stage 01): the 20 samples are merged, log-normalised (scale 10⁴), reduced on 3,000 variable genes and 50 PCs, integrated with Harmony on donor and sort fraction, and clustered (Louvain, resolution 0.6) on Harmony dimensions 1–30.
2. **Cleaning** (stages 02–03): six joint clusters that were neither clean EC nor clean immune (7, 12, 15, 16, 17, 19) are removed; doublets are flagged as the top 5 % of z(genes) + z(UMIs); ECs above the 90th percentile of an immune module score and immune cells above the 90th percentile of an EC module score are removed. The result is the "Clean" object the immune branch starts from. Two more aggressive variants were tried and kept under `scripts/exploratory/`.
3. **Quality control** (stage 04): after every step the same single-cell heatmap of 37 lymphatic / arterial / venous / core-EC / immune markers, cells ordered EC then immune, shows whether either compartment still carries the other's program (`R/marker_panel_heatmap.R`).
4. **Endothelial subtyping** (stages 10–12): ECs are re-processed (2,000 variable genes, Harmony on donor and group, resolution 0.8) and assigned to ArtECs, VenECs, CapECs or LECs by the highest of four module scores (arterial GJA5/EFNB2/DLL4/SOX17/HEY2/NOTCH4/BMX/SEMA3G; venous NR2F2/EPHB4/ACKR1/VCAM1/SELP/SELE/VWF/PLVAP; capillary CA4/EDNRB/EMCN/ESAM/KDR/RGCC/AQP1/HSPG2; lymphatic PROX1/LYVE1/PDPN/FLT4/VEGFC/CCL21/RELN/TBX1/FOXC2/NFATC1). EC clusters 0, 1, 2, 3, 6, 7, 10 and 20 are kept; a subtype × group heatmap of 30 zonation genes checks the assignment.

### Experiments

**E1: LEC subclustering** (stage 20) – LECs are re-processed on their own (2,000 variable genes, Harmony on donor and group) and clustered at resolution 0.7, giving LEC0–LEC7; composition by disease group, lymphatic marker feature plots.

**E2: LEC state signatures** (stages 21–23) – violin and feature plots of candidate state genes (SPARCL1, IL3RA, POSTN, KLF4, TBX1, EBF1, COL15A1, STC1, HIF1A, EDNRB, FLT4); dot plots of a hypoxia / inflammation (HIF1A, EPAS1, ICAM1, SELP, TNFSF4, CCL2, CXCL1, CXCL2), antigen-presentation (CD74, HLA-DPA1/DPB1/DQA1/DRB1) and lymphatic-identity (PROX1, FLT4, TEK, S1PR1, EBF1, MEOX1) signature across the subclusters, in two orders (LEC0–7 and the merged LEC0, LEC6, LEC1, LEC4, LEC2, LEC3/5/7); a 45-gene, 8-category heatmap of subcluster means.

**E3: Trajectory** (stage 24) – Slingshot on the LEC UMAP, automatic and guided (start LEC0, ends LEC6 and LEC7), with Monocle3 as an optional second method. Subclusters are labelled by composition as Normal (LEC1, LEC6; control-enriched), Proinflammatory (LEC0, LEC2, LEC3, LEC5; disease-enriched) or Intermediate (LEC4, LEC7), and PROX1, EPAS1, HIF1A and SPARCL1 are followed along pseudotime by state and by disease group.

**E4: Differential expression** (stages 25–26) – donor-level pseudobulk DESeq2. Within LEC1 and within LEC3: IPF vs CTRL, SSc-ILD vs CTRL, IPF vs SSc-ILD (`~ group`). Between subclusters, with the homeostatic LEC1 as reference: LEC1 vs LEC3, LEC1 vs LEC6, LEC1 vs LEC5 (`~ cluster`), each run on all donors and on IPF + SSc-ILD donors only, with two-sided volcano plots (`R/pseudobulk_de.R`).

**E5: Pathway enrichment** (stage 27) – `clusterProfiler::enrichGO` (BP / MF / CC, org.Hs.eg.db) on the padj < 0.05, |log2FC| > 0.25 gene sets of LEC1 vs LEC3 and LEC1 vs LEC6, against the universe of tested genes, with the top biological-process terms compared across the four gene sets of each comparison (`R/go_enrichment.R`).

**E6: Immune compartment** (stages 30–33) – immune cells from the Clean object are re-processed (3,000 variable genes, Harmony) and clustered at resolution 0.8; `FindAllMarkers` (only.pos, min.pct 0.25, logFC 0.25); clusters 13 (epithelial contamination), 20 (unknown), 22 (proliferating B cells) and 24 (mast / basophils) are removed and the remaining 21 are renumbered 0–20 and annotated from their top-10 markers as lymphoid (T, NK, B, Treg, proliferating, cytotoxic T), myeloid (macrophage, monocyte, M2-like, interferon-responsive, lipid-associated) and mixed (neutrophil, mast, plasma, dendritic, tissue-resident CTL, plasmacytoid DC) sets, with marker heatmaps for each set.

**E7: LEC–immune communication** (stage 40) – CellChat (CellChatDB.human, raw counts on shared genes, min. 10 cells) between all eight LEC subclusters and the 21 immune clusters ("Complete"), and with the disease-enriched LEC0, LEC3 and LEC5 as senders ("Targeted"); interaction counts and strengths, top pathways, and LEC → immune signalling matrices (`R/cellchat_lec_immune.R`).

## Key Results

### Endothelial Atlas

From the 200,500-cell dataset (3 control, 4 IPF, 3 SSc-ILD donors), 24,329 ECs passed QC: 2,825 control, 13,066 IPF and 8,438 SSc-ILD. Module-score subtyping resolved them into four vascular beds:

| Subtype | Cells | Markers used |
|---|---|---|
| Arterial (ArtECs) | 3,118 | EFNB2, GJA5, SEMA3G |
| Capillary (CapECs) | 6,114 | EMCN, ESAM, RGCC |
| Venous (VenECs) | 11,565 | ACKR1, NR2F2, SELP, SELE |
| Lymphatic (LECs) | 3,532 | PROX1, PDPN, LYVE1, CCL21, FLT4 |

On the joint UMAP, IPF and SSc-ILD ECs shift away from control and resemble each other. Within the CD31+ compartment the LEC share is higher in IPF and lower in SSc-ILD than in control, in line with earlier reports of enhanced lymphangiogenesis in IPF and reduced lymphangiogenesis in SSc-ILD.

### LEC Subclusters

The 3,532 LECs (1,146 control, 1,548 IPF, 838 SSc-ILD) re-cluster at resolution 0.7 into eight subsets, numbered by size. LECs from fibrotic lungs shift position on the LEC-only UMAP and show reduced PROX1, PDPN and FLT4 density.

| Subcluster | CTRL | IPF | SSc-ILD | % disease | Role |
|---|---|---|---|---|---|
| LEC0 | 289 | 516 | 222 | 72 % | trajectory root; normal-like sender in CellChat |
| LEC1 | 527 | 104 | 83  | 26 % | normal LEC, end of trajectory 1; DE reference |
| LEC2 | 106 | 322 | 216 | 84 % | trajectory 2 |
| LEC3 | 21  | 269 | 139 | 95 % | hypoxia-reprogrammed LEC, end of trajectory 2 |
| LEC4 | 93  | 135 | 80  | 70 % | not on either trajectory |
| LEC5 | 9   | 151 | 50  | 96 % | intermediate, trajectory 2 |
| LEC6 | 100 | 29  | 44  | 42 % | normal LEC |
| LEC7 | 1   | 22  | 4   | 96 % | disease-only, very small |

LEC3, LEC5 and LEC7 are ≥ 95 % disease-derived and together make up nearly 20 % of all LECs. LEC1 is the control-dominated, homeostatic state and is the reference for every between-subcluster comparison. (The composition-based labels used inside stage 24 – Normal / Proinflammatory / Intermediate – are coarser than the roles above; see `docs/pipeline.md`.)

### Trajectories

Slingshot infers two paths out of LEC0: **LEC0 → LEC1** (trajectory 1; LEC1 is 74 % control, so this is read as normal lymphangiogenesis) and **LEC0 → LEC2 → LEC5 → LEC3** (trajectory 2; LEC3 is 95 % disease, read as pathological LEC reprogramming). Along trajectory 2, PROX1 and FLT4 fall while EPAS1 and HIF1A rise. In the 8-category signature heatmap, LEC3, LEC5 and LEC7 show reduced lymphatic markers (PDPN, LYVE1, MMRN1, FLT4), elevated hypoxia-associated factors (EPAS1, HIF1A) and the highest MHC class II expression (HLA-DPA1, HLA-DPB1, HLA-DRA, HLA-DRB1, HLA-DQA1, CD74).

### Differential Expression: Diseased LEC3 vs LEC1

Donor-level pseudobulk DESeq2 on the IPF + SSc-ILD donors tested 2,848 genes: 263 higher in LEC3 and 114 higher in LEC1 (padj < 0.05, |log2FC| > 0.25). The genes up in LEC3 include SPARCL1 and CXCL12 (proinflammatory), POSTN and ACKR1 (venous), FLT1 (shared venous / lymphatic) and COL15A1 (extracellular-matrix stabilisation), i.e. the reprogrammed LECs adopt a proinflammatory, vein-like phenotype.

### Immune Compartment

Twenty-one immune clusters after cleaning (renumbered 0–20), annotated from their top-10 markers:

| Set | Clusters |
|---|---|
| Lymphoid | 1 T cells · 3 NK cells · 5 B cells · 15 regulatory T cells · 16 proliferating cells · 18 cytotoxic T cells |
| Monocyte / macrophage | 0 M2d macrophages · 2 macrophages · 4 monocytes · 8 M1/M2 macrophages · 9 M2a macrophages · 10 M2c macrophages · 11 complement-producing macrophages · 13 M1-like macrophages · 14 lipid-associated macrophages |
| Granulocytes / DCs | 6 neutrophils · 7 mast cells · 12 cDC2 · 17 mature DCs · 19 cDC1 · 20 plasmacytoid DCs |

IPF and SSc-ILD immune cells largely overlap and differ markedly from control. Monocytes (4), M1/M2 (8), M2c (10), M2d (0), complement-producing (11) and M1-like (13) macrophages, plasmacytoid (20) and mature (17) DCs, and T cells (1) are all expanded in diseased lungs. (The working labels written into the stage-33 script and `docs/pipeline.md` differ for a few clusters, e.g. 12, 17 and 19; the table above is the final annotation.)

### LEC–Immune Communication

With LEC0 (normal-like), LEC5 (intermediate) and LEC3 (hypoxia-reprogrammed) as signal senders and all 21 immune clusters as receivers, LEC5 and LEC3 engage in more ligand–receptor interactions than LEC0. LEC3 has the most interactions with M1-like macrophages (cluster 13) and the fewest with M2a and general macrophages, which are themselves depleted in diseased tissue. Together with the DE result this gives the working model: reduced tissue oxygenation reprograms LECs into cells with higher POSTN / ACKR1 (reduced immune-cell clearance), MHC II and SPARCL1 / CXCL12 (immune stimulation), which enhances immune activation and promotes disease progression.

### Where the Numbers Come From

The pseudobulk DE tables, GO term tables, trajectory pseudotimes and CellChat interaction tables are written as CSV / RDS files by stages 24–27 and 40 (file names in `docs/pipeline.md`); the summary statistics each script prints (numbers of significant genes per comparison, pathways per gene set, interactions per direction) are the place to read the quantitative outcome of a given run.

## Technical Implementation

### Dependencies

- **Single-cell**: Seurat ≥ 5 (v4 code paths are kept where the API differs), harmony, hdf5r (for `Read10X_h5`), Matrix
- **Statistics**: DESeq2, slingshot, SingleCellExperiment, clusterProfiler, org.Hs.eg.db, enrichplot (Bioconductor); monocle3 (optional, GitHub `cole-trapnell-lab/monocle3`)
- **Cell–cell communication**: CellChat (GitHub `sqjin/CellChat`)
- **Data handling and plotting**: dplyr, tidyr, tibble, reshape2, ggplot2, ggrepel, patchwork, gridExtra, pheatmap, RColorBrewer, viridis, scales
- **Hardware**: a single machine, no GPU; the integrated 20-sample object and all downstream steps run in memory on a laptop

### Architecture

```
single-cell/
├── config.R                           # data / results paths; sources R/; cd to results
├── R/
│   ├── colors.R                       # shared palettes
│   ├── manifest.R                     # GEO manifest, sample discovery, per-sample loader + QC
│   ├── marker_panel_heatmap.R         # EC / immune marker-panel heatmap used after each cleaning step
│   ├── pseudobulk_de.R                # donor-level pseudobulk DESeq2 between LEC subclusters + volcano
│   ├── go_enrichment.R                # enrichGO on the DE gene lists
│   └── cellchat_lec_immune.R          # LEC <-> immune CellChat + directional heatmaps
├── scripts/
│   ├── 01_load_integrate.R            # load 20 samples, QC, merge, Harmony, cluster
│   ├── 02_remove_clusters.R           # drop six contaminating joint clusters
│   ├── 03_doublet_contamination_qc.R  # doublet flag + cross-compartment module-score filter
│   ├── 04_marker_panel_heatmaps.R     # QC heatmaps for Combined / Filtered / Clean
│   ├── 10_ec_subtypes.R               # EC re-clustering + arterial/venous/capillary/lymphatic scores
│   ├── 11_ec_cluster_filter.R         # keep bona fide EC clusters
│   ├── 12_ec_subtype_heatmap.R        # subtype x group zonation heatmap
│   ├── 20_lec_subcluster.R            # LEC0-LEC7
│   ├── 21_lec_marker_plots.R          # candidate state genes
│   ├── 22_lec_signature_dotplots.R    # signature dot plots (two cluster orders)
│   ├── 23_lec_heatmap.R               # 44-gene subcluster heatmap
│   ├── 24_lec_trajectory.R            # Slingshot (+ optional Monocle3), state labels, genes vs pseudotime
│   ├── 25_lec_within_cluster_disease_de.R   # disease contrasts inside LEC1 / LEC3
│   ├── 26_lec_cluster_vs_cluster_de.R       # LEC1 vs LEC3 / LEC6 / LEC5
│   ├── 27_lec_go_enrichment.R         # GO on the stage-26 gene lists
│   ├── 30_immune_subcluster.R         # immune re-clustering (res 0.8)
│   ├── 31_immune_markers.R            # FindAllMarkers + marker heatmaps
│   ├── 32_immune_filter_renumber.R    # drop 13/20/22/24, renumber 0-20
│   ├── 33_immune_renamed_markers_heatmaps.R # lymphoid / myeloid / mixed heatmaps, M1/M2 markers
│   ├── 40_cellchat_lec_immune.R       # Complete and Targeted CellChat runs
│   ├── exploratory/                   # cleaning variants and an alternative immune start
│   └── run_all.sh                     # run everything, or just one branch
├── notebooks/original_analysis.Rmd    # the working notebook the scripts were extracted from
├── docs/
│   ├── pipeline.md                    # stage table: inputs, outputs, parameters
│   ├── changes.md                     # every edit made when turning the notebook into scripts
│   └── exploratory/research.nb.html   # rendered output of an earlier version of the analysis
└── data/                              # GSE159354_RAW goes here (git-ignored)
```

Every `.rds`, `.pdf` and `.csv` the scripts write goes to the results folder and is git-ignored; each stage re-reads the object of the stage before it, so any stage can be re-run alone once its input exists.

### Usage

```bash
# 1. put the GEO archive under data/ (or point GSE159354_RAW at it)
tar -xf GSE159354_RAW.tar -C data/GSE159354_RAW

# 2. run the whole pipeline, or one branch
bash scripts/run_all.sh            # all stages, ~hours; CellChat and FindAllMarkers dominate
bash scripts/run_all.sh ec         # stages 01-04 + endothelial / lymphatic branch
bash scripts/run_all.sh immune     # stages 01-04 + immune branch

# 3. or run / re-run a single stage from the repository root
Rscript scripts/26_lec_cluster_vs_cluster_de.R
```

```r
# use the helpers interactively
source("config.R")                               # also cd's into results/
sce_LECs <- readRDS("LECs_subclustered.rds")
de <- run_cluster_vs_cluster_de(sce_LECs, ref = "LEC1", test = "LEC3")
head(de$diseased[de$diseased$padj < 0.05, c("gene", "log2FoldChange", "padj")])
go <- run_go_for_comparison("LEC1", "LEC3")
plot_marker_panel_heatmap(readRDS("GSE159354_Clean_EC-Immune_harmony.rds"), label = "Clean")
```

### Documentation

- [docs/pipeline.md](docs/pipeline.md): every stage with the files it reads and writes and the parameters it uses, plus the immune cluster annotations
- [docs/changes.md](docs/changes.md): what was changed between the working notebook and this repository, including the bugs fixed
- [notebooks/original_analysis.Rmd](notebooks/original_analysis.Rmd): the original notebook, kept for provenance

## Research Impact

- **A clean lymphatic compartment from sorted lung data**: three cleaning passes, each verified with the same marker-panel heatmap, turn the CD31+ fraction into an endothelial object whose lymphatic subset is large enough to subcluster, which unsorted lung atlases rarely allow
- **A hypoxia-reprogrammed LEC state**: LEC3, LEC5 and LEC7 (≈ 20 % of LECs, ≥ 95 % disease-derived) sit at the end of a trajectory along which PROX1 and FLT4 fall and EPAS1 and HIF1A rise, and express venous (POSTN, ACKR1), proinflammatory (SPARCL1, CXCL12) and MHC class II programs – evidence that fibrotic LECs lose their identity and adopt a vein-like, immunostimulatory phenotype
- **LECs as direct modulators of immunity**: receptor–ligand analysis places the reprogrammed LECs in the most interactions with M1-like macrophages, linking the lymphatic vasculature to the macrophage-driven inflammation that is already known to drive fibrosis
- **Donor-level statistics**: all between-state and between-disease tests are run on donor × cluster pseudobulks with DESeq2, so the results are not inflated by thousands of cells from the same lung
- **A reusable pipeline**: the integration, cleaning, subtyping, pseudobulk DE, GO and CellChat steps are parameterised functions that can be pointed at another sorted-fraction dataset or another subcluster comparison by changing one argument

## Future Directions

- Replace the module-score argmax EC subtyping with reference mapping to a published lung endothelial atlas, and compare the subtype calls
- Test the LEC state labels with a composition model (e.g. scCODA or a donor-level proportion test) instead of raw counts
- Run the within-cluster disease contrasts (stage 25) for every LEC subcluster, not only LEC1 and LEC3, and the GO enrichment for LEC1 vs LEC5
- Repeat the CellChat analysis per disease group to ask which LEC → immune interactions are gained in IPF versus SSc-ILD
- Validate the top LEC3 / LEC5 markers and ligands on lung sections (RNAscope, immunofluorescence against PDPN / PROX1, or spatial transcriptomics), since none of the findings have been confirmed *in vivo*
- Confirm the proposed LEC reprogramming with lineage tracing in a preclinical fibrosis model, and test the roles of SPARCL1 and MHC II molecules in regulating the macrophage and T-cell response in cell culture and animal models
