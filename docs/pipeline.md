# Pipeline stages

All file names are relative to the results folder (`SC_RESULTS`, default `results/`).
Every stage starts by `readRDS()`-ing the object of the stage it depends on, so a stage
can be re-run on its own once its input exists.

```
01 load + integrate ─► 02 remove clusters ─► 03 doublet / contamination ─► 04 marker heatmaps
        │                                              │
        │ (ThreeGroup.rds)                             │ (Clean.rds)
        ▼                                              ▼
 10 EC subtypes ─► 11 EC cluster filter ─► 12      30 immune subcluster ─► 31 markers ─► 32 filter + renumber ─► 33
        │                                                                                      │
        ▼                                                                                      │
 20 LEC subcluster ─► 21, 22, 23 (plots) ─► 24 trajectory ─► 25, 26 DE ─► 27 GO                │
        │                                                                                      │
        └──────────────────────────────────► 40 CellChat LEC ◄─► immune ◄────────────────────┘
```

## Core

| Stage | Script | Reads | Writes | Parameters |
|---|---|---|---|---|
| 01 | `01_load_integrate.R` | 20 `*_filtered_feature_bc_matrix.h5` (CD31+, CD45+ of 10 donors) | `GSE159354_ThreeGroup_EC-Immune_harmony.rds`, `ThreeGroup_*.pdf/csv` | min.cells 3; 200 < nFeature < 8000; mt < 20%; LogNormalize 1e4; 3000 HVG; 50 PCs; Harmony(donor, fraction); UMAP/SNN dims 1-30; res 0.6 |
| 02 | `02_remove_clusters.R` | ThreeGroup | `GSE159354_FilteredClusters_EC-Immune_harmony.rds`, `Filtered_UMAP_Analysis.pdf` | drops joint clusters 7, 12, 15, 16, 17, 19 |
| 03 | `03_doublet_contamination_qc.R` | FilteredClusters | `GSE159354_Clean_EC-Immune_harmony.rds`, `Cell_Quality_Diagnostics_v5.pdf` | doublets = top 5% of z(nFeature)+z(nCount); EC with immune score > 90th pct of ECs, immune with EC score > 90th pct of immune cells removed |
| 04 | `04_marker_panel_heatmaps.R` | ThreeGroup, FilteredClusters, Clean | `<Label>_Expression_Heatmap.pdf`, `Enhanced_<Label>_Expression_Heatmap.pdf` | 37-gene LEC / arterial / venous / core-EC / immune panel, row z-score clipped at ±2.5 |

Exploratory variants of stage 03 (`scripts/exploratory/03b_*`, `03c_*`) produce
`GSE159354_UltraClean_*` and `GSE159354_ClusterFiltered_*`; nothing downstream reads them.

## Endothelial / lymphatic branch (starts from ThreeGroup)

| Stage | Script | Reads | Writes | Parameters |
|---|---|---|---|---|
| 10 | `10_ec_subtypes.R` | ThreeGroup | `EC_subtypes_assigned_all_groups.rds`, `1-6_EC_Subtypes_*.pdf`, `EC_subtype_*.csv` | ECs re-processed: 2000 HVG, 50 PCs, Harmony(donor, group), res 0.8; subtype = argmax of AddModuleScore over arterial (GJA5, EFNB2, DLL4, SOX17, HEY2, NOTCH4, BMX, SEMA3G), venous (NR2F2, EPHB4, ACKR1, VCAM1, SELP, SELE, VWF, PLVAP), capillary (CA4, EDNRB, EMCN, ESAM, KDR, RGCC, AQP1, HSPG2) and lymphatic (PROX1, LYVE1, PDPN, FLT4, VEGFC, CCL21, RELN, TBX1, FOXC2, NFATC1) signatures |
| 11 | `11_ec_cluster_filter.R` | EC_subtypes_assigned | `sce_EC_filtered_clusters.rds`, `*_Filtered_UMAP.pdf`, `EC_filtering_comparison_summary.csv` | keeps EC clusters 0, 1, 2, 3, 6, 7, 10, 20 |
| 12 | `12_ec_subtype_heatmap.R` | sce_EC_filtered_clusters | `EC_Filtered_Pheatmap_Heatmap*.pdf`, `EC_filtered_heatmap_*.csv` | 30 EC-zonation genes averaged per subtype × group, row-scaled |
| 20 | `20_lec_subcluster.R` | sce_EC_filtered_clusters | `LECs_subclustered.rds`, `LEC_*_UMAP.pdf`, `LEC_Composition_Bars.pdf`, `LEC_subcluster_*.csv` | LECs re-processed: 2000 HVG, Harmony(donor, group), res 0.7 → LEC0–LEC7 |
| 21 | `21_lec_marker_plots.R` | LECs_subclustered | `LEC_Candidate_Gene_Plots.pdf` | violin (LEC1–7) + feature plots of SPARCL1, IL3RA, POSTN, KLF4, TBX1, EBF1, COL15A1, STC1, HIF1A, EDNRB, FLT4 |
| 22 | `22_lec_signature_dotplots.R` | LECs_subclustered | `LEC_Signature_DotPlot_*.pdf`, `LEC_*_cluster_*.csv` | scaled mean expression × % expressing; two cluster orders (LEC0–7; LEC0, 6, 1, 4, 2, 3/5/7) |
| 23 | `23_lec_heatmap.R` | LECs_subclustered | `LEC_Comprehensive_Heatmap_*.pdf`, `LEC_Core_Genes_Heatmap_*.pdf`, `LEC_heatmap_*.csv` | 45 genes in 8 categories, per-subcluster means, row z-score |
| 24 | `24_lec_trajectory.R` | LECs_subclustered | `LEC_Slingshot_Trajectories.rds`, `LEC_Slingshot_Pseudotime.csv`, `LEC_*Trajector*.pdf`, `LEC_Pseudotime_*.pdf`, `LEC_Gene_*.pdf`, `LEC_Trajectory_Gene_Analysis.csv` (+ `LEC_Monocle3_*` if monocle3 is installed) | Slingshot on UMAP, auto and guided (start LEC0, ends LEC6/LEC7, 150 approx points); LEC state labels Normal (LEC1, LEC6) / Proinflammatory (LEC0, 2, 3, 5) / Intermediate (LEC4, 7); PROX1, EPAS1, HIF1A, SPARCL1 along pseudotime |
| 25 | `25_lec_within_cluster_disease_de.R` | LECs_subclustered | `LEC1_DE_*.csv`, `LEC3_DE_*.csv`, `LEC1_vs_LEC3_*_Summary.csv`, `*_Volcano.pdf`, `LEC1_vs_LEC3_Top_Genes_Heatmap.pdf` | pseudobulk per donor × group × subcluster (≥ 5 cells); DESeq2 `~ group`; IPF vs CTRL, SSc-ILD vs CTRL, IPF vs SSc-ILD inside LEC1 and inside LEC3 |
| 26 | `26_lec_cluster_vs_cluster_de.R` | LECs_subclustered | `ALL_/DISEASED_<ref>_vs_<test>_DE_results.csv`, `*_Volcano.pdf`, `TOP_50_*.csv`, `*_upregulated_genes.csv`, `<ref>_vs_<test>_Volcano_Summary_Stats.csv`, `LEC_cluster_vs_cluster_DE_results.rds` | pseudobulk per donor × subcluster (≥ 5 cells; genes ≥ 10 counts in ≥ 3 samples); DESeq2 `~ cluster`; LEC1 vs LEC3, LEC1 vs LEC6, LEC1 vs LEC5; each on all donors and on IPF + SSc-ILD donors |
| 27 | `27_lec_go_enrichment.R` | stage-26 CSVs | `<set>_<BP\|MF\|CC>_GO_results.csv`, `*_dotplot.pdf`, `*_barplot.pdf`, `GO_Pathway_Comparison_*.pdf`, `GO_*_Summary_*.csv`, `LEC_GO_enrichment_results.rds` | clusterProfiler enrichGO, org.Hs.eg.db, universe = tested genes, padj < 0.05 & \|log2FC\| > 0.25 gene lists, BH < 0.05, q < 0.2; LEC1 vs LEC3 and LEC1 vs LEC6 |

## Immune branch (starts from Clean)

| Stage | Script | Reads | Writes | Parameters |
|---|---|---|---|---|
| 30 | `30_immune_subcluster.R` | Clean | `Immune_only_seurat_processed.rds`, `Immune_UMAP_*.pdf` | immune cells re-processed: 3000 HVG, 50 PCs, Harmony over donor / group / fraction / gsm (those with ≥ 2 levels), res 0.3–1.2 computed, 0.8 used |
| 31 | `31_immune_markers.R` | Immune_only_seurat_processed | `Immune_all_cluster_markers.csv`, `Immune_top10_markers_per_cluster.csv`, `Immune_Cluster_Marker_Heatmaps.pdf`, `Comprehensive_Immune_DE_Analysis.pdf`, `Readable_*`/`Transposed_*` heatmaps | FindAllMarkers only.pos, min.pct 0.25, logfc 0.25; top-5 / top-10 marker heatmaps |
| 32 | `32_immune_filter_renumber.R` | Immune_only_seurat_processed | `Immune_Renamed_0to20_Final.rds`, `Final_*.pdf`, `Clean_Immune_Cluster_Disease_Composition.pdf/.png` | removes res-0.8 clusters 13 (epithelial contamination), 20 (unknown), 22 (proliferating B), 24 (mast/basophil); remaining 21 renumbered 0–20 (`renamed_clusters`) |
| 33 | `33_immune_renamed_markers_heatmaps.R` | Immune_Renamed_0to20_Final | `renamed_markers.rds`, `Renamed_Clusters_*_Markers.csv`, `Selected_Clusters_*.pdf`, `Myeloid_Clusters_*.pdf`, `Mixed_Immune_*.pdf`, `M1_M2_Markers_FeaturePlot.pdf` | lymphoid set 1, 3, 5, 15, 16, 18; myeloid set 0, 2, 4, 8, 9, 10, 11, 13, 14; mixed set 6, 7, 12, 16, 17, 19, 20 |

`scripts/exploratory/30b_*` re-clusters the immune cells directly from ThreeGroup
(`immune_cells_clustered.rds`, IMM0… labels); nothing downstream reads it.

## Joint

| Stage | Script | Reads | Writes | Parameters |
|---|---|---|---|---|
| 40 | `40_cellchat_lec_immune.R` | LECs_subclustered, Immune_Renamed_0to20_Final | `Complete_*`, `Targeted_*` (CellChat .rds, network / pathway PDFs, communication CSVs), `Targeted_LEC_Signaling_Heatmaps.pdf`, `Targeted_LEC_Signaling_Matrix.pdf` | raw counts on shared genes, CellChatDB.human, raw.use = TRUE, min.cells 10; "Complete" = LEC0–7, "Targeted" = LEC0 / LEC3 / LEC5 as senders |

## Immune cluster annotations (renumbered ids, from stage 33)

| Set | Cluster → annotation |
|---|---|
| Lymphoid | 1 T cells · 3 NK cells · 5 B cells · 15 regulatory T cells · 16 proliferating cells · 18 cytotoxic T cells |
| Myeloid | 0, 8, 9, 10, 11 M2-like macrophages · 2 macrophages · 4 monocytes · 13 interferon-responsive myeloid · 14 lipid-associated macrophages |
| Mixed | 6 S100A8+ neutrophils · 7 TPSAB1+ mast cells · 12 JCHAIN+ plasma cells · 16 MKI67+ proliferating · 17 CLEC10A+ dendritic cells · 19 PRF1+ tissue-resident CTL · 20 LILRA4+ plasmacytoid DC |
