## Load the CD31+ (EC) and CD45+ (immune) fractions of GSE159354, QC, merge, integrate and cluster
## Stage 01 of the GSE159354 EC / immune pipeline.
##
## Inputs:  <GSE159354_RAW>/**/GSM*_filtered_feature_bc_matrix.h5 (20 files: 10 donors x 2 fractions)
## Outputs: GSE159354_ThreeGroup_EC-Immune_harmony.rds, ThreeGroup_Markers.pdf,
##          ThreeGroup_Summary_Plots.pdf, ThreeGroup_sample_summary_EC-Immune.csv,
##          ThreeGroup_cell_counts_by_group_celltype.csv
## Per-sample QC: genes in >= 3 cells, 200 < nFeature_RNA < 8000, percent.mt < 20. Merged
## object: LogNormalize (1e4), 3000 variable genes (vst), 50 PCs, Harmony on donor + sorted
## fraction, UMAP / SNN on Harmony dims 1-30, Louvain clustering at resolution 0.6.
##
## Run from the repository root:  Rscript scripts/01_load_integrate.R

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(tibble)
  library(tidyr)
})

## ===== Sample table: CD31+ and CD45+ fractions of all three groups =====
filtered_manifest <- build_sample_table(cfg$data_dir, fractions = c("CD31p", "CD45p"))
cat("Samples found:", nrow(filtered_manifest), "\n")
print(filtered_manifest %>% count(group, cellset))

## ===== Load, QC and merge =====
sample_objects <- lapply(seq_len(nrow(filtered_manifest)), function(i) {
  row <- filtered_manifest[i, ]
  cat("Loading:", row$sample_id, "\n")
  obj <- load_sample(row)
  cat("  Kept", ncol(obj), "cells after QC\n")
  obj
})
names(sample_objects) <- filtered_manifest$sample_id
sample_objects <- Filter(function(o) ncol(o) > 0, sample_objects)
stopifnot(length(sample_objects) > 0)

combined_seurat <- Reduce(function(x, y) merge(x, y), sample_objects)
cat("Total cells after merge:", ncol(combined_seurat), "\n")

## ===== Standard workflow =====
DefaultAssay(combined_seurat) <- "RNA"
combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
combined_seurat <- FindVariableFeatures(combined_seurat, nfeatures = 3000, selection.method = "vst", verbose = FALSE)
combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(combined_seurat, npcs = 50, verbose = FALSE)

## Harmony over donor and sorted fraction
combined_seurat <- RunHarmony(combined_seurat, group.by.vars = c("donor", "fraction"), verbose = TRUE)

## UMAP / neighbors / clusters
combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:30)
combined_seurat <- FindNeighbors(combined_seurat, reduction = "harmony", dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.6)

## ===== Save the integrated object =====
saveRDS(combined_seurat, file = "GSE159354_ThreeGroup_EC-Immune_harmony.rds")

## ===== UMAP plots =====
celltype_plot <- DimPlot(combined_seurat, group.by = "cellset", reduction = "umap", cols = cellset_colors) +
  ggtitle("UMAP: Endothelial Cells vs Immune Cells") +
  theme(legend.title = element_text(size = 12))
group_plot <- DimPlot(combined_seurat, group.by = "group", reduction = "umap", cols = group_colors) +
  ggtitle("UMAP: CTRL vs IPF vs SSc-ILD")
cluster_plot <- DimPlot(combined_seurat, label = TRUE, reduction = "umap") +
  ggtitle("UMAP: Clusters")
print(celltype_plot); print(group_plot); print(cluster_plot)

## ===== Summary statistics =====
cell_count_summary <- combined_seurat@meta.data %>%
  group_by(group, cellset) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = cellset, values_from = n_cells, values_fill = 0) %>%
  mutate(Total = EC + Immune)
group_totals    <- combined_seurat@meta.data %>% count(group,   name = "n_cells")
celltype_totals <- combined_seurat@meta.data %>% count(cellset, name = "n_cells")

cat("\n=== CELL COUNT SUMMARY ===\n");        print(cell_count_summary)
cat("\n=== TOTAL CELLS BY GROUP ===\n");      print(group_totals)
cat("\n=== TOTAL CELLS BY CELL TYPE ===\n");  print(celltype_totals)
cat("\n=== OVERALL STATISTICS ===\n")
cat("Total cells:", ncol(combined_seurat), "\n")
cat("Total genes:", nrow(combined_seurat), "\n")
cat("Number of samples:", length(sample_objects), "\n")
cat("Number of donors:", length(unique(combined_seurat$donor)), "\n")
cat("Number of clusters:", length(levels(combined_seurat$seurat_clusters)), "\n")

## ===== Marker plots =====
pdf("ThreeGroup_Markers.pdf", height = 10, width = 20)
print(FeaturePlot(combined_seurat,
                  features = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "PTPRC", "MS4A1", "LYZ"),
                  ncol = 4, reduction = "umap"))
dev.off()

## ===== Per-sample summary tables =====
detailed_summary <- data.frame(
  sample_id = names(sample_objects),
  donor     = sapply(sample_objects, \(o) unique(o$donor)),
  group     = sapply(sample_objects, \(o) as.character(unique(o$group))),
  cellset   = sapply(sample_objects, \(o) unique(o$cellset)),
  n_cells   = sapply(sample_objects, ncol),
  n_genes   = sapply(sample_objects, nrow),
  row.names = NULL
)
write.csv(detailed_summary,   "ThreeGroup_sample_summary_EC-Immune.csv",       row.names = FALSE)
write.csv(cell_count_summary, "ThreeGroup_cell_counts_by_group_celltype.csv", row.names = FALSE)

## ===== Composition plots =====
pdf("ThreeGroup_Summary_Plots.pdf", height = 8, width = 12)
count_plot <- combined_seurat@meta.data %>%
  count(group, cellset, name = "n_cells") %>%
  ggplot(aes(x = group, y = n_cells, fill = cellset)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cellset_colors) +
  labs(title = "Cell Counts by Disease Group and Cell Type", x = "Disease Group", y = "Number of Cells", fill = "Cell Type") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(count_plot)
prop_plot <- combined_seurat@meta.data %>%
  count(group, cellset, name = "n_cells") %>%
  group_by(group) %>% mutate(proportion = n_cells / sum(n_cells)) %>%
  ggplot(aes(x = group, y = proportion, fill = cellset)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cellset_colors) +
  labs(title = "Cell Type Proportions by Disease Group", x = "Disease Group", y = "Proportion", fill = "Cell Type") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(prop_plot)
dev.off()

cat("\nSaved GSE159354_ThreeGroup_EC-Immune_harmony.rds\n")
