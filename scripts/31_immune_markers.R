## Immune branch: cluster markers (FindAllMarkers) and marker heatmaps
## Stage 31 of the GSE159354 EC / immune pipeline.
##
## Inputs:  Immune_only_seurat_processed.rds (stage 30)
## Outputs: Immune_all_cluster_markers.csv, Immune_top10_markers_per_cluster.csv, Immune_Cluster_Marker_Heatmaps.pdf, Comprehensive_Immune_DE_Analysis.pdf, Readable_*_Heatmap.pdf, Transposed_*_Heatmap.pdf, Immune_*_summary.csv
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(pheatmap); library(gridExtra) })

immune_seurat <- readRDS("Immune_only_seurat_processed.rds")
Idents(immune_seurat) <- immune_seurat$RNA_snn_res.0.8
p_by_group <- DimPlot(immune_seurat, reduction = "umap", group.by = "group") + ggtitle("Immune cells: UMAP by disease group")
p_clusters <- DimPlot(immune_seurat, reduction = "umap", label = TRUE) + ggtitle("Immune cells: UMAP by clusters")
immune_seurat$disease_status <- factor(ifelse(immune_seurat$group == "CTRL", "CTRL", "Disease"), levels = c("CTRL", "Disease"))
immune_ctrl_disease_plot <- DimPlot(immune_seurat, group.by = "disease_status", reduction = "umap",
                                    cols = disease_status_colors) + ggtitle("Immune Cells: CTRL vs Disease") + theme_minimal()

# ── 8) DIFFERENTIAL EXPRESSION ANALYSIS ────────────────────────────────────────
cat("\n=== RUNNING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
cat("Finding markers for each immune cell cluster...\n")

# Find all markers for each cluster
# This may take several minutes depending on the number of cells
immune_markers <- FindAllMarkers(immune_seurat,
                                only.pos = TRUE,           # Only positive markers
                                min.pct = 0.25,            # Gene must be expressed in 25% of cells
                                logfc.threshold = 0.25,    # Minimum log fold change
                                verbose = TRUE)

# Check the results
cat("Total marker genes found:", nrow(immune_markers), "\n")
cat("Clusters analyzed:", length(unique(immune_markers$cluster)), "\n")

# Display top markers per cluster
top_markers_summary <- immune_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    top_gene = gene[which.max(avg_log2FC)],
    max_logfc = max(avg_log2FC),
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(cluster))

cat("\n=== TOP MARKER SUMMARY BY CLUSTER ===\n")
print(top_markers_summary)

# ── 9) SELECT TOP 10 MARKERS PER CLUSTER ───────────────────────────────────────
cat("\nSelecting top 10 upregulated genes per cluster...\n")

# Get top 10 markers per cluster based on avg_log2FC
top10_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Create a unique list of genes for the heatmap
heatmap_genes <- unique(top10_markers$gene)
cat("Total unique genes for heatmap:", length(heatmap_genes), "\n")

# Save the marker results
write.csv(immune_markers, "Immune_all_cluster_markers.csv", row.names = FALSE)
write.csv(top10_markers, "Immune_top10_markers_per_cluster.csv", row.names = FALSE)

# ── 10) CREATE HEATMAP OF TOP CLUSTER MARKERS ──────────────────────────────────
cat("\nCreating heatmap of top cluster markers...\n")

# Calculate average expression of marker genes per cluster
avg_exp <- AverageExpression(immune_seurat,
                            features = heatmap_genes,
                            group.by = "RNA_snn_res.0.8",
                            slot = "data")$RNA

# Z-score normalize across clusters (rows = genes, cols = clusters)
avg_exp_scaled <- t(scale(t(avg_exp)))

# Handle any NaN values (replace with 0)
avg_exp_scaled[!is.finite(avg_exp_scaled)] <- 0

# Create cluster annotation
cluster_names <- colnames(avg_exp_scaled)
cluster_annotation <- data.frame(
  Cluster = paste("Cluster", cluster_names),
  row.names = cluster_names
)

# Create a discrete color palette for clusters
n_clusters <- ncol(avg_exp_scaled)
cluster_colors <- rainbow(n_clusters, start = 0, end = 0.8)
names(cluster_colors) <- paste("Cluster", cluster_names)

annotation_colors <- list(Cluster = cluster_colors)

# ── 11) PLOT HEATMAP USING PHEATMAP ────────────────────────────────────────────
cat("Generating heatmap...\n")

# Create the heatmap
marker_heatmap <- pheatmap(
  avg_exp_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,          # Cluster genes
  cluster_cols = FALSE,         # Keep cluster order
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = "grey60",
  annotation_col = cluster_annotation,
  annotation_colors = annotation_colors,
  main = "Top 10 Marker Genes per Immune Cell Cluster",
  angle_col = 45
)

# ── 12) ALTERNATIVE HEATMAP: TOP 5 MARKERS ONLY ───────────────────────────────
cat("Creating simplified heatmap with top 5 markers per cluster...\n")

# Get top 5 markers per cluster for a cleaner visualization
top5_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

heatmap_genes_top5 <- unique(top5_markers$gene)

# Calculate average expression for top 5 markers
avg_exp_top5 <- AverageExpression(immune_seurat,
                                 features = heatmap_genes_top5,
                                 group.by = "RNA_snn_res.0.8",
                                 slot = "data")$RNA

# Z-score normalize
avg_exp_top5_scaled <- t(scale(t(avg_exp_top5)))
avg_exp_top5_scaled[!is.finite(avg_exp_top5_scaled)] <- 0

# Create cleaner heatmap with top 5 markers
marker_heatmap_top5 <- pheatmap(
  avg_exp_top5_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = "grey60",
  annotation_col = cluster_annotation,
  annotation_colors = annotation_colors,
  main = "Top 5 Marker Genes per Immune Cell Cluster",
  angle_col = 45
)

# ── 13) FEATURE PLOTS OF TOP MARKERS ───────────────────────────────────────────
cat("Creating feature plots of top marker genes...\n")

# Get the top marker gene for each cluster
top_marker_per_cluster <- immune_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 1) %>%
  ungroup()

# Create feature plots for top markers
top_marker_plots <- list()
for (i in 1:nrow(top_marker_per_cluster)) {
  gene <- top_marker_per_cluster$gene[i]
  cluster <- top_marker_per_cluster$cluster[i]

  top_marker_plots[[paste0("Cluster_", cluster)]] <- FeaturePlot(
    immune_seurat,
    features = gene,
    reduction = "umap"
  ) + ggtitle(paste("Cluster", cluster, "top marker:", gene))
}

# Print top marker plots
for (plot in top_marker_plots) {
  print(plot)
}

# ── 14) SUMMARY STATISTICS FOR DE ANALYSIS ────────────────────────────────────
cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")

# Markers per cluster summary
markers_per_cluster <- immune_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    avg_logfc_median = round(median(avg_log2FC), 3),
    avg_logfc_max = round(max(avg_log2FC), 3),
    top_marker = gene[which.max(avg_log2FC)],
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(cluster))

print(markers_per_cluster)

# ── 15) SAVE ALL RESULTS ───────────────────────────────────────────────────────
cat("\nSaving all results...\n")

# Save the processed immune object
saveRDS(immune_seurat, file = "Immune_only_seurat_processed.rds")

# Save heatmaps as PDF
pdf("Immune_Cluster_Marker_Heatmaps.pdf", height = 12, width = 10)
pheatmap(
  avg_exp_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = "grey60",
  annotation_col = cluster_annotation,
  annotation_colors = annotation_colors,
  main = "Top 10 Marker Genes per Immune Cell Cluster",
  angle_col = 45
)

pheatmap(
  avg_exp_top5_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = "grey60",
  annotation_col = cluster_annotation,
  annotation_colors = annotation_colors,
  main = "Top 5 Marker Genes per Immune Cell Cluster",
  angle_col = 45
)
dev.off()

# Create comprehensive analysis PDF
pdf("Comprehensive_Immune_DE_Analysis.pdf", height = 10, width = 15)

# Basic UMAP plots
print(p_by_group)
print(p_clusters)
print(immune_ctrl_disease_plot)

# Resolution comparison plots
resolutions <- c(0.3, 0.5, 0.8, 1.0, 1.2)
res_plots <- list()
for (res in resolutions) {
  cluster_col <- paste0("RNA_snn_res.", res)
  res_plots[[as.character(res)]] <- DimPlot(immune_seurat,
                                           group.by = cluster_col,
                                           label = TRUE,
                                           reduction = "umap") +
                                   ggtitle(paste("Resolution", res)) +
                                   theme_minimal() +
                                   theme(legend.position = "none")
}

library(gridExtra)
grid.arrange(grobs = res_plots, ncol = 3, nrow = 2,
            top = "Immune Cell Clustering at Different Resolutions")

# Heatmaps
pheatmap(
  avg_exp_top5_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = "grey60",
  annotation_col = cluster_annotation,
  annotation_colors = annotation_colors,
  main = "Top 5 Marker Genes per Immune Cell Cluster",
  angle_col = 45
)

# Top marker feature plots (first 9 clusters)
if (length(top_marker_plots) >= 9) {
  grid.arrange(grobs = top_marker_plots[1:9], ncol = 3, nrow = 3,
              top = "Top Marker Genes by Cluster (Feature Plots)")
}

dev.off()

# ── 16) DETAILED CLUSTER CHARACTERIZATION ──────────────────────────────────────
cat("\n=== DETAILED CLUSTER CHARACTERIZATION ===\n")

# Create detailed cluster summary with top markers
detailed_cluster_summary <- immune_markers %>%
  group_by(cluster) %>%
  summarise(
    n_cells = length(unique(WhichCells(immune_seurat, idents = cluster))),
    n_markers = n(),
    top_3_markers = paste(head(gene[order(-avg_log2FC)], 3), collapse = ", "),
    avg_logfc_range = paste0(round(min(avg_log2FC), 2), " to ", round(max(avg_log2FC), 2)),
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(cluster))

print(detailed_cluster_summary)

# ── 17) CLUSTER COMPOSITION ANALYSIS ───────────────────────────────────────────
cat("\n=== CLUSTER COMPOSITION BY DISEASE GROUP ===\n")

cluster_composition <- immune_seurat@meta.data %>%
  group_by(RNA_snn_res.0.8, group) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  group_by(RNA_snn_res.0.8) %>%
  mutate(
    total_in_cluster = sum(n_cells),
    percentage = round((n_cells / total_in_cluster) * 100, 1)
  ) %>%
  select(cluster = RNA_snn_res.0.8, group, n_cells, percentage) %>%
  arrange(as.numeric(cluster))

print(cluster_composition)

# Create composition visualization
composition_plot <- immune_seurat@meta.data %>%
  group_by(RNA_snn_res.0.8, group) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  ggplot(aes(x = RNA_snn_res.0.8, y = n_cells, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CTRL" = "#2E8B57", "IPF" = "#9370DB", "SSc-ILD" = "#DC143C")) +
  labs(title = "Immune Cluster Composition by Disease Group",
       x = "Cluster", y = "Number of Cells", fill = "Disease Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(composition_plot)

# ── 18) SAVE ALL SUMMARY TABLES ────────────────────────────────────────────────
write.csv(detailed_cluster_summary, "Immune_detailed_cluster_summary.csv", row.names = FALSE)
write.csv(cluster_composition, "Immune_cluster_composition_analysis.csv", row.names = FALSE)

# ── 19) FINAL SUMMARY ──────────────────────────────────────────────────────────
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("=== IMMUNE CELL DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

cat("\nDATASET OVERVIEW:\n")
cat("  Total immune cells:", ncol(immune_seurat), "\n")
cat("  Number of clusters (res 0.8):", length(unique(immune_seurat$RNA_snn_res.0.8)), "\n")
cat("  Total marker genes identified:", nrow(immune_markers), "\n")
cat("  Unique genes in heatmap:", length(heatmap_genes), "\n")

cat("\nFILES GENERATED:\n")
cat("  1. Immune_only_seurat_processed.rds - Main Seurat object\n")
cat("  2. Immune_all_cluster_markers.csv - All cluster markers\n")
cat("  3. Immune_top10_markers_per_cluster.csv - Top 10 markers per cluster\n")
cat("  4. Immune_detailed_cluster_summary.csv - Cluster characteristics\n")
cat("  5. Immune_cluster_composition_analysis.csv - Disease group composition\n")
cat("  6. Immune_Cluster_Marker_Heatmaps.pdf - Heatmap visualizations\n")
cat("  7. Comprehensive_Immune_DE_Analysis.pdf - Complete analysis report\n")

cat("\nHEATMAP SUMMARY:\n")
cat("  Two heatmaps created:\n")
cat("    - Top 10 markers per cluster (", length(heatmap_genes), " total genes)\n")
cat("    - Top 5 markers per cluster (", length(unique(top5_markers$gene)), " total genes)\n")
cat("  Color scheme: Navy (low) -> White (medium) -> Red (high expression)\n")

cat("\nNEXT STEPS:\n")
cat("  1. Review the heatmaps to identify potential immune cell types\n")
cat("  2. Use known immune markers to annotate clusters:\n")
cat("     - T cells: CD3D+, CD4+ or CD8A+\n")
cat("     - B cells: MS4A1+, CD79A+\n")
cat("     - NK cells: NKG7+, GNLY+\n")
cat("     - Macrophages: CD68+, CD163+\n")
cat("     - Monocytes: CD14+, FCGR3A+\n")
cat("  3. Validate annotations with known immune biology\n")
cat("  4. Compare cluster distributions across disease groups\n")

cat("\n=== ANALYSIS READY FOR IMMUNE CELL TYPE ANNOTATION ===\n")


# ===== CREATE READABLE HEATMAP WITH BETTER SPACING =====
library(pheatmap)

# Get top 5 markers per cluster (assuming immune_markers already exists)
top5_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

heatmap_genes <- unique(top5_markers$gene)

# Calculate average expression per cluster
avg_exp <- AverageExpression(immune_seurat,
                            features = heatmap_genes,
                            group.by = "RNA_snn_res.0.8",
                            slot = "data")$RNA

# Z-score normalize
avg_exp_scaled <- t(scale(t(avg_exp)))
avg_exp_scaled[!is.finite(avg_exp_scaled)] <- 0

# Create cluster annotation
cluster_annotation <- data.frame(
  Cluster = paste("Cluster", colnames(avg_exp_scaled)),
  row.names = colnames(avg_exp_scaled)
)

# ===== SAVE AS READABLE PDF =====
# Make PDF much taller and adjust font sizes for readability
pdf("Readable_Immune_Markers_Heatmap.pdf", height = 30, width = 12)

pheatmap(
  avg_exp_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 12,        # Larger font for gene names
  fontsize_col = 14,        # Larger font for cluster names
  border_color = "grey60",
  annotation_col = cluster_annotation,
  main = "Top 5 Marker Genes per Immune Cell Cluster",
  angle_col = 0,           # Horizontal cluster labels
  cellwidth = 25,          # Wider cells
  cellheight = 15          # Taller cells for better spacing
)

dev.off()


# ===== CREATE READABLE HEATMAP WITH TOP 10 MARKERS =====
library(pheatmap)

# Get top 10 markers per cluster (assuming immune_markers already exists)
top10_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

heatmap_genes_top10 <- unique(top10_markers$gene)

# Calculate average expression per cluster
avg_exp_top10 <- AverageExpression(immune_seurat,
                                  features = heatmap_genes_top10,
                                  group.by = "RNA_snn_res.0.8",
                                  slot = "data")$RNA

# Z-score normalize
avg_exp_top10_scaled <- t(scale(t(avg_exp_top10)))
avg_exp_top10_scaled[!is.finite(avg_exp_top10_scaled)] <- 0

# Create cluster annotation
cluster_annotation <- data.frame(
  Cluster = paste("Cluster", colnames(avg_exp_top10_scaled)),
  row.names = colnames(avg_exp_top10_scaled)
)

# ===== SAVE AS READABLE PDF =====
# Make PDF even taller for more genes and adjust spacing
pdf("Readable_Top10_Immune_Markers_Heatmap.pdf", height = 55, width = 12)

pheatmap(
  avg_exp_top10_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 12,        # Larger font for gene names
  fontsize_col = 14,        # Larger font for cluster names
  border_color = "grey60",
  annotation_col = cluster_annotation,
  main = "Top 10 Marker Genes per Immune Cell Cluster",
  angle_col = 0,           # Horizontal cluster labels
  cellwidth = 25,          # Wider cells
  cellheight = 15          # Taller cells for better spacing
)

dev.off()

cat("Readable top 10 heatmap saved as: Readable_Top10_Immune_Markers_Heatmap.pdf\n")
cat("Dimensions: 30 inches tall x 12 inches wide for easy reading\n")


# ===== TRANSPOSED HEATMAPS: GENES ON X-AXIS, CLUSTERS ON Y-AXIS =====
library(pheatmap)

# ===== TOP 5 MARKERS TRANSPOSED HEATMAP =====
# Get top 5 markers per cluster
top5_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

heatmap_genes_top5 <- unique(top5_markers$gene)

# Calculate average expression per cluster
avg_exp_top5 <- AverageExpression(immune_seurat,
                                 features = heatmap_genes_top5,
                                 group.by = "RNA_snn_res.0.8",
                                 slot = "data")$RNA

# Z-score normalize
avg_exp_top5_scaled <- t(scale(t(avg_exp_top5)))
avg_exp_top5_scaled[!is.finite(avg_exp_top5_scaled)] <- 0

# TRANSPOSE: genes become columns (x-axis), clusters become rows (y-axis)
avg_exp_top5_transposed <- t(avg_exp_top5_scaled)

# Create row annotation for clusters
cluster_annotation_rows <- data.frame(
  Cluster = paste("Cluster", rownames(avg_exp_top5_transposed)),
  row.names = rownames(avg_exp_top5_transposed)
)

# Save top 5 transposed heatmap
pdf("Transposed_Top5_Immune_Markers_Heatmap.pdf", height = 12, width = 40)

pheatmap(
  avg_exp_top5_transposed,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = FALSE,      # Keep cluster order
  cluster_cols = TRUE,       # Cluster genes
  show_colnames = TRUE,      # Show gene names (now on x-axis)
  show_rownames = TRUE,      # Show cluster names (now on y-axis)
  fontsize_row = 14,         # Cluster names font
  fontsize_col = 12,         # Gene names font
  border_color = "grey60",
  annotation_row = cluster_annotation_rows,
  main = "Top 5 Marker Genes per Immune Cell Cluster (Transposed)",
  angle_col = 45,           # Angle gene names for readability
  cellwidth = 20,          # Width for each gene
  cellheight = 25          # Height for each cluster
)

dev.off()

# ===== TOP 10 MARKERS TRANSPOSED HEATMAP =====
# Get top 10 markers per cluster
top10_markers <- immune_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

heatmap_genes_top10 <- unique(top10_markers$gene)

# Calculate average expression per cluster
avg_exp_top10 <- AverageExpression(immune_seurat,
                                  features = heatmap_genes_top10,
                                  group.by = "RNA_snn_res.0.8",
                                  slot = "data")$RNA

# Z-score normalize
avg_exp_top10_scaled <- t(scale(t(avg_exp_top10)))
avg_exp_top10_scaled[!is.finite(avg_exp_top10_scaled)] <- 0

# TRANSPOSE: genes become columns (x-axis), clusters become rows (y-axis)
avg_exp_top10_transposed <- t(avg_exp_top10_scaled)

# Create row annotation for clusters
cluster_annotation_rows_top10 <- data.frame(
  Cluster = paste("Cluster", rownames(avg_exp_top10_transposed)),
  row.names = rownames(avg_exp_top10_transposed)
)

# Save top 10 transposed heatmap (extra wide for more genes)
pdf("Transposed_Top10_Immune_Markers_Heatmap.pdf", height = 12, width = 60)

pheatmap(
  avg_exp_top10_transposed,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = FALSE,      # Keep cluster order
  cluster_cols = TRUE,       # Cluster genes
  show_colnames = TRUE,      # Show gene names (now on x-axis)
  show_rownames = TRUE,      # Show cluster names (now on y-axis)
  fontsize_row = 14,         # Cluster names font
  fontsize_col = 12,         # Gene names font
  border_color = "grey60",
  annotation_row = cluster_annotation_rows_top10,
  main = "Top 10 Marker Genes per Immune Cell Cluster (Transposed)",
  angle_col = 45,           # Angle gene names for readability
  cellwidth = 18,          # Width for each gene
  cellheight = 25          # Height for each cluster
)

dev.off()

cat("Transposed heatmaps created:\n")
cat("- Top 5: Transposed_Top5_Immune_Markers_Heatmap.pdf (12h x 25w inches)\n")
cat("- Top 10: Transposed_Top10_Immune_Markers_Heatmap.pdf (12h x 40w inches)\n")
cat("Gene names are now horizontal and much easier to read!\n")
