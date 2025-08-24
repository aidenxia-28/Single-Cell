## Cluster-based removal of EC cells from immune-dominated clusters (exploratory)
## Exploratory variant of stage 03; not used by the main branches.
##
## Inputs:  GSE159354_Clean_EC-Immune_harmony.rds (stage 03)
## Outputs: GSE159354_ClusterFiltered_EC-Immune_harmony.rds, Cluster_Filtering_Diagnostics.pdf, ClusterFiltered_Expression_Heatmap.pdf
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

clean_lung_seurat <- readRDS("GSE159354_Clean_EC-Immune_harmony.rds")

# ===== CLUSTER-BASED FILTERING APPROACH =====

# Target clusters to remove EC cells from
target_clusters <- c("1", "2", "5", "9", "13", "14")

cat("=== CLUSTER-BASED FILTERING ===\n")
cat("Target clusters for EC removal:", paste(target_clusters, collapse = ", "), "\n")

# Check current cluster distribution
cat("\n--- Current cluster distribution ---\n")
cluster_table <- table(clean_lung_seurat$seurat_clusters, clean_lung_seurat$cellset)
print(cluster_table)

# Get EC cells in target clusters
ec_cells_all <- WhichCells(clean_lung_seurat, expression = cellset == "EC")
ec_cells_in_target_clusters <- WhichCells(clean_lung_seurat,
                                         expression = cellset == "EC" &
                                         seurat_clusters %in% target_clusters)

cat("\n--- EC cells analysis ---\n")
cat("Total EC cells:", length(ec_cells_all), "\n")
cat("EC cells in target clusters:", length(ec_cells_in_target_clusters), "\n")
cat("Percentage of EC cells to remove:",
    round((length(ec_cells_in_target_clusters)/length(ec_cells_all))*100, 1), "%\n")

# Break down by individual cluster
cat("\nEC cells per target cluster:\n")
for (cluster in target_clusters) {
  ec_in_cluster <- sum(clean_lung_seurat$seurat_clusters == cluster &
                       clean_lung_seurat$cellset == "EC")
  total_in_cluster <- sum(clean_lung_seurat$seurat_clusters == cluster)
  cat("  Cluster", cluster, ":", ec_in_cluster, "EC cells out of", total_in_cluster, "total\n")
}

# ===== CHECK WHAT'S IN THESE CLUSTERS =====
# Let's see what cell types are actually in these target clusters
cat("\n--- Target cluster composition ---\n")
for (cluster in target_clusters) {
  cluster_cells <- WhichCells(clean_lung_seurat, expression = seurat_clusters == cluster)
  cluster_celltype <- table(clean_lung_seurat$cellset[cluster_cells])
  cat("Cluster", cluster, "composition:", paste(names(cluster_celltype), cluster_celltype, collapse = ", "), "\n")
}

# ===== VISUALIZE TARGET CLUSTERS =====
# Mark the target clusters on UMAP
clean_lung_seurat$cluster_filter_status <- "Keep"
clean_lung_seurat$cluster_filter_status[clean_lung_seurat$seurat_clusters %in% target_clusters] <- "Target_Cluster"

cluster_viz_plot <- DimPlot(clean_lung_seurat, group.by = "cluster_filter_status",
                           reduction = "umap",
                           cols = c("Keep" = "lightgray", "Target_Cluster" = "orange")) +
                           ggtitle("Target Clusters for EC Removal\n(Orange = Target Clusters)")

print(cluster_viz_plot)

# Show original cluster plot for reference
cluster_plot <- DimPlot(clean_lung_seurat, group.by = "seurat_clusters",
                       reduction = "umap", label = TRUE) +
                       ggtitle("Original Clusters (Reference)")

print(cluster_plot)

# ===== CREATE CLUSTER-FILTERED OBJECT =====
# Remove EC cells from target clusters, keep all immune cells
cells_to_remove <- ec_cells_in_target_clusters
cells_to_keep <- !colnames(clean_lung_seurat) %in% cells_to_remove

cluster_filtered_seurat <- subset(clean_lung_seurat,
                                 cells = colnames(clean_lung_seurat)[cells_to_keep])

cat("\n=== CLUSTER FILTERING RESULTS ===\n")
cat("Original clean object:", ncol(clean_lung_seurat), "cells\n")
cat("Cluster filtered object:", ncol(cluster_filtered_seurat), "cells\n")
cat("EC cells removed:", length(cells_to_remove), "\n")
cat("Retention rate:", round((ncol(cluster_filtered_seurat)/ncol(clean_lung_seurat))*100, 1), "%\n")

# Cell type breakdown
original_ec <- sum(clean_lung_seurat$cellset == "EC")
original_immune <- sum(clean_lung_seurat$cellset == "Immune")
filtered_ec <- sum(cluster_filtered_seurat$cellset == "EC")
filtered_immune <- sum(cluster_filtered_seurat$cellset == "Immune")

cat("\nCell type changes:\n")
cat("EC cells:", original_ec, "→", filtered_ec,
    "(removed:", original_ec - filtered_ec, ",",
    round(((original_ec - filtered_ec)/original_ec)*100, 1), "%)\n")
cat("Immune cells:", original_immune, "→", filtered_immune, "(no change)\n")

# ===== CHECK REMAINING CLUSTER DISTRIBUTION =====
cat("\n--- Remaining cluster distribution ---\n")
remaining_cluster_table <- table(cluster_filtered_seurat$seurat_clusters,
                               cluster_filtered_seurat$cellset)
print(remaining_cluster_table)

# ===== VISUALIZE RESULTS =====
# Mark what was removed
clean_lung_seurat$ec_removal_status <- ifelse(clean_lung_seurat$cellset == "EC" &
                                              clean_lung_seurat$seurat_clusters %in% target_clusters,
                                              "EC_Removed",
                                              ifelse(clean_lung_seurat$cellset == "EC", "EC_Kept", "Immune"))

removal_plot <- DimPlot(clean_lung_seurat, group.by = "ec_removal_status",
                       reduction = "umap",
                       cols = c("EC_Kept" = "#E31A1C", "Immune" = "#1F78B4", "EC_Removed" = "black")) +
                       ggtitle("Cluster-Based EC Removal\n(Black = Removed EC cells)")

print(removal_plot)

# Final clean result
final_cluster_plot <- DimPlot(cluster_filtered_seurat, group.by = "cellset",
                            reduction = "umap",
                            cols = c("EC" = "#E31A1C", "Immune" = "#1F78B4")) +
                            ggtitle("Cluster-Filtered Object: EC vs Immune")

print(final_cluster_plot)

# ===== TEST FOR IMMUNE CONTAMINATION IN REMAINING EC CELLS =====
cat("\n=== CHECKING CONTAMINATION IN REMAINING EC CELLS ===\n")
remaining_ec <- WhichCells(cluster_filtered_seurat, expression = cellset == "EC")
problematic_genes <- c("MRC1", "CD163", "CD68", "LYZ")
available_genes <- problematic_genes[problematic_genes %in% rownames(cluster_filtered_seurat)]

for (gene in available_genes) {
  gene_expr <- FetchData(cluster_filtered_seurat, vars = gene, cells = remaining_ec, layer = "data")
  above_1 <- sum(gene_expr[[gene]] > 1.0)
  above_2 <- sum(gene_expr[[gene]] > 2.0)
  max_expr <- max(gene_expr[[gene]])

  cat("Gene", gene, "in remaining EC cells:\n")
  cat("  Expression > 1.0:", above_1, "/", length(remaining_ec),
      "(", round((above_1/length(remaining_ec))*100, 1), "%)\n")
  cat("  Expression > 2.0:", above_2, "/", length(remaining_ec),
      "(", round((above_2/length(remaining_ec))*100, 1), "%)\n")
  cat("  Max expression:", round(max_expr, 2), "\n\n")
}

# Save the cluster-filtered object
saveRDS(cluster_filtered_seurat, file = "GSE159354_ClusterFiltered_EC-Immune_harmony.rds")

# Save diagnostic plots
pdf("Cluster_Filtering_Diagnostics.pdf", height = 12, width = 16)
print(cluster_viz_plot)
print(cluster_plot)
print(removal_plot)
print(final_cluster_plot)
dev.off()

cat("=== CLUSTER-BASED FILTERING COMPLETE ===\n")
cat("New object: 'cluster_filtered_seurat'\n")
cat("Saved as: GSE159354_ClusterFiltered_EC-Immune_harmony.rds\n")
cat("This approach uses clustering to identify problematic EC populations\n")
cat("Run your heatmap with this object to see if contamination is reduced!\n")

# Marker-panel heatmap on the cluster-filtered object
plot_marker_panel_heatmap(cluster_filtered_seurat, label = "ClusterFiltered")
