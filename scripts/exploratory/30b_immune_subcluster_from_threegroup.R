## Immune branch, alternative start: re-cluster CD45+ cells straight from the integrated object (exploratory)
## Exploratory alternative to stage 30; the main branch starts from the cleaned object.
##
## Inputs:  GSE159354_ThreeGroup_EC-Immune_harmony.rds (stage 01)
## Outputs: immune_cells_clustered.rds, Immune_*.pdf, immune_cluster_*.csv
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# =============================================================================
# Immune Cell Extraction and Clustering Analysis
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# LOAD DATA AND EXTRACT IMMUNE CELLS
# =============================================================================

cat("Loading combined dataset and extracting immune cells...\n")

# Load the combined dataset
combined_seurat <- readRDS("GSE159354_ThreeGroup_EC-Immune_harmony.rds")

cat("Total cells in combined dataset:", ncol(combined_seurat), "\n")
cat("Cell types in dataset:\n")
print(table(combined_seurat$cellset))

# Extract only immune cells
immune_cells <- subset(combined_seurat, subset = cellset == "Immune")

cat("Immune cells extracted:", ncol(immune_cells), "\n")
cat("Immune cells by disease group:\n")
print(table(immune_cells$group))

cat("Immune cells by donor:\n")
print(table(immune_cells$donor))

cat("Immune cells by fraction:\n")
print(table(immune_cells$fraction))

# Check if we have sufficient cells for analysis
if(ncol(immune_cells) < 100) {
  stop("Too few immune cells for meaningful clustering analysis")
}

# =============================================================================
# RE-PROCESS IMMUNE CELLS FOR CLUSTERING
# =============================================================================

cat("\n=== RE-PROCESSING IMMUNE CELLS FOR CLUSTERING ===\n")

# Clear existing reductions to start fresh
immune_cells@reductions <- list()

# Set default assay
DefaultAssay(immune_cells) <- "RNA"

# Join layers for Seurat v5 compatibility
if(packageVersion("Seurat") >= "5.0.0") {
  immune_cells[["RNA"]] <- JoinLayers(immune_cells[["RNA"]])
}

# Normalize data
cat("Normalizing immune cell data...\n")
immune_cells <- NormalizeData(immune_cells, normalization.method = "LogNormalize",
                             scale.factor = 1e4, verbose = FALSE)

# Find variable features
cat("Finding variable features...\n")
immune_cells <- FindVariableFeatures(immune_cells, selection.method = "vst",
                                    nfeatures = 3000, verbose = FALSE)

# Scale data
cat("Scaling immune cell data...\n")
immune_cells <- ScaleData(immune_cells, verbose = FALSE)

# Run PCA
cat("Running PCA...\n")
n_pcs <- min(50, ncol(immune_cells) - 1, nrow(immune_cells) - 1)
immune_cells <- RunPCA(immune_cells, npcs = n_pcs, verbose = FALSE)

cat("Using", n_pcs, "principal components\n")

# =============================================================================
# HARMONY BATCH CORRECTION FOR IMMUNE CELLS
# =============================================================================

cat("Running Harmony batch correction for immune cells...\n")

# Run Harmony to correct for donor and fraction effects
immune_cells <- RunHarmony(immune_cells,
                          group.by.vars = c("donor", "fraction"),
                          verbose = TRUE)

# =============================================================================
# UMAP AND CLUSTERING
# =============================================================================

cat("Running UMAP and clustering...\n")

# Run UMAP
immune_cells <- RunUMAP(immune_cells, reduction = "harmony",
                       dims = 1:min(30, n_pcs), verbose = FALSE)

# Find neighbors
immune_cells <- FindNeighbors(immune_cells, reduction = "harmony",
                             dims = 1:min(30, n_pcs), verbose = FALSE)

# Test different resolutions for clustering
resolutions <- c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3)
cat("Testing different clustering resolutions:\n")

for(res in resolutions) {
  immune_cells <- FindClusters(immune_cells, resolution = res, verbose = FALSE)
  n_clusters <- length(unique(Idents(immune_cells)))
  cat("Resolution", res, ":", n_clusters, "clusters\n")
}

# Choose optimal resolution (adjust based on your preference - aim for 8-15 clusters typically)
selected_resolution <- 0.8
immune_cells <- FindClusters(immune_cells, resolution = selected_resolution, verbose = FALSE)

# Rename clusters for clarity
n_clusters <- length(unique(Idents(immune_cells)))
cluster_names <- paste0("IMM", 0:(n_clusters-1))
names(cluster_names) <- levels(Idents(immune_cells))
immune_cells <- RenameIdents(immune_cells, cluster_names)
immune_cells$immune_cluster <- Idents(immune_cells)

cat("Immune cell clustering complete. Found", n_clusters, "clusters:",
    paste(levels(immune_cells$immune_cluster), collapse = ", "), "\n")

# =============================================================================
# CREATE VISUALIZATIONS
# =============================================================================

cat("\n=== CREATING IMMUNE CELL VISUALIZATIONS ===\n")

# Define colors
group_colors <- c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")

# Generate colors for immune clusters
if(n_clusters <= 12) {
  immune_colors <- brewer.pal(max(3, n_clusters), "Set3")[1:n_clusters]
} else {
  immune_colors <- rainbow(n_clusters)
}
names(immune_colors) <- levels(immune_cells$immune_cluster)

# 1. UMAP colored by immune clusters
p1_immune <- DimPlot(immune_cells, reduction = "umap", group.by = "immune_cluster",
                    cols = immune_colors, label = TRUE, label.size = 5, pt.size = 1.0) +
  ggtitle("Immune Cell Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12))

print(p1_immune)

# 2. UMAP colored by disease group
p2_immune <- DimPlot(immune_cells, reduction = "umap", group.by = "group",
                    cols = group_colors, pt.size = 1.0) +
  ggtitle("Immune Cells by Disease Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12))

print(p2_immune)

# 3. Split by disease group
p3_immune <- DimPlot(immune_cells, reduction = "umap", group.by = "immune_cluster",
                    split.by = "group", cols = immune_colors,
                    label = TRUE, label.size = 4, pt.size = 0.8, ncol = 3) +
  ggtitle("Immune Clusters by Disease Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p3_immune)

# 4. UMAP colored by donor
p4_immune <- DimPlot(immune_cells, reduction = "umap", group.by = "donor",
                    pt.size = 1.0) +
  ggtitle("Immune Cells by Donor") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p4_immune)

# =============================================================================
# CLUSTER COMPOSITION ANALYSIS
# =============================================================================

cat("\n=== IMMUNE CLUSTER COMPOSITION ANALYSIS ===\n")

# Show cluster sizes
cat("Immune cluster sizes:\n")
immune_cluster_table <- table(immune_cells$immune_cluster)
print(immune_cluster_table)

# Show composition by group
cat("Immune clusters by disease group:\n")
immune_group_table <- table(immune_cells$immune_cluster, immune_cells$group)
print(immune_group_table)

# Calculate proportions within each cluster
immune_group_props <- prop.table(immune_group_table, margin = 1) * 100
cat("Percentages within each immune cluster:\n")
print(round(immune_group_props, 1))

# Show composition by donor
cat("Immune clusters by donor:\n")
immune_donor_table <- table(immune_cells$immune_cluster, immune_cells$donor)
print(immune_donor_table)

# Create composition bar plot
immune_prop_df <- as.data.frame(immune_group_props)
colnames(immune_prop_df) <- c("Immune_Cluster", "Group", "Percentage")

p5_immune <- ggplot(immune_prop_df, aes(x = Immune_Cluster, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Immune Cluster Composition by Disease Group",
       x = "Immune Cluster", y = "Percentage", fill = "Disease Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  geom_text(aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3)

print(p5_immune)

# =============================================================================
# IMMUNE MARKER EXPRESSION
# =============================================================================

cat("\n=== IMMUNE MARKER EXPRESSION ANALYSIS ===\n")

# Check expression of key immune markers across clusters
immune_markers <- c(
  # General immune
  "PTPRC", "CD68", "CD3E", "CD4", "CD8A",
  # T cells
  "CD3D", "CD3G", "IL7R", "CCR7", "CD40LG",
  # B cells
  "MS4A1", "CD19", "CD20", "IGHM", "IGHG1",
  # Myeloid
  "CD14", "FCGR3A", "LYZ", "S100A8", "S100A9",
  # NK cells
  "GNLY", "NKG7", "GZMB", "PRF1", "KLRD1",
  # Dendritic cells
  "CLEC4C", "IRF8", "CLEC9A", "CD1C",
  # Macrophages
  "CD68", "CD163", "MRC1", "MARCO", "LYVE1"
)

available_immune_markers <- intersect(immune_markers, rownames(immune_cells))

cat("Available immune markers:", paste(available_immune_markers, collapse = ", "), "\n")

if(length(available_immune_markers) > 0) {
  # Create feature plots for immune markers (first 6)
  if(length(available_immune_markers) >= 6) {
    p_markers <- FeaturePlot(immune_cells, features = available_immune_markers[1:6],
                           ncol = 3, pt.size = 0.5)
    print(p_markers)
  } else {
    p_markers <- FeaturePlot(immune_cells, features = available_immune_markers,
                           ncol = 2, pt.size = 0.5)
    print(p_markers)
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== IMMUNE CLUSTERING SUMMARY ===\n")

cat("Total immune cells analyzed:", ncol(immune_cells), "\n")
cat("Number of immune clusters:", n_clusters, "\n")
cat("Immune clusters:", paste(levels(immune_cells$immune_cluster), collapse = ", "), "\n")
cat("Resolution used:", selected_resolution, "\n")

# Show distribution across groups
for(group in c("CTRL", "IPF", "SSc-ILD")) {
  group_immune <- sum(immune_cells$group == group)
  group_pct <- round(100 * group_immune / ncol(immune_cells), 1)
  cat(group, "immune cells:", group_immune, "(", group_pct, "%)\n")
}

# Show distribution across donors
cat("Distribution by donor:\n")
for(donor in unique(immune_cells$donor)) {
  donor_immune <- sum(immune_cells$donor == donor)
  donor_pct <- round(100 * donor_immune / ncol(immune_cells), 1)
  cat("  ", donor, ":", donor_immune, "(", donor_pct, "%)\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING IMMUNE CLUSTERING RESULTS ===\n")

# Save the immune cell object
saveRDS(immune_cells, "immune_cells_clustered.rds")

# Save plots
ggsave("Immune_Clusters_UMAP.pdf", p1_immune, width = 10, height = 8)
ggsave("Immune_by_Group_UMAP.pdf", p2_immune, width = 10, height = 8)
ggsave("Immune_Split_by_Group_UMAP.pdf", p3_immune, width = 15, height = 5)
ggsave("Immune_by_Donor_UMAP.pdf", p4_immune, width = 12, height = 8)
ggsave("Immune_Composition_Bars.pdf", p5_immune, width = 10, height = 6)

if(exists("p_markers")) {
  ggsave("Immune_Markers_FeaturePlots.pdf", p_markers, width = 15, height = 10)
}

# Save data tables
write.csv(immune_group_table, "immune_cluster_counts_by_group.csv")
write.csv(immune_group_props, "immune_cluster_percentages_by_group.csv")
write.csv(immune_donor_table, "immune_cluster_counts_by_donor.csv")

# Save cluster summary
immune_summary <- data.frame(
  Immune_Cluster = names(immune_cluster_table),
  Cell_Count = as.numeric(immune_cluster_table),
  Percentage = round(100 * as.numeric(immune_cluster_table) / sum(immune_cluster_table), 1)
)
write.csv(immune_summary, "immune_cluster_summary.csv", row.names = FALSE)

# Create detailed metadata summary
immune_metadata_summary <- immune_cells@meta.data %>%
  group_by(immune_cluster, group, donor) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  arrange(immune_cluster, group, donor)

write.csv(immune_metadata_summary, "immune_detailed_metadata_summary.csv", row.names = FALSE)

cat("Files saved:\n")
cat("- immune_cells_clustered.rds\n")
cat("- Immune_Clusters_UMAP.pdf\n")
cat("- Immune_by_Group_UMAP.pdf\n")
cat("- Immune_Split_by_Group_UMAP.pdf\n")
cat("- Immune_by_Donor_UMAP.pdf\n")
cat("- Immune_Composition_Bars.pdf\n")
cat("- Immune_Markers_FeaturePlots.pdf\n")
cat("- CSV files with composition data\n")

# =============================================================================
# IMMUNE CELL TYPE IDENTIFICATION MARKERS
# =============================================================================

cat("\n=== IMMUNE CELL TYPE IDENTIFICATION ===\n")

# Define marker gene sets for major immune cell types
immune_marker_sets <- list(
  "T_cells" = c("CD3D", "CD3E", "CD3G"),
  "CD4_T_cells" = c("CD4", "IL7R", "CCR7"),
  "CD8_T_cells" = c("CD8A", "CD8B", "GZMK"),
  "NK_cells" = c("GNLY", "NKG7", "GZMB", "PRF1", "KLRD1"),
  "B_cells" = c("MS4A1", "CD19", "IGHM", "IGHD"),
  "Plasma_cells" = c("IGHG1", "IGHG3", "MZB1", "SDC1"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Macrophages" = c("CD68", "CD163", "MRC1", "MARCO"),
  "Dendritic_cells" = c("CLEC4C", "IRF8", "CLEC9A", "CD1C"),
  "Neutrophils" = c("FCGR3B", "CEACAM8", "CSF3R")
)

# Calculate average expression for each marker set
marker_expression <- data.frame()

for(cell_type in names(immune_marker_sets)) {
  markers <- immune_marker_sets[[cell_type]]
  available_markers <- intersect(markers, rownames(immune_cells))

  if(length(available_markers) > 0) {
    # Calculate average expression across clusters
    avg_exp <- AverageExpression(immune_cells,
                               features = available_markers,
                               group.by = "immune_cluster",
                               slot = "data")$RNA

    # Average across markers for this cell type
    if(length(available_markers) == 1) {
      cell_type_score <- avg_exp[1, ]
    } else {
      cell_type_score <- colMeans(avg_exp)
    }

    # Add to summary
    for(cluster in names(cell_type_score)) {
      marker_expression <- rbind(marker_expression, data.frame(
        cluster = cluster,
        cell_type = cell_type,
        avg_expression = cell_type_score[[cluster]],
        markers_available = paste(available_markers, collapse = ", ")
      ))
    }
  }
}

# Find predicted cell type for each cluster
if(nrow(marker_expression) > 0) {
  cluster_predictions <- marker_expression %>%
    group_by(cluster) %>%
    slice_max(avg_expression, n = 1) %>%
    select(cluster, predicted_type = cell_type, expression_score = avg_expression)

  cat("Predicted immune cell types by cluster:\n")
  print(cluster_predictions)

  # Save predictions
  write.csv(marker_expression, "immune_cluster_marker_expression.csv", row.names = FALSE)
  write.csv(cluster_predictions, "immune_cluster_predictions.csv", row.names = FALSE)
}

# =============================================================================
# CREATE COMBINED LAYOUT PLOTS
# =============================================================================

# Basic layout
combined_immune_basic <- p1_immune + p2_immune
print(combined_immune_basic)

# Full layout
combined_immune_full <- (p1_immune + p2_immune) / p5_immune
print(combined_immune_full)

# Save combined plots
ggsave("Immune_Combined_Basic.pdf", combined_immune_basic, width = 16, height = 8)
ggsave("Immune_Combined_Full.pdf", combined_immune_full, width = 16, height = 12)

cat("\n=== IMMUNE CELL CLUSTERING COMPLETE ===\n")

# =============================================================================
# NEXT STEPS RECOMMENDATIONS
# =============================================================================

cat("\n=== NEXT STEPS RECOMMENDATIONS ===\n")
cat("IMMUNE CELL ANALYSIS COMPLETE:\n\n")

cat("1. CELL TYPE ANNOTATION:\n")
cat("   - Review predicted cell types based on marker expression\n")
cat("   - Validate predictions with known immune markers\n")
cat("   - Consider manual annotation based on expert knowledge\n\n")

cat("2. DIFFERENTIAL EXPRESSION:\n")
cat("   - Compare immune clusters between disease groups\n")
cat("   - Identify disease-specific immune responses\n")
cat("   - Focus on clusters enriched in IPF or SSc-ILD\n\n")

cat("3. FUNCTIONAL ANALYSIS:\n")
cat("   - Pathway analysis for each immune cluster\n")
cat("   - Cytokine and chemokine expression analysis\n")
cat("   - Activation state assessment\n\n")

cat("4. INTERACTION ANALYSIS:\n")
cat("   - Ligand-receptor analysis between immune clusters and LECs\n")
cat("   - Communication pathway mapping\n")
cat("   - Disease-specific interactions\n\n")

cat("IMMUNE CLUSTERING READY FOR DOWNSTREAM ANALYSIS\n")
