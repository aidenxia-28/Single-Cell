## Lymphatic branch: subcluster the LECs
## Stage 20 of the GSE159354 EC / immune pipeline.
##
## Inputs:  sce_EC_filtered_clusters.rds (stage 11)
## Outputs: LECs_subclustered.rds, LEC_*_UMAP.pdf, LEC_Composition_Bars.pdf, LEC_subcluster_*.csv
## LECs are re-processed on their own (2000 HVGs, Harmony on donor + group) and
## clustered at resolution 0.7, giving subclusters LEC0-LEC7.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

sce_EC_filtered <- readRDS("sce_EC_filtered_clusters.rds")


# =============================================================================
# LEC SUBCLUSTERING ANALYSIS
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)

# Define colors
group_colors <- c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")

# =============================================================================
# EXTRACT LECS FROM FILTERED OBJECT
# =============================================================================

cat("=== EXTRACTING LECS FOR SUBCLUSTERING ===\n")

# Extract only LECs from the filtered object
sce_LECs <- subset(sce_EC_filtered, subset = EC_subtype == "LECs")

cat("Total LECs extracted:", ncol(sce_LECs), "\n")
cat("LEC counts by group:\n")
print(table(sce_LECs$group))

if(ncol(sce_LECs) < 50) {
  cat("Warning: Very few LECs found. May need to check EC subtype assignments.\n")
  cat("Available EC subtypes in filtered data:\n")
  print(table(sce_EC_filtered$EC_subtype))
  stop("Insufficient LECs for subclustering")
}

# =============================================================================
# RE-PROCESS LECS FOR SUBCLUSTERING
# =============================================================================

cat("\n=== RE-PROCESSING LECS FOR SUBCLUSTERING ===\n")

# Clear existing reductions to start fresh
sce_LECs@reductions <- list()

# Normalize data
cat("Normalizing LEC data...\n")
sce_LECs <- NormalizeData(sce_LECs, verbose = FALSE)

# Find variable features
cat("Finding variable features...\n")
sce_LECs <- FindVariableFeatures(sce_LECs, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Scale data
cat("Scaling LEC data...\n")
sce_LECs <- ScaleData(sce_LECs, verbose = FALSE)

# Run PCA
cat("Running PCA...\n")
n_pcs <- min(50, ncol(sce_LECs) - 1, nrow(sce_LECs) - 1)
sce_LECs <- RunPCA(sce_LECs, npcs = n_pcs, verbose = FALSE)

cat("Using", n_pcs, "principal components\n")

# Run Harmony for batch correction
cat("Running Harmony batch correction...\n")
sce_LECs <- RunHarmony(sce_LECs, group.by.vars = c("donor", "group"), verbose = FALSE)

# Run UMAP
cat("Running UMAP...\n")
sce_LECs <- RunUMAP(sce_LECs, reduction = "harmony", dims = 1:min(30, n_pcs), verbose = FALSE)

# Find neighbors and clusters
cat("Finding neighbors and clusters...\n")
sce_LECs <- FindNeighbors(sce_LECs, reduction = "harmony", dims = 1:min(30, n_pcs), verbose = FALSE)

# Test different resolutions to find optimal clustering
resolutions <- c(0.3, 0.5, 0.7, 0.9, 1.1)
for(res in resolutions) {
  sce_LECs <- FindClusters(sce_LECs, resolution = res, verbose = FALSE)
  n_clusters <- length(unique(Idents(sce_LECs)))
  cat("Resolution", res, ":", n_clusters, "clusters\n")
}

# Choose resolution (adjust as needed - aim for 3-8 clusters)
selected_resolution <- 0.7
sce_LECs <- FindClusters(sce_LECs, resolution = selected_resolution, verbose = FALSE)

# Rename clusters
n_clusters <- length(unique(Idents(sce_LECs)))
cluster_names <- paste0("LEC", 0:(n_clusters-1))
names(cluster_names) <- levels(Idents(sce_LECs))
sce_LECs <- RenameIdents(sce_LECs, cluster_names)
sce_LECs$LEC_subcluster <- Idents(sce_LECs)

cat("LEC subclustering complete. Found", n_clusters, "subclusters:", paste(levels(sce_LECs$LEC_subcluster), collapse = ", "), "\n")

# =============================================================================
# CREATE UMAP VISUALIZATIONS
# =============================================================================

cat("\n=== CREATING LEC SUBCLUSTER VISUALIZATIONS ===\n")

# Generate distinct colors for LEC subclusters
library(RColorBrewer)
if(n_clusters <= 8) {
  lec_colors <- brewer.pal(max(3, n_clusters), "Set2")[1:n_clusters]
} else {
  lec_colors <- rainbow(n_clusters)
}
names(lec_colors) <- levels(sce_LECs$LEC_subcluster)

# 1. UMAP colored by LEC subclusters
p1_lec <- DimPlot(sce_LECs, reduction = "umap", group.by = "LEC_subcluster",
                  cols = lec_colors, label = TRUE, label.size = 5, pt.size = 1.2) +
  ggtitle("LEC Subclusters - All Groups") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p1_lec)

# 2. UMAP colored by disease group
p2_lec <- DimPlot(sce_LECs, reduction = "umap", group.by = "group",
                  cols = group_colors, pt.size = 1.2) +
  ggtitle("LECs by Disease Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p2_lec)

# 3. Split by disease group
p3_lec <- DimPlot(sce_LECs, reduction = "umap", group.by = "LEC_subcluster",
                  split.by = "group", cols = lec_colors,
                  label = TRUE, label.size = 4, pt.size = 1.0, ncol = 3) +
  ggtitle("LEC Subclusters by Disease Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p3_lec)

# =============================================================================
# CLUSTER COMPOSITION ANALYSIS
# =============================================================================

cat("\n=== LEC SUBCLUSTER COMPOSITION ANALYSIS ===\n")

# Show cluster sizes
cat("LEC subcluster sizes:\n")
lec_cluster_table <- table(sce_LECs$LEC_subcluster)
print(lec_cluster_table)

# Show composition by group
cat("LEC subclusters by disease group:\n")
lec_group_table <- table(sce_LECs$LEC_subcluster, sce_LECs$group)
print(lec_group_table)

# Calculate proportions
lec_group_props <- prop.table(lec_group_table, margin = 1) * 100
cat("Percentages within each LEC subcluster:\n")
print(round(lec_group_props, 1))

# Create composition bar plot
lec_prop_df <- as.data.frame(lec_group_props)
colnames(lec_prop_df) <- c("LEC_Subcluster", "Group", "Percentage")

p4_lec <- ggplot(lec_prop_df, aes(x = LEC_Subcluster, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "LEC Subcluster Composition by Disease Group",
       x = "LEC Subcluster", y = "Percentage", fill = "Disease Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  geom_text(aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3)

print(p4_lec)

# =============================================================================
# COMBINED LAYOUTS
# =============================================================================

# Create combined layouts
combined_lec_basic <- p1_lec + p2_lec
print(combined_lec_basic)

combined_lec_full <- (p1_lec + p2_lec) / p4_lec
print(combined_lec_full)

# =============================================================================
# LYMPHATIC MARKER EXPRESSION
# =============================================================================

cat("\n=== LYMPHATIC MARKER EXPRESSION ANALYSIS ===\n")

# Check expression of key lymphatic markers across subclusters
lymphatic_markers <- c("PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "VEGFC")
available_lymphatic <- intersect(lymphatic_markers, rownames(sce_LECs))

cat("Available lymphatic markers:", paste(available_lymphatic, collapse = ", "), "\n")

if(length(available_lymphatic) > 0) {
  # Create feature plots for lymphatic markers
  if(length(available_lymphatic) >= 4) {
    p_markers <- FeaturePlot(sce_LECs, features = available_lymphatic[1:4],
                            ncol = 2, pt.size = 0.8)
    print(p_markers)
  } else {
    p_markers <- FeaturePlot(sce_LECs, features = available_lymphatic,
                            ncol = 2, pt.size = 0.8)
    print(p_markers)
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== LEC SUBCLUSTERING SUMMARY ===\n")

cat("Total LECs analyzed:", ncol(sce_LECs), "\n")
cat("Number of LEC subclusters:", n_clusters, "\n")
cat("Subclusters:", paste(levels(sce_LECs$LEC_subcluster), collapse = ", "), "\n")
cat("Resolution used:", selected_resolution, "\n")

# Show distribution across groups
for(group in c("CTRL", "IPF", "SSc-ILD")) {
  group_lecs <- sum(sce_LECs$group == group)
  group_pct <- round(100 * group_lecs / ncol(sce_LECs), 1)
  cat(group, "LECs:", group_lecs, "(", group_pct, "%)\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING LEC SUBCLUSTERING RESULTS ===\n")

# Save the LEC object
saveRDS(sce_LECs, "LECs_subclustered.rds")

# Save plots
ggsave("LEC_Subclusters_UMAP.pdf", p1_lec, width = 10, height = 8)
ggsave("LEC_by_Group_UMAP.pdf", p2_lec, width = 10, height = 8)
ggsave("LEC_Split_by_Group_UMAP.pdf", p3_lec, width = 15, height = 5)
ggsave("LEC_Composition_Bars.pdf", p4_lec, width = 10, height = 6)
ggsave("LEC_Combined_Basic.pdf", combined_lec_basic, width = 16, height = 8)
ggsave("LEC_Combined_Full.pdf", combined_lec_full, width = 16, height = 12)

if(exists("p_markers")) {
  ggsave("LEC_Lymphatic_Markers.pdf", p_markers, width = 12, height = 10)
}

# Save data tables
write.csv(lec_group_table, "LEC_subcluster_counts_by_group.csv")
write.csv(lec_group_props, "LEC_subcluster_percentages_by_group.csv")

# Save cluster summary
lec_summary <- data.frame(
  LEC_Subcluster = names(lec_cluster_table),
  Cell_Count = as.numeric(lec_cluster_table),
  Percentage = round(100 * as.numeric(lec_cluster_table) / sum(lec_cluster_table), 1)
)
write.csv(lec_summary, "LEC_subcluster_summary.csv", row.names = FALSE)

cat("Files saved:\n")
cat("- LECs_subclustered.rds\n")
cat("- LEC_Subclusters_UMAP.pdf\n")
cat("- LEC_by_Group_UMAP.pdf\n")
cat("- LEC_Split_by_Group_UMAP.pdf\n")
cat("- LEC_Composition_Bars.pdf\n")
cat("- LEC_Combined_Basic.pdf\n")
cat("- LEC_Combined_Full.pdf\n")
cat("- LEC_Lymphatic_Markers.pdf (if markers available)\n")
cat("- CSV files with composition data\n")

cat("\n=== LEC SUBCLUSTERING COMPLETE ===\n")
