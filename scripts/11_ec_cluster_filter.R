## Endothelial branch: keep the EC clusters that are bona fide ECs
## Stage 11 of the GSE159354 EC / immune pipeline.
##
## Inputs:  EC_subtypes_assigned_all_groups.rds (stage 10)
## Outputs: sce_EC_filtered_clusters.rds, *_Filtered_UMAP.pdf, EC_*_filtered.csv, EC_filtering_comparison_summary.csv
## EC clusters 0, 1, 2, 3, 6, 7, 10 and 20 are kept; the rest were judged
## contaminating on the stage-10 cluster / subtype table.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

sce_EC_work <- readRDS("EC_subtypes_assigned_all_groups.rds")

# =============================================================================
# EC ANALYSIS - FILTERED CLUSTERS ONLY (0,1,2,3,6,7,10,20)
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Define colors
group_colors <- c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")
ec_subtype_colors <- c("ArtECs" = "#E74C3C", "CapECs" = "#9B59B6",
                      "VenECs" = "#27AE60", "LECs" = "#F1C40F")

# =============================================================================
# 1. IDENTIFY AND FILTER CLUSTERS
# =============================================================================

cat("=== FILTERING TO KEEP ONLY SPECIFIC CLUSTERS ===\n")

# Clusters to keep
clusters_to_keep <- c("0", "1", "2", "3", "6", "7", "10", "20")
cat("Clusters to keep:", paste(clusters_to_keep, collapse = ", "), "\n")

# Identify the clustering column to use
if("seurat_clusters" %in% colnames(sce_EC_work@meta.data)) {
  cluster_column <- "seurat_clusters"
} else if("RNA_snn_res.0.8" %in% colnames(sce_EC_work@meta.data)) {
  cluster_column <- "RNA_snn_res.0.8"
} else {
  res_cols <- grep("RNA_snn_res\\.", colnames(sce_EC_work@meta.data), value = TRUE)
  cluster_column <- res_cols[1]
}

cat("Using clustering column:", cluster_column, "\n")

# Check which clusters exist in the data
existing_clusters <- unique(as.character(sce_EC_work@meta.data[[cluster_column]]))
cat("Existing clusters:", paste(sort(existing_clusters), collapse = ", "), "\n")

# Find which requested clusters actually exist
available_clusters <- intersect(clusters_to_keep, existing_clusters)
missing_clusters <- setdiff(clusters_to_keep, existing_clusters)

cat("Available clusters to keep:", paste(available_clusters, collapse = ", "), "\n")
if(length(missing_clusters) > 0) {
  cat("Missing clusters (not in data):", paste(missing_clusters, collapse = ", "), "\n")
}

# Filter the data to keep only specified clusters
cat("Filtering data...\n")
cells_to_keep <- sce_EC_work@meta.data[[cluster_column]] %in% available_clusters
sce_EC_filtered <- sce_EC_work[, cells_to_keep]

cat("Original EC cells:", ncol(sce_EC_work), "\n")
cat("Filtered EC cells:", ncol(sce_EC_filtered), "\n")
cat("Cells removed:", ncol(sce_EC_work) - ncol(sce_EC_filtered), "\n")

# Show cluster composition after filtering
cat("\nFiltered cluster composition:\n")
filtered_cluster_table <- table(sce_EC_filtered@meta.data[[cluster_column]])
print(filtered_cluster_table)

# =============================================================================
# 2. RE-ASSIGN EC SUBTYPES FOR FILTERED DATA
# =============================================================================

cat("\n=== RE-ASSIGNING EC SUBTYPES FOR FILTERED DATA ===\n")

# Check if EC_subtype column exists, if not, reassign
if(!"EC_subtype" %in% colnames(sce_EC_filtered@meta.data)) {
  cat("EC_subtype column not found, reassigning subtypes...\n")

  # Use the same module score approach as before
  arterial_markers <- c("GJA5","EFNB2","DLL4","SOX17","HEY2","NOTCH4","BMX","SEMA3G")
  venous_markers <- c("NR2F2","EPHB4","ACKR1","VCAM1","SELP","SELE","VWF","PLVAP")
  capillary_markers <- c("CA4","EDNRB","EMCN","ESAM","KDR","RGCC","AQP1","HSPG2")
  lymphatic_markers <- c("PROX1","LYVE1","PDPN","FLT4","VEGFC","CCL21","RELN","TBX1","FOXC2","NFATC1")

  # Calculate module scores (simplified assignment)
  available_arterial <- intersect(arterial_markers, rownames(sce_EC_filtered))
  available_venous <- intersect(venous_markers, rownames(sce_EC_filtered))
  available_capillary <- intersect(capillary_markers, rownames(sce_EC_filtered))
  available_lymphatic <- intersect(lymphatic_markers, rownames(sce_EC_filtered))

  if(length(available_arterial) > 0) {
    sce_EC_filtered <- AddModuleScore(sce_EC_filtered, features = list(available_arterial), name = "Art_Score", verbose = FALSE)
  }
  if(length(available_venous) > 0) {
    sce_EC_filtered <- AddModuleScore(sce_EC_filtered, features = list(available_venous), name = "Ven_Score", verbose = FALSE)
  }
  if(length(available_capillary) > 0) {
    sce_EC_filtered <- AddModuleScore(sce_EC_filtered, features = list(available_capillary), name = "Cap_Score", verbose = FALSE)
  }
  if(length(available_lymphatic) > 0) {
    sce_EC_filtered <- AddModuleScore(sce_EC_filtered, features = list(available_lymphatic), name = "Lym_Score", verbose = FALSE)
  }

  # Assign subtypes based on highest score
  score_columns <- grep("Score1$", colnames(sce_EC_filtered@meta.data), value = TRUE)
  if(length(score_columns) > 0) {
    score_matrix <- as.matrix(sce_EC_filtered@meta.data[, score_columns])
    subtype_names <- c("ArtECs", "VenECs", "CapECs", "LECs")[1:ncol(score_matrix)]
    colnames(score_matrix) <- subtype_names

    max_scores <- apply(score_matrix, 1, which.max)
    sce_EC_filtered$EC_subtype <- colnames(score_matrix)[max_scores]
  }
}

cat("EC subtype composition in filtered data:\n")
print(table(sce_EC_filtered$EC_subtype))

# =============================================================================
# 3. UMAP VISUALIZATIONS - FILTERED DATA
# =============================================================================

# Set identities
Idents(sce_EC_filtered) <- sce_EC_filtered$EC_subtype

cat("\n=== CREATING UMAPS FOR FILTERED DATA ===\n")

# 1. All groups UMAP
p1_filtered <- DimPlot(sce_EC_filtered, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 0.8) +
  ggtitle("EC Subtypes - Filtered Clusters (All Groups)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p1_filtered)

# 2. CTRL only
sce_EC_ctrl_filtered <- subset(sce_EC_filtered, subset = group == "CTRL")
cat("CTRL ECs (filtered):", ncol(sce_EC_ctrl_filtered), "\n")
print(table(sce_EC_ctrl_filtered$EC_subtype))

p2_filtered <- DimPlot(sce_EC_ctrl_filtered, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 1.0) +
  ggtitle("EC Subtypes - Filtered Clusters (CTRL Only)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p2_filtered)

# 3. IPF only
sce_EC_ipf_filtered <- subset(sce_EC_filtered, subset = group == "IPF")
cat("IPF ECs (filtered):", ncol(sce_EC_ipf_filtered), "\n")
print(table(sce_EC_ipf_filtered$EC_subtype))

p3_filtered <- DimPlot(sce_EC_ipf_filtered, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 1.0) +
  ggtitle("EC Subtypes - Filtered Clusters (IPF Only)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p3_filtered)

# 4. SSc-ILD only
sce_EC_ssc_filtered <- subset(sce_EC_filtered, subset = group == "SSc-ILD")
cat("SSc-ILD ECs (filtered):", ncol(sce_EC_ssc_filtered), "\n")
print(table(sce_EC_ssc_filtered$EC_subtype))

p4_filtered <- DimPlot(sce_EC_ssc_filtered, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 1.0) +
  ggtitle("EC Subtypes - Filtered Clusters (SSc-ILD Only)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p4_filtered)

# 5. Diseased (IPF + SSc-ILD)
sce_EC_diseased_filtered <- subset(sce_EC_filtered, subset = group %in% c("IPF", "SSc-ILD"))
cat("Diseased ECs (filtered):", ncol(sce_EC_diseased_filtered), "\n")
print(table(sce_EC_diseased_filtered$EC_subtype))

p5_filtered <- DimPlot(sce_EC_diseased_filtered, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 1.0) +
  ggtitle("EC Subtypes - Filtered Clusters (Diseased)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p5_filtered)

# =============================================================================
# 4. PERCENTAGE COMPOSITION ANALYSIS - FILTERED DATA
# =============================================================================

cat("\n=== CREATING PERCENTAGE COMPOSITION GRAPHS ===\n")

# Cross-tabulation of EC subtype vs group (filtered data)
proportion_data_filtered <- table(sce_EC_filtered$EC_subtype, sce_EC_filtered$group)
cat("EC subtype composition by group (filtered):\n")
print(proportion_data_filtered)

# Convert to percentages within each EC subtype
proportion_matrix_filtered <- prop.table(proportion_data_filtered, margin = 1) * 100
cat("Percentages within each EC subtype (filtered):\n")
print(round(proportion_matrix_filtered, 1))

# Convert to long format for plotting
proportion_df_filtered <- as.data.frame(proportion_matrix_filtered)
colnames(proportion_df_filtered) <- c("EC_Subtype", "Group", "Percentage")

# Stacked percentage bar plot
p6_filtered <- ggplot(proportion_df_filtered, aes(x = EC_Subtype, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Percentage of Sample Groups within Each EC Subtype (Filtered)",
       x = "EC Subtype", y = "Percentage", fill = "Sample Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold")

print(p6_filtered)

# Dodged percentage bar plot
p6_dodged_filtered <- ggplot(proportion_df_filtered, aes(x = EC_Subtype, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Percentage of Sample Groups within Each EC Subtype (Side-by-Side, Filtered)",
       x = "EC Subtype", y = "Percentage", fill = "Sample Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3)

print(p6_dodged_filtered)

# =============================================================================
# 5. COMBINED LAYOUTS - FILTERED DATA
# =============================================================================

cat("\n=== CREATING COMBINED LAYOUTS ===\n")

# All UMAPs combined
combined_umaps_filtered <- (p1_filtered + p2_filtered) / (p3_filtered + p4_filtered) / (p5_filtered + p6_filtered)
print(combined_umaps_filtered)

# UMAPs only
umaps_only_filtered <- (p1_filtered + p2_filtered) / (p3_filtered + p4_filtered) / p5_filtered
print(umaps_only_filtered)

# =============================================================================
# 6. COMPARISON: BEFORE VS AFTER FILTERING
# =============================================================================

cat("\n=== BEFORE VS AFTER FILTERING COMPARISON ===\n")

cat("Cell counts before filtering:\n")
print(table(sce_EC_work$group))
cat("EC subtype counts before filtering:\n")
print(table(sce_EC_work$EC_subtype))

cat("\nCell counts after filtering:\n")
print(table(sce_EC_filtered$group))
cat("EC subtype counts after filtering:\n")
print(table(sce_EC_filtered$EC_subtype))

cat("\nClusters kept:", paste(available_clusters, collapse = ", "), "\n")
if(length(missing_clusters) > 0) {
  cat("Requested clusters not found:", paste(missing_clusters, collapse = ", "), "\n")
}

# =============================================================================
# 7. SAVE ALL RESULTS - FILTERED DATA
# =============================================================================

cat("\n=== SAVING FILTERED DATA RESULTS ===\n")

# Save the filtered object
saveRDS(sce_EC_filtered, "sce_EC_filtered_clusters.rds")

# Save individual plots
ggsave("1_EC_Subtypes_All_Groups_Filtered_UMAP.pdf", p1_filtered, width = 10, height = 8)
ggsave("2_EC_Subtypes_CTRL_Only_Filtered_UMAP.pdf", p2_filtered, width = 10, height = 8)
ggsave("3_EC_Subtypes_IPF_Only_Filtered_UMAP.pdf", p3_filtered, width = 10, height = 8)
ggsave("4_EC_Subtypes_SSc-ILD_Only_Filtered_UMAP.pdf", p4_filtered, width = 10, height = 8)
ggsave("5_EC_Subtypes_Diseased_Only_Filtered_UMAP.pdf", p5_filtered, width = 10, height = 8)
ggsave("6_EC_Subtypes_Percentages_Stacked_Filtered.pdf", p6_filtered, width = 12, height = 8)
ggsave("6_EC_Subtypes_Percentages_Dodged_Filtered.pdf", p6_dodged_filtered, width = 12, height = 8)

# Save combined layouts
ggsave("EC_Subtypes_All_Visualizations_Filtered_Combined.pdf", combined_umaps_filtered, width = 16, height = 18)
ggsave("EC_Subtypes_UMAPs_Only_Filtered_Combined.pdf", umaps_only_filtered, width = 16, height = 12)

# Save data tables
write.csv(proportion_data_filtered, "EC_subtype_counts_by_group_filtered.csv")
write.csv(proportion_matrix_filtered, "EC_subtype_percentages_by_group_filtered.csv")
write.csv(filtered_cluster_table, "EC_filtered_cluster_composition.csv")

# Create summary comparison table
comparison_summary <- data.frame(
  Metric = c("Total_Cells", "CTRL_Cells", "IPF_Cells", "SSc-ILD_Cells",
            "ArtECs", "CapECs", "VenECs", "LECs"),
  Before_Filtering = c(ncol(sce_EC_work),
                      sum(sce_EC_work$group == "CTRL"),
                      sum(sce_EC_work$group == "IPF"),
                      sum(sce_EC_work$group == "SSc-ILD"),
                      sum(sce_EC_work$EC_subtype == "ArtECs"),
                      sum(sce_EC_work$EC_subtype == "CapECs"),
                      sum(sce_EC_work$EC_subtype == "VenECs"),
                      sum(sce_EC_work$EC_subtype == "LECs")),
  After_Filtering = c(ncol(sce_EC_filtered),
                     sum(sce_EC_filtered$group == "CTRL"),
                     sum(sce_EC_filtered$group == "IPF"),
                     sum(sce_EC_filtered$group == "SSc-ILD"),
                     sum(sce_EC_filtered$EC_subtype == "ArtECs"),
                     sum(sce_EC_filtered$EC_subtype == "CapECs"),
                     sum(sce_EC_filtered$EC_subtype == "VenECs"),
                     sum(sce_EC_filtered$EC_subtype == "LECs"))
)
comparison_summary$Cells_Removed <- comparison_summary$Before_Filtering - comparison_summary$After_Filtering
comparison_summary$Percent_Retained <- round(100 * comparison_summary$After_Filtering / comparison_summary$Before_Filtering, 1)

write.csv(comparison_summary, "EC_filtering_comparison_summary.csv", row.names = FALSE)

cat("Files saved for filtered analysis:\n")
cat("- sce_EC_filtered_clusters.rds\n")
cat("- All UMAP PDFs (1-5) with 'Filtered' in name\n")
cat("- Percentage bar plots with percentages labeled\n")
cat("- Combined layouts\n")
cat("- Data tables and comparison summary\n")

cat("\n=== FILTERED CLUSTER ANALYSIS COMPLETE ===\n")
cat("Clusters kept:", paste(available_clusters, collapse = ", "), "\n")
cat("Total cells retained:", ncol(sce_EC_filtered), "out of", ncol(sce_EC_work), "\n")
