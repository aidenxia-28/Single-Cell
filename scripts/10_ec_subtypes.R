## Endothelial branch: re-cluster ECs and assign arterial / venous / capillary / lymphatic subtypes
## Stage 10 of the GSE159354 EC / immune pipeline.
##
## Inputs:  GSE159354_ThreeGroup_EC-Immune_harmony.rds (stage 01)
## Outputs: EC_subtypes_assigned_all_groups.rds, 1-6_EC_Subtypes_*.pdf, EC_subtype_*.csv
## ECs are taken from the integrated object, re-normalised (2000 HVGs, 50 PCs,
## Harmony on donor + group, Louvain at resolution 0.8) and labelled by the
## highest of four AddModuleScore signatures. The second half draws the UMAP
## with the original cluster numbers and tabulates cluster -> subtype, which is
## what stage 11 uses to decide which clusters to keep.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

combined_seurat <- readRDS("GSE159354_ThreeGroup_EC-Immune_harmony.rds")

# Extract only ECs into a NEW object so your existing sce_EC is untouched
sce_EC_work <- subset(combined_seurat, subset = cellset == "EC")

cat("Total ECs extracted:", ncol(sce_EC_work), "\n")
cat("EC counts by group:\n")
print(table(sce_EC_work$group))

# =============================================================================
# RE-PROCESS ECs FOR SUBCLUSTERING
# =============================================================================

cat("\n=== RE-PROCESSING ECs FOR SUBCLUSTERING ===\n")

# Clear existing reductions
sce_EC_work@reductions <- list()

# Normalize, find variable features, scale
sce_EC_work <- NormalizeData(sce_EC_work, verbose = FALSE)
sce_EC_work <- FindVariableFeatures(sce_EC_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce_EC_work <- ScaleData(sce_EC_work, verbose = FALSE)

# Run PCA
sce_EC_work <- RunPCA(sce_EC_work, npcs = 50, verbose = FALSE)

# Run Harmony for batch correction
sce_EC_work <- RunHarmony(sce_EC_work, group.by.vars = c("donor", "group"), verbose = FALSE)

# Run UMAP
sce_EC_work <- RunUMAP(sce_EC_work, reduction = "harmony", dims = 1:30, verbose = FALSE)

# Find neighbors and clusters
sce_EC_work <- FindNeighbors(sce_EC_work, reduction = "harmony", dims = 1:30, verbose = FALSE)
sce_EC_work <- FindClusters(sce_EC_work, resolution = 0.8, verbose = FALSE)

cat("Initial EC clustering complete. Clusters found:", length(unique(Idents(sce_EC_work))), "\n")

# =============================================================================
# DEFINE EC SUBTYPE MARKER GENES
# =============================================================================

# Define marker genes for each EC subtype

arterial_markers <- c("GJA5","EFNB2","DLL4","SOX17","HEY2","NOTCH4","BMX","SEMA3G")
venous_markers <- c("NR2F2","EPHB4","ACKR1","VCAM1","SELP","SELE","VWF","PLVAP")
capillary_markers <- c("CA4","EDNRB","EMCN","ESAM","KDR","RGCC","AQP1","HSPG2")
lymphatic_markers <- c("PROX1","LYVE1","PDPN","FLT4","VEGFC","CCL21","RELN","TBX1","FOXC2","NFATC1")

# =============================================================================
# EC Subtype Assignment and 6 Visualizations
# =============================================================================

# =============================================================================
# EC Subtype Assignment and 6 Visualizations
# =============================================================================

# Define colors for groups and EC subtypes
group_colors <- c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")
ec_subtype_colors <- c("ArtECs" = "#E74C3C",   # Bright red
                      "CapECs" = "#9B59B6",   # Purple
                      "VenECs" = "#27AE60",   # Green
                      "LECs" = "#F1C40F")     # Yellow

# =============================================================================
# ASSIGN EC SUBTYPES USING MODULE SCORES
# =============================================================================


# Check which markers are available
all_markers <- c(arterial_markers, venous_markers, capillary_markers, lymphatic_markers)
available_markers <- intersect(all_markers, rownames(sce_EC_work))
cat("Available marker genes:", length(available_markers), "out of", length(all_markers), "\n")

missing_markers <- setdiff(all_markers, rownames(sce_EC_work))
if(length(missing_markers) > 0) {
  cat("Missing markers:", paste(missing_markers, collapse = ", "), "\n")
}

# Calculate module scores for each EC subtype
cat("Calculating module scores...\n")

# Arterial ECs
available_arterial <- intersect(arterial_markers, rownames(sce_EC_work))
if(length(available_arterial) > 0) {
  sce_EC_work <- AddModuleScore(sce_EC_work, features = list(available_arterial),
                               name = "Arterial_Score", verbose = FALSE)
  cat("Arterial markers used:", paste(available_arterial, collapse = ", "), "\n")
}

# Venous ECs
available_venous <- intersect(venous_markers, rownames(sce_EC_work))
if(length(available_venous) > 0) {
  sce_EC_work <- AddModuleScore(sce_EC_work, features = list(available_venous),
                               name = "Venous_Score", verbose = FALSE)
  cat("Venous markers used:", paste(available_venous, collapse = ", "), "\n")
}

# Capillary ECs
available_capillary <- intersect(capillary_markers, rownames(sce_EC_work))
if(length(available_capillary) > 0) {
  sce_EC_work <- AddModuleScore(sce_EC_work, features = list(available_capillary),
                               name = "Capillary_Score", verbose = FALSE)
  cat("Capillary markers used:", paste(available_capillary, collapse = ", "), "\n")
}

# Lymphatic ECs
available_lymphatic <- intersect(lymphatic_markers, rownames(sce_EC_work))
if(length(available_lymphatic) > 0) {
  sce_EC_work <- AddModuleScore(sce_EC_work, features = list(available_lymphatic),
                               name = "Lymphatic_Score", verbose = FALSE)
  cat("Lymphatic markers used:", paste(available_lymphatic, collapse = ", "), "\n")
}

# =============================================================================
# ASSIGN EC SUBTYPES BASED ON HIGHEST MODULE SCORE
# =============================================================================

cat("Assigning EC subtypes based on highest module scores...\n")

# Get the score columns
score_columns <- grep("Score1$", colnames(sce_EC_work@meta.data), value = TRUE)
cat("Score columns found:", paste(score_columns, collapse = ", "), "\n")
cat("Number of score columns:", length(score_columns), "\n")

# Debug: check all column names in metadata
cat("All metadata columns containing 'Score':\n")
score_cols_all <- grep("Score", colnames(sce_EC_work@meta.data), value = TRUE)
print(score_cols_all)

if(length(score_columns) > 0) {
  # Create score matrix
  score_matrix <- as.matrix(sce_EC_work@meta.data[, score_columns])
  cat("Score matrix dimensions:", dim(score_matrix), "\n")

  # Assign appropriate column names based on actual number of columns
  subtype_names <- c("ArtECs", "VenECs", "CapECs", "LECs")

  if(ncol(score_matrix) == 4) {
    colnames(score_matrix) <- subtype_names
  } else if(ncol(score_matrix) == 3) {
    colnames(score_matrix) <- subtype_names[1:3]
    cat("Warning: Only 3 module scores found, missing one subtype\n")
  } else if(ncol(score_matrix) == 2) {
    colnames(score_matrix) <- subtype_names[1:2]
    cat("Warning: Only 2 module scores found, missing two subtypes\n")
  } else {
    colnames(score_matrix) <- subtype_names[1:ncol(score_matrix)]
    cat("Warning: Unexpected number of module scores:", ncol(score_matrix), "\n")
  }

  cat("Final score matrix column names:", colnames(score_matrix), "\n")

  # Assign subtypes based on highest score
  max_scores <- apply(score_matrix, 1, which.max)
  sce_EC_work$EC_subtype <- colnames(score_matrix)[max_scores]

  cat("EC subtype assignment complete:\n")
  print(table(sce_EC_work$EC_subtype))

} else {
  cat("Error: No module score columns found\n")
  cat("Attempting alternative assignment based on individual marker expression...\n")

  # Fallback: use individual marker expression
  # Calculate mean expression for each marker set
  arterial_exp <- if(length(available_arterial) > 0) {
    Matrix::rowMeans(GetAssayData(sce_EC_work, slot = "data")[available_arterial, , drop=FALSE])
  } else { rep(0, ncol(sce_EC_work)) }

  venous_exp <- if(length(available_venous) > 0) {
    Matrix::rowMeans(GetAssayData(sce_EC_work, slot = "data")[available_venous, , drop=FALSE])
  } else { rep(0, ncol(sce_EC_work)) }

  capillary_exp <- if(length(available_capillary) > 0) {
    Matrix::rowMeans(GetAssayData(sce_EC_work, slot = "data")[available_capillary, , drop=FALSE])
  } else { rep(0, ncol(sce_EC_work)) }

  lymphatic_exp <- if(length(available_lymphatic) > 0) {
    Matrix::rowMeans(GetAssayData(sce_EC_work, slot = "data")[available_lymphatic, , drop=FALSE])
  } else { rep(0, ncol(sce_EC_work)) }

  # Create expression matrix
  exp_matrix <- cbind(arterial_exp, venous_exp, capillary_exp, lymphatic_exp)
  colnames(exp_matrix) <- c("ArtECs", "VenECs", "CapECs", "LECs")

  # Assign subtypes based on highest expression
  max_exp <- apply(exp_matrix, 1, which.max)
  sce_EC_work$EC_subtype <- colnames(exp_matrix)[max_exp]

  cat("Fallback EC subtype assignment complete:\n")
  print(table(sce_EC_work$EC_subtype))
}

# Set identities
Idents(sce_EC_work) <- sce_EC_work$EC_subtype

# =============================================================================
# 1. UMAP OF ALL ECs (ALL GROUPS) WITH SUBTYPES COLORED
# =============================================================================

cat("\n=== CREATING VISUALIZATION 1: ALL ECs, ALL GROUPS ===\n")

p1 <- DimPlot(sce_EC_work, reduction = "umap", group.by = "EC_subtype",
              cols = ec_subtype_colors, label = TRUE, label.size = 4, pt.size = 0.8) +
  ggtitle("EC Subtypes - All Groups (CTRL, IPF, SSc-ILD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p1)

# =============================================================================
# 2. UMAP OF ECs FROM CTRL ONLY
# =============================================================================

cat("\n=== CREATING VISUALIZATION 2: ECs FROM CTRL ONLY ===\n")

sce_EC_ctrl <- subset(sce_EC_work, subset = group == "CTRL")
cat("CTRL ECs:", ncol(sce_EC_ctrl), "\n")
print(table(sce_EC_ctrl$EC_subtype))

p2 <- DimPlot(sce_EC_ctrl, reduction = "umap", group.by = "EC_subtype",
              cols = ec_subtype_colors, label = TRUE, label.size = 4, pt.size = 1.0) +
  ggtitle("EC Subtypes - CTRL Only") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p2)

# =============================================================================
# 3. UMAP OF ECs FROM IPF ONLY
# =============================================================================

cat("\n=== CREATING VISUALIZATION 3: ECs FROM IPF ONLY ===\n")

sce_EC_ipf <- subset(sce_EC_work, subset = group == "IPF")
cat("IPF ECs:", ncol(sce_EC_ipf), "\n")
print(table(sce_EC_ipf$EC_subtype))

p3 <- DimPlot(sce_EC_ipf, reduction = "umap", group.by = "EC_subtype",
              cols = ec_subtype_colors, label = TRUE, label.size = 4, pt.size = 1.0) +
  ggtitle("EC Subtypes - IPF Only") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p3)

# =============================================================================
# 4. UMAP OF ECs FROM SSc-ILD ONLY
# =============================================================================

cat("\n=== CREATING VISUALIZATION 4: ECs FROM SSc-ILD ONLY ===\n")

sce_EC_ssc <- subset(sce_EC_work, subset = group == "SSc-ILD")
cat("SSc-ILD ECs:", ncol(sce_EC_ssc), "\n")
print(table(sce_EC_ssc$EC_subtype))

p4 <- DimPlot(sce_EC_ssc, reduction = "umap", group.by = "EC_subtype",
              cols = ec_subtype_colors, label = TRUE, label.size = 4, pt.size = 1.0) +
  ggtitle("EC Subtypes - SSc-ILD Only") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p4)

# =============================================================================
# 5. UMAP OF ECs FROM DISEASED (IPF + SSc-ILD)
# =============================================================================

cat("\n=== CREATING VISUALIZATION 5: ECs FROM DISEASED (IPF + SSc-ILD) ===\n")

sce_EC_diseased <- subset(sce_EC_work, subset = group %in% c("IPF", "SSc-ILD"))
cat("Diseased ECs:", ncol(sce_EC_diseased), "\n")
print(table(sce_EC_diseased$EC_subtype))

p5 <- DimPlot(sce_EC_diseased, reduction = "umap", group.by = "EC_subtype",
              cols = ec_subtype_colors, label = TRUE, label.size = 4, pt.size = 1.0) +
  ggtitle("EC Subtypes - Diseased (IPF + SSc-ILD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p5)

# =============================================================================
# 6. BAR GRAPH SHOWING PROPORTION OF EC SUBTYPES BY SAMPLE GROUP
# =============================================================================

cat("\n=== CREATING VISUALIZATION 6: PROPORTION BAR GRAPH ===\n")

# Create cross-tabulation of EC subtype vs group
proportion_data <- table(sce_EC_work$EC_subtype, sce_EC_work$group)
cat("EC subtype composition by group:\n")
print(proportion_data)

# Convert to proportions within each EC subtype
proportion_matrix <- prop.table(proportion_data, margin = 1)
cat("Proportions within each EC subtype:\n")
print(round(proportion_matrix, 3))

# Convert to long format for plotting
proportion_df <- as.data.frame(proportion_matrix)
colnames(proportion_df) <- c("EC_Subtype", "Group", "Proportion")

# Create the bar graph
p6 <- ggplot(proportion_df, aes(x = EC_Subtype, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Proportion of Sample Groups within Each EC Subtype",
       x = "EC Subtype", y = "Proportion", fill = "Sample Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(labels = scales::percent_format())

print(p6)

# Alternative version with bars side by side (dodged)
p6_dodged <- ggplot(proportion_df, aes(x = EC_Subtype, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Proportion of Sample Groups within Each EC Subtype (Side-by-Side)",
       x = "EC Subtype", y = "Proportion", fill = "Sample Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(labels = scales::percent_format())

print(p6_dodged)

# =============================================================================
# CREATE COMBINED LAYOUTS
# =============================================================================

cat("\n=== CREATING COMBINED LAYOUTS ===\n")

# Layout 1: All individual UMAPs
combined_umaps <- (p1 + p2) / (p3 + p4) / (p5 + p6)
print(combined_umaps)

# Layout 2: UMAPs only
umaps_only <- (p1 + p2) / (p3 + p4) / p5
print(umaps_only)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== EC SUBTYPE ANALYSIS SUMMARY ===\n")

cat("Total ECs analyzed:", ncol(sce_EC_work), "\n")
cat("EC subtypes identified:", paste(unique(sce_EC_work$EC_subtype), collapse = ", "), "\n")

# Detailed breakdown by group
cat("\nDetailed breakdown by group:\n")
for(group in c("CTRL", "IPF", "SSc-ILD")) {
  group_ecs <- subset(sce_EC_work, subset = group == group)
  cat("\n", group, ":\n")
  cat("  Total ECs:", ncol(group_ecs), "\n")
  group_table <- table(group_ecs$EC_subtype)
  for(subtype in names(group_table)) {
    cat("  ", subtype, ":", group_table[subtype], "\n")
  }
}

# Overall proportions
cat("\nOverall EC subtype proportions:\n")
overall_props <- prop.table(table(sce_EC_work$EC_subtype))
print(round(overall_props, 3))

# =============================================================================
# SAVE ALL RESULTS
# =============================================================================

cat("\n=== SAVING ALL VISUALIZATIONS ===\n")

# Save the processed EC object
saveRDS(sce_EC_work, "EC_subtypes_assigned_all_groups.rds")

# Save individual plots
ggsave("1_EC_Subtypes_All_Groups_UMAP.pdf", p1, width = 10, height = 8)
ggsave("2_EC_Subtypes_CTRL_Only_UMAP.pdf", p2, width = 10, height = 8)
ggsave("3_EC_Subtypes_IPF_Only_UMAP.pdf", p3, width = 10, height = 8)
ggsave("4_EC_Subtypes_SSc-ILD_Only_UMAP.pdf", p4, width = 10, height = 8)
ggsave("5_EC_Subtypes_Diseased_Only_UMAP.pdf", p5, width = 10, height = 8)
ggsave("6_EC_Subtypes_Proportions_Stacked_Bar.pdf", p6, width = 10, height = 8)
ggsave("6_EC_Subtypes_Proportions_Dodged_Bar.pdf", p6_dodged, width = 10, height = 8)

# Save combined layouts
ggsave("EC_Subtypes_All_Visualizations_Combined.pdf", combined_umaps, width = 16, height = 18)
ggsave("EC_Subtypes_UMAPs_Only_Combined.pdf", umaps_only, width = 16, height = 12)

# Save data tables
write.csv(proportion_data, "EC_subtype_counts_by_group.csv")
write.csv(proportion_matrix, "EC_subtype_proportions_by_group.csv")

# Save summary statistics
summary_stats <- data.frame(
  EC_Subtype = names(table(sce_EC_work$EC_subtype)),
  Total_Cells = as.numeric(table(sce_EC_work$EC_subtype)),
  Proportion = as.numeric(prop.table(table(sce_EC_work$EC_subtype)))
)
write.csv(summary_stats, "EC_subtype_summary_statistics.csv", row.names = FALSE)

cat("Files saved:\n")
cat("- EC_subtypes_assigned_all_groups.rds\n")
cat("- 1_EC_Subtypes_All_Groups_UMAP.pdf\n")
cat("- 2_EC_Subtypes_CTRL_Only_UMAP.pdf\n")
cat("- 3_EC_Subtypes_IPF_Only_UMAP.pdf\n")
cat("- 4_EC_Subtypes_SSc-ILD_Only_UMAP.pdf\n")
cat("- 5_EC_Subtypes_Diseased_Only_UMAP.pdf\n")
cat("- 6_EC_Subtypes_Proportions_Stacked_Bar.pdf\n")
cat("- 6_EC_Subtypes_Proportions_Dodged_Bar.pdf\n")
cat("- EC_Subtypes_All_Visualizations_Combined.pdf\n")
cat("- EC_Subtypes_UMAPs_Only_Combined.pdf\n")
cat("- EC_subtype_counts_by_group.csv\n")
cat("- EC_subtype_proportions_by_group.csv\n")
cat("- EC_subtype_summary_statistics.csv\n")

cat("\n=== ALL 6 VISUALIZATIONS COMPLETE ===\n")


# =============================================================================
# EC UMAP WITH ORIGINAL CLUSTER NUMBERS
# =============================================================================


# =============================================================================
# CREATE UMAP SHOWING ORIGINAL CLUSTER NUMBERS
# =============================================================================

cat("Creating UMAP with original cluster numbers...\n")

# Check what clustering information is available
cat("Available metadata columns with 'cluster' or 'snn':\n")
cluster_cols <- grep("cluster|snn", colnames(sce_EC_work@meta.data), value = TRUE, ignore.case = TRUE)
print(cluster_cols)

# Use the most likely cluster column (try seurat_clusters first, then RNA_snn_res.0.8)
if("seurat_clusters" %in% colnames(sce_EC_work@meta.data)) {
  cluster_column <- "seurat_clusters"
  cat("Using 'seurat_clusters' for cluster numbers\n")
} else if("RNA_snn_res.0.8" %in% colnames(sce_EC_work@meta.data)) {
  cluster_column <- "RNA_snn_res.0.8"
  cat("Using 'RNA_snn_res.0.8' for cluster numbers\n")
} else {
  # Find the first available clustering resolution
  res_cols <- grep("RNA_snn_res\\.", colnames(sce_EC_work@meta.data), value = TRUE)
  if(length(res_cols) > 0) {
    cluster_column <- res_cols[1]
    cat("Using", cluster_column, "for cluster numbers\n")
  } else {
    stop("No clustering information found in metadata")
  }
}

# Check the clusters
cat("Clusters found in", cluster_column, ":\n")
cluster_table <- table(sce_EC_work@meta.data[[cluster_column]])
print(cluster_table)

# Create UMAP with cluster numbers
p1_clusters <- DimPlot(sce_EC_work,
                      reduction = "umap",
                      group.by = cluster_column,
                      label = TRUE,
                      label.size = 6,
                      pt.size = 0.8,
                      repel = TRUE) +
  ggtitle("EC Clusters - All Groups (CTRL, IPF, SSc-ILD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(color = "Cluster")

print(p1_clusters)

# =============================================================================
# ALTERNATIVE VERSION WITH DIFFERENT STYLING
# =============================================================================

# Version with larger labels and no legend (cleaner look)
p1_clusters_clean <- DimPlot(sce_EC_work,
                            reduction = "umap",
                            group.by = cluster_column,
                            label = TRUE,
                            label.size = 7,
                            pt.size = 0.8,
                            repel = TRUE) +
  ggtitle("EC Clusters - All Groups (CTRL, IPF, SSc-ILD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none") +  # Remove legend for cleaner look
  NoLegend()

print(p1_clusters_clean)

# =============================================================================
# COMPARISON: SIDE-BY-SIDE WITH EC SUBTYPES
# =============================================================================

cat("\nCreating side-by-side comparison...\n")

# Recreate the EC subtype plot for comparison
ec_subtype_colors <- c("ArtECs" = "#E74C3C", "CapECs" = "#9B59B6",
                      "VenECs" = "#27AE60", "LECs" = "#F1C40F")

p1_subtypes <- DimPlot(sce_EC_work, reduction = "umap", group.by = "EC_subtype",
                      cols = ec_subtype_colors, label = FALSE, pt.size = 0.8) +
  ggtitle("EC Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Side-by-side comparison
comparison <- p1_clusters_clean + p1_subtypes
print(comparison)

# =============================================================================
# CLUSTER-SUBTYPE MAPPING TABLE
# =============================================================================

cat("\nCluster to EC subtype mapping:\n")
cluster_subtype_table <- table(sce_EC_work@meta.data[[cluster_column]],
                              sce_EC_work$EC_subtype)
print(cluster_subtype_table)

# Show which cluster each subtype predominantly comes from
cat("\nPredominant cluster for each EC subtype:\n")
subtype_clusters <- apply(cluster_subtype_table, 2, function(x) {
  cluster_name <- rownames(cluster_subtype_table)[which.max(x)]
  max_count <- max(x)
  total_count <- sum(x)
  percentage <- round(100 * max_count / total_count, 1)
  return(paste0("Cluster ", cluster_name, " (", max_count, "/", total_count, " = ", percentage, "%)"))
})

for(subtype in names(subtype_clusters)) {
  cat(subtype, ":", subtype_clusters[subtype], "\n")
}
