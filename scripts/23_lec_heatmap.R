## Lymphatic branch: LEC subcluster heatmap over a curated gene set
## Stage 23 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: LEC_Comprehensive_Heatmap_8Clusters_*.pdf, LEC_Core_Genes_Heatmap_8Clusters_pheatmap.pdf, LEC_heatmap_*.csv
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# =============================================================================
# LEC Heatmap - 8 Clusters (LEC0-LEC7) with Comprehensive Gene Set using pheatmap
# =============================================================================

library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# =============================================================================
# LOAD LEC SUBCLUSTERED OBJECT
# =============================================================================

cat("Loading LEC subclustered object...\n")
sce_LECs <- readRDS("LECs_subclustered.rds")

cat("Original 8 clusters in dataset:\n")
print(table(sce_LECs$LEC_subcluster))

cat("Total cells in dataset:", ncol(sce_LECs), "\n")
cat("Groups in dataset:", paste(unique(sce_LECs$group), collapse = ", "), "\n")

# =============================================================================
# DEFINE COMPREHENSIVE GENE SET FOR LEC HEATMAP
# =============================================================================

# Core lymphatic identity markers
lymphatic_core <- c("PROX1", "LYVE1", "FLT4", "PDPN", "CCL21", "VEGFC", "MMRN1")

# Lymphatic transcription factors and development
lymphatic_tf <- c("PROX1", "SOX18", "EBF1", "MEOX1", "FOXC2", "GATA2")

# Vessel structure and ECM
vessel_structure <- c("TEK", "CDH5", "PECAM1", "VWF", "COL4A1", "LAMA4", "ITGA9")

# Hypoxia and stress response
hypoxia_stress <- c("HIF1A", "EPAS1", "VEGFA", "ICAM1", "SELP", "TNFSF4")

# Immune and antigen presentation
immune_antigen <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1")

# Chemokines and signaling
chemokine_signaling <- c("CCL2", "CXCL1", "CXCL2", "CCL21", "S1PR1", "ACKR4")

# Cell adhesion and junction
adhesion_junction <- c("PECAM1", "CDH5", "VCAM1", "ICAM1", "OCLN", "TJP1")

# Metabolic and functional
metabolic_functional <- c("FABP4", "ANGPT2", "NOS3", "KDR", "FLT1", "NOTCH1")

# Combine all gene sets
all_heatmap_genes <- c(
  lymphatic_core,
  lymphatic_tf,
  vessel_structure,
  hypoxia_stress,
  immune_antigen,
  chemokine_signaling,
  adhesion_junction,
  metabolic_functional
)

# Remove duplicates and check availability
all_heatmap_genes <- unique(all_heatmap_genes)
available_heatmap_genes <- intersect(all_heatmap_genes, rownames(sce_LECs))

cat("Total genes selected for heatmap:", length(all_heatmap_genes), "\n")
cat("Available genes in dataset:", length(available_heatmap_genes), "\n")
cat("Missing genes:", paste(setdiff(all_heatmap_genes, available_heatmap_genes), collapse = ", "), "\n")

# Create gene categories for annotation
gene_categories <- c()
gene_categories[intersect(lymphatic_core, available_heatmap_genes)] <- "Lymphatic Core"
gene_categories[intersect(lymphatic_tf, available_heatmap_genes)] <- "Transcription Factors"
gene_categories[intersect(vessel_structure, available_heatmap_genes)] <- "Vessel Structure"
gene_categories[intersect(hypoxia_stress, available_heatmap_genes)] <- "Hypoxia/Stress"
gene_categories[intersect(immune_antigen, available_heatmap_genes)] <- "Antigen Presentation"
gene_categories[intersect(chemokine_signaling, available_heatmap_genes)] <- "Chemokine/Signaling"
gene_categories[intersect(adhesion_junction, available_heatmap_genes)] <- "Cell Adhesion"
gene_categories[intersect(metabolic_functional, available_heatmap_genes)] <- "Metabolic/Functional"

# =============================================================================
# DEFINE CLUSTER ORDER
# =============================================================================

# Use custom cluster order: 1, 0, 2, 6, 4, 5, 7, 3
cluster_order <- c("LEC1", "LEC6", "LEC0", "LEC4", "LEC2", "LEC5", "LEC7", "LEC3")
available_clusters <- intersect(cluster_order, unique(sce_LECs$LEC_subcluster))

cat("Custom cluster order for heatmap:", paste(available_clusters, collapse = ", "), "\n")

# Set identities
Idents(sce_LECs) <- sce_LECs$LEC_subcluster

# =============================================================================
# JOIN LAYERS FOR SEURAT V5 COMPATIBILITY
# =============================================================================

# Join layers if using Seurat v5
if(packageVersion("Seurat") >= "5.0.0") {
  sce_LECs[["RNA"]] <- JoinLayers(sce_LECs[["RNA"]])
}

# =============================================================================
# CALCULATE AVERAGE EXPRESSION FOR HEATMAP
# =============================================================================

cat("Calculating average expression for heatmap...\n")

# Calculate average expression
avg_exp_heatmap <- AverageExpression(sce_LECs,
                                     features = available_heatmap_genes,
                                     group.by = "LEC_subcluster",
                                     slot = "data")$RNA

# Ensure all clusters are present
for(cluster in available_clusters) {
  if(!cluster %in% colnames(avg_exp_heatmap)) {
    avg_exp_heatmap[, cluster] <- 0
    cat("Added missing cluster:", cluster, "\n")
  }
}

# Reorder columns to match cluster order
avg_exp_heatmap <- avg_exp_heatmap[, available_clusters, drop = FALSE]

# Scale the data (z-score normalization per gene)
avg_exp_heatmap_scaled <- t(scale(t(avg_exp_heatmap)))

# Check for and handle problematic values
cat("Checking for problematic values...\n")
cat("NA values:", sum(is.na(avg_exp_heatmap_scaled)), "\n")
cat("Inf values:", sum(is.infinite(avg_exp_heatmap_scaled)), "\n")
cat("NaN values:", sum(is.nan(avg_exp_heatmap_scaled)), "\n")

# Remove genes with zero variance (will create NaN after scaling)
gene_vars <- apply(avg_exp_heatmap, 1, var)
zero_var_genes <- names(gene_vars[gene_vars == 0 | is.na(gene_vars)])

if(length(zero_var_genes) > 0) {
  cat("Removing", length(zero_var_genes), "genes with zero variance:", paste(zero_var_genes, collapse = ", "), "\n")
  avg_exp_heatmap <- avg_exp_heatmap[!rownames(avg_exp_heatmap) %in% zero_var_genes, ]
  available_heatmap_genes <- rownames(avg_exp_heatmap)

  # Update gene categories to match remaining genes
  gene_categories <- gene_categories[available_heatmap_genes]
}

# Re-scale the cleaned data
avg_exp_heatmap_scaled <- t(scale(t(avg_exp_heatmap)))

# Replace any remaining problematic values
avg_exp_heatmap_scaled[is.na(avg_exp_heatmap_scaled)] <- 0
avg_exp_heatmap_scaled[is.infinite(avg_exp_heatmap_scaled)] <- 0
avg_exp_heatmap_scaled[is.nan(avg_exp_heatmap_scaled)] <- 0

# Final check
cat("After cleaning - NA values:", sum(is.na(avg_exp_heatmap_scaled)), "\n")
cat("After cleaning - Inf values:", sum(is.infinite(avg_exp_heatmap_scaled)), "\n")
cat("After cleaning - NaN values:", sum(is.nan(avg_exp_heatmap_scaled)), "\n")

cat("Final heatmap matrix dimensions:", dim(avg_exp_heatmap_scaled), "\n")

# =============================================================================
# CREATE ANNOTATIONS FOR PHEATMAP
# =============================================================================

# Calculate group composition for each cluster
cluster_composition <- table(sce_LECs$LEC_subcluster, sce_LECs$group)
cluster_composition_props <- prop.table(cluster_composition, margin = 1)

# Create annotation data frame for clusters (columns)
col_annotation <- data.frame(
  Total_cells = cluster_composition[available_clusters, "CTRL"] +
                cluster_composition[available_clusters, "IPF"] +
                cluster_composition[available_clusters, "SSc-ILD"],
  CTRL_pct = cluster_composition_props[available_clusters, "CTRL"] * 100,
  IPF_pct = cluster_composition_props[available_clusters, "IPF"] * 100,
  SSc_ILD_pct = cluster_composition_props[available_clusters, "SSc-ILD"] * 100
)
rownames(col_annotation) <- available_clusters

# Create annotation data frame for genes (rows) - using remaining genes after cleaning
row_annotation <- data.frame(
  Category = gene_categories[available_heatmap_genes]
)
rownames(row_annotation) <- available_heatmap_genes

# Define annotation colors
category_colors <- brewer.pal(length(unique(row_annotation$Category)), "Set3")
names(category_colors) <- unique(row_annotation$Category)

annotation_colors <- list(
  Category = category_colors,
  Total_cells = colorRampPalette(c("lightblue", "darkblue"))(50),
  CTRL_pct = colorRampPalette(c("white", "#2CA02C"))(50),
  IPF_pct = colorRampPalette(c("white", "#FF7F00"))(50),
  SSc_ILD_pct = colorRampPalette(c("white", "#9467BD"))(50)
)

# =============================================================================
# CREATE PHEATMAP HEATMAP
# =============================================================================

cat("Creating comprehensive LEC heatmap with pheatmap...\n")

# Check if we have enough genes for meaningful analysis
if(nrow(avg_exp_heatmap_scaled) < 3) {
  stop("Not enough genes remaining after cleaning for heatmap analysis")
}

# Define color palette for expression values
color_palette <- colorRampPalette(c("#2166AC", "white", "#762A83"))(100)

# Create the main heatmap with clustering - with error handling
tryCatch({
  pheatmap(
    avg_exp_heatmap_scaled,
    color = color_palette,
    scale = "none",  # already scaled

    # Annotations
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,

    # Clustering - specify distance and linkage methods
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    clustering_distance_cols = "euclidean",

    # Handle potential clustering issues
    show_rownames = TRUE,
    show_colnames = TRUE,

    # Appearance
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 10,
    angle_col = 45,

    # Legends and labels
    main = "LEC Expression Heatmap - 8 Clusters",

    # Size and output
    filename = "LEC_Comprehensive_Heatmap_8Clusters_pheatmap2.pdf",
    width = 12,
    height = 14
  )

  cat("Main clustered heatmap saved successfully\n")
}, error = function(e) {
  cat("Error in clustering, creating heatmap without clustering:\n")
  print(e)

  # Create heatmap without clustering as fallback
  pheatmap(
    avg_exp_heatmap_scaled,
    color = color_palette,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 10,
    angle_col = 45,
    main = "LEC Expression Heatmap - 8 Clusters (No Clustering)",
    filename = "LEC_Comprehensive_Heatmap_8Clusters_pheatmap2.pdf",
    width = 12,
    height = 14
  )
  cat("Fallback heatmap without clustering saved\n")
})

# =============================================================================
# CREATE ORDERED VERSION WITHOUT CLUSTERING
# =============================================================================

cat("Creating ordered heatmap without clustering...\n")

# Order genes by category - using remaining genes after cleaning
ordered_genes <- c()
for(cat in c("Lymphatic Core", "Transcription Factors", "Vessel Structure",
            "Hypoxia/Stress", "Antigen Presentation", "Chemokine/Signaling",
            "Cell Adhesion", "Metabolic/Functional")) {
  cat_genes <- names(gene_categories[gene_categories == cat & !is.na(gene_categories)])
  cat_genes <- intersect(cat_genes, available_heatmap_genes)
  ordered_genes <- c(ordered_genes, cat_genes)
}

# Make sure we have genes to plot
if(length(ordered_genes) == 0) {
  ordered_genes <- available_heatmap_genes
  cat("No genes found in categories, using all available genes\n")
}

# Create ordered heatmap without clustering
pheatmap(
  avg_exp_heatmap_scaled[ordered_genes, available_clusters],
  color = color_palette,
  scale = "none",

  # Turn off clustering
  cluster_rows = FALSE,
  cluster_cols = FALSE,

  # Annotations
  annotation_col = col_annotation,
  annotation_row = row_annotation[ordered_genes, , drop = FALSE],
  annotation_colors = annotation_colors,

  # Appearance
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 10,
  angle_col = 45,

  # Legends and labels
  main = "LEC Expression Heatmap - Ordered by Gene Category",

  # Size and output
  filename = "LEC_Comprehensive_Heatmap_8Clusters_Ordered_pheatmap2.pdf",
  width = 12,
  height = 14
)

cat("Ordered heatmap saved as LEC_Comprehensive_Heatmap_8Clusters_Ordered_pheatmap.pdf\n")

# =============================================================================
# CREATE SIMPLIFIED VERSION WITH CORE GENES ONLY
# =============================================================================

cat("Creating simplified heatmap with core LEC genes...\n")

# Select core lymphatic genes for a cleaner view
core_lec_genes <- c(
  # Core lymphatic
  "PROX1", "LYVE1", "FLT4", "PDPN", "CCL21", "VEGFC",
  # Key TFs
  "SOX18", "EBF1", "MEOX1", "FOXC2",
  # Vessel markers
  "TEK", "CDH5", "PECAM1",
  # Hypoxia/stress
  "HIF1A", "EPAS1", "ICAM1", "SELP",
  # Antigen presentation
  "CD74", "HLA-DRA", "HLA-DRB1",
  # Chemokines
  "CCL2", "CXCL1", "S1PR1",
  # Functional
  "FABP4", "ANGPT2", "KDR"
)

available_core_genes <- intersect(core_lec_genes, available_heatmap_genes)

cat("Available core genes for simplified heatmap:", length(available_core_genes), "\n")

# Only create simplified heatmap if we have core genes
if(length(available_core_genes) > 2) {
  # Create simplified heatmap
  pheatmap(
    avg_exp_heatmap_scaled[available_core_genes, available_clusters],
    color = color_palette,
    scale = "none",

    # Clustering
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    cluster_cols = FALSE,  # Keep cluster order

    # Annotations
    annotation_col = col_annotation,
    annotation_row = row_annotation[available_core_genes, , drop = FALSE],
    annotation_colors = annotation_colors,

    # Appearance
    fontsize = 9,
    fontsize_row = 8,
    fontsize_col = 11,
    angle_col = 45,

    # Legends and labels
    main = "LEC Core Gene Expression - 8 Clusters",

    # Size and output
    filename = "LEC_Core_Genes_Heatmap_8Clusters_pheatmap.pdf",
    width = 10,
    height = 8
  )

  cat("Simplified core genes heatmap saved as LEC_Core_Genes_Heatmap_8Clusters_pheatmap.pdf\n")
} else {
  cat("Not enough core genes available for simplified heatmap\n")
}

# =============================================================================
# CLUSTER-SPECIFIC MARKER ANALYSIS
# =============================================================================

cat("\n=== TOP EXPRESSING CLUSTERS PER GENE ===\n")

# Find which cluster has highest expression for each gene
top_clusters <- apply(avg_exp_heatmap_scaled, 1, function(x) {
  colnames(avg_exp_heatmap_scaled)[which.max(x)]
})

# Group genes by their top expressing cluster
for(cluster in available_clusters) {
  cluster_top_genes <- names(top_clusters[top_clusters == cluster])
  if(length(cluster_top_genes) > 0) {
    cat(cluster, "top genes:", paste(cluster_top_genes, collapse = ", "), "\n")
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== HEATMAP SUMMARY ===\n")

cat("Total cells analyzed:", ncol(sce_LECs), "\n")
cat("Clusters analyzed:", paste(available_clusters, collapse = ", "), "\n")
cat("Genes in final heatmap:", length(available_heatmap_genes), "\n")
cat("Gene categories:", length(unique(gene_categories[!is.na(gene_categories)])), "\n")

# Show cluster sizes
cat("\nCluster sizes:\n")
print(table(sce_LECs$LEC_subcluster))

# Show category breakdown
cat("\nGenes per category (after cleaning):\n")
clean_gene_categories <- gene_categories[!is.na(gene_categories)]
print(table(clean_gene_categories))

# =============================================================================
# SAVE DATA AND RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save expression matrices
write.csv(avg_exp_heatmap, "LEC_heatmap_average_expression.csv")
write.csv(avg_exp_heatmap_scaled, "LEC_heatmap_scaled_expression.csv")

# Save annotations
write.csv(row_annotation, "LEC_heatmap_gene_annotations.csv")
write.csv(col_annotation, "LEC_heatmap_cluster_annotations.csv")

# Save gene lists by category (cleaned data)
clean_gene_categories <- gene_categories[!is.na(gene_categories)]
for(cat in unique(clean_gene_categories)) {
  cat_genes <- names(clean_gene_categories[clean_gene_categories == cat])
  write.csv(data.frame(Gene = cat_genes),
           paste0("LEC_heatmap_genes_", gsub("[/ ]", "_", cat), ".csv"),
           row.names = FALSE)
}

cat("Files saved:\n")
cat("- LEC_Comprehensive_Heatmap_8Clusters_pheatmap.pdf (clustered version)\n")
cat("- LEC_Comprehensive_Heatmap_8Clusters_Ordered_pheatmap.pdf (ordered by category)\n")
cat("- LEC_Core_Genes_Heatmap_8Clusters_pheatmap.pdf (simplified core genes)\n")
cat("- LEC_heatmap_average_expression.csv\n")
cat("- LEC_heatmap_scaled_expression.csv\n")
cat("- LEC_heatmap_gene_annotations.csv\n")
cat("- LEC_heatmap_cluster_annotations.csv\n")
cat("- Gene category files (CSV format)\n")

cat("\n=== LEC HEATMAP ANALYSIS COMPLETE ===\n")
