## Endothelial branch: EC-subtype x disease-group marker heatmap
## Stage 12 of the GSE159354 EC / immune pipeline.
##
## Inputs:  sce_EC_filtered_clusters.rds (stage 11)
## Outputs: EC_Filtered_Pheatmap_Heatmap*.pdf, EC_filtered_heatmap_*.csv
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# =============================================================================
# EC FILTERED DATA PHEATMAP HEATMAP
# =============================================================================

library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

sce_EC_filtered <- readRDS("sce_EC_filtered_clusters.rds")

# =============================================================================
# DEFINE GENES AND CHECK AVAILABILITY
# =============================================================================

cat("=== PREPARING HEATMAP DATA ===\n")

# Genes from your list
heatmap_genes <- c("DKK2", "BMX", "GJA5", "EFNB2", "IGFBP3",
                   "EDNRB", "CA4", "FCN3", "BTNL9", "CX3CL1",
                   "TNFAIP3", "AKAP12", "ICAM1", "HDAC9", "ACKR1",
                   "SELP", "NRG1", "CPE", "LRRC1",
                   "POSTN", "SPRY1", "ZNF385D", "PLVAP", "INSR", "COL15A1", "MKI67", "MMRN1", "PDPN", "PROX1", "EFNA5")

cat("Total genes requested:", length(heatmap_genes), "\n")

# Check which genes are available in the filtered data
available_genes <- intersect(heatmap_genes, rownames(sce_EC_filtered))
missing_genes <- setdiff(heatmap_genes, rownames(sce_EC_filtered))

cat("Available genes:", length(available_genes), "\n")
cat("Missing genes:", length(missing_genes), "\n")

if(length(missing_genes) > 0) {
  cat("Missing genes list:\n")
  for(i in 1:length(missing_genes)) {
    cat("  ", missing_genes[i], "\n")
  }
}

cat("Available genes list:\n")
for(i in 1:length(available_genes)) {
  cat("  ", available_genes[i], "\n")
}

# =============================================================================
# PREPARE DATA FOR HEATMAP
# =============================================================================

cat("\n=== CALCULATING AVERAGE EXPRESSION ===\n")

# Join layers if using Seurat v5
if(packageVersion("Seurat") >= "5.0.0") {
  sce_EC_filtered[["RNA"]] <- JoinLayers(sce_EC_filtered[["RNA"]])
}

# Calculate average expression for each combination of EC subtype and sample group
cat("Calculating average expression by EC subtype and group...\n")
avg_exp_matrix <- AverageExpression(sce_EC_filtered,
                                   features = available_genes,
                                   group.by = c("EC_subtype", "group"),
                                   slot = "data")$RNA

cat("Expression matrix dimensions:", dim(avg_exp_matrix), "\n")
cat("Column names (EC_subtype_group combinations):\n")
print(colnames(avg_exp_matrix))

# =============================================================================
# CREATE ANNOTATIONS
# =============================================================================

cat("\n=== CREATING ANNOTATIONS ===\n")

# Extract EC subtype and group from column names
col_info <- data.frame(
  sample_name = colnames(avg_exp_matrix),
  stringsAsFactors = FALSE
)

# Parse column names (format: ECsubtype_group)
col_info$EC_subtype <- sub("_.*", "", col_info$sample_name)
col_info$group <- sub(".*_", "", col_info$sample_name)

cat("Parsed column information:\n")
print(col_info)

# Create annotation dataframe for pheatmap
annotation_col <- data.frame(
  row.names = col_info$sample_name,
  `Cell type` = col_info$EC_subtype,
  `Sample type` = col_info$group,
  check.names = FALSE
)

cat("Annotation dataframe:\n")
print(annotation_col)

# =============================================================================
# DEFINE COLORS FOR ANNOTATIONS
# =============================================================================

# Define colors for cell types (EC subtypes)
cell_type_colors <- c("ArtECs" = "#E74C3C",   # Red
                     "CapECs" = "#9B59B6",   # Purple
                     "VenECs" = "#27AE60",   # Green
                     "LECs" = "#F1C40F")     # Yellow

# Define colors for sample types
sample_type_colors <- c("CTRL" = "#2CA02C",     # Green
                       "IPF" = "#FF7F00",       # Orange
                       "SSc-ILD" = "#89CFF0")   # baby blue

# Create annotation colors list
annotation_colors <- list(
  `Cell type` = cell_type_colors[names(cell_type_colors) %in% unique(annotation_col$`Cell type`)],
  `Sample type` = sample_type_colors[names(sample_type_colors) %in% unique(annotation_col$`Sample type`)]
)

cat("Annotation colors:\n")
print(annotation_colors)

# =============================================================================
# CREATE PHEATMAP HEATMAP
# =============================================================================

cat("\n=== CREATING PHEATMAP HEATMAP ===\n")

# Create yellow to purple color palette (matching your reference image)
yellow_purple_colors <- colorRampPalette(c("#FFFF99", "#9370DB"))(100)

# Alternative: More similar to the attached image
# yellow_purple_colors <- colorRampPalette(c("#FFF2CC", "#4A148C"))(100)

cat("Creating heatmap...\n")

# Create the heatmap
heatmap_plot <- pheatmap(avg_exp_matrix,
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors,
                        color = yellow_purple_colors,
                        scale = "row",  # Scale by genes (rows)
                        cluster_rows = TRUE,   # Cluster genes
                        cluster_cols = FALSE,  # Don't cluster samples (keep grouping)
                        show_rownames = TRUE,
                        show_colnames = FALSE,  # Hide column names (show in annotation instead)
                        fontsize_row = 8,
                        fontsize_col = 8,
                        main = "EC Subtype Marker Expression (Filtered Clusters)",
                        border_color = "white",
                        cellwidth = 20,
                        cellheight = 12,
                        filename = "EC_Filtered_Pheatmap_Heatmap2.pdf",
                        width = 10,
                        height = 12)

# Display the heatmap
print(heatmap_plot)

# =============================================================================
# ALTERNATIVE VERSION WITH DIFFERENT CLUSTERING
# =============================================================================

cat("\n=== CREATING ALTERNATIVE VERSION ===\n")

# Version with no clustering (preserve gene order from your list)
heatmap_plot_ordered <- pheatmap(avg_exp_matrix[available_genes, ],  # Use original gene order
                                annotation_col = annotation_col,
                                annotation_colors = annotation_colors,
                                color = yellow_purple_colors,
                                scale = "row",
                                cluster_rows = FALSE,  # Keep original gene order
                                cluster_cols = FALSE,
                                show_rownames = TRUE,
                                show_colnames = FALSE,
                                fontsize_row = 8,
                                fontsize_col = 8,
                                main = "EC Marker Expression - Original Gene Order (Filtered)",
                                border_color = "white",
                                cellwidth = 20,
                                cellheight = 12,
                                filename = "EC_Filtered_Pheatmap_Heatmap_Ordered.pdf",
                                width = 10,
                                height = 12)

print(heatmap_plot_ordered)

# =============================================================================
# CREATE VERSION WITH ENHANCED CONTRAST
# =============================================================================

cat("\n=== CREATING HIGH CONTRAST VERSION ===\n")

# More contrasting yellow-purple palette (closer to your reference image)
high_contrast_colors <- colorRampPalette(c("#FFFFCC", "#8E24AA"))(100)

heatmap_plot_contrast <- pheatmap(avg_exp_matrix,
                                 annotation_col = annotation_col,
                                 annotation_colors = annotation_colors,
                                 color = high_contrast_colors,
                                 scale = "row",
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 show_rownames = TRUE,
                                 show_colnames = FALSE,
                                 fontsize_row = 8,
                                 main = "EC Marker Expression - High Contrast (Filtered)",
                                 border_color = "white",
                                 cellwidth = 20,
                                 cellheight = 12,
                                 filename = "EC_Filtered_Pheatmap_Heatmap_HighContrast.pdf",
                                 width = 10,
                                 height = 12)

print(heatmap_plot_contrast)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== HEATMAP SUMMARY ===\n")

cat("Filtered data used:\n")
cat("  Total cells:", ncol(sce_EC_filtered), "\n")
cat("  EC subtypes:", paste(unique(sce_EC_filtered$EC_subtype), collapse = ", "), "\n")
cat("  Sample groups:", paste(unique(sce_EC_filtered$group), collapse = ", "), "\n")

cat("Heatmap details:\n")
cat("  Genes plotted:", length(available_genes), "out of", length(heatmap_genes), "requested\n")
cat("  Sample combinations:", ncol(avg_exp_matrix), "\n")
cat("  Color scheme: Yellow (low) to Purple (high)\n")

# Show sample composition
cat("Sample composition:\n")
sample_composition <- table(sce_EC_filtered$EC_subtype, sce_EC_filtered$group)
print(sample_composition)

# =============================================================================
# SAVE ADDITIONAL DATA
# =============================================================================

cat("\n=== SAVING ADDITIONAL DATA ===\n")

# Save expression matrix
write.csv(avg_exp_matrix, "EC_filtered_heatmap_expression_matrix.csv")

# Save annotation information
write.csv(annotation_col, "EC_filtered_heatmap_annotations.csv")

# Save gene availability information
gene_availability <- data.frame(
  Gene = heatmap_genes,
  Available = heatmap_genes %in% rownames(sce_EC_filtered)
)
write.csv(gene_availability, "EC_filtered_heatmap_gene_availability.csv", row.names = FALSE)

# Save sample composition
write.csv(sample_composition, "EC_filtered_heatmap_sample_composition.csv")

cat("Files saved:\n")
cat("- EC_Filtered_Pheatmap_Heatmap.pdf (clustered)\n")
cat("- EC_Filtered_Pheatmap_Heatmap_Ordered.pdf (original gene order)\n")
cat("- EC_Filtered_Pheatmap_Heatmap_HighContrast.pdf (enhanced colors)\n")
cat("- EC_filtered_heatmap_expression_matrix.csv\n")
cat("- EC_filtered_heatmap_annotations.csv\n")
cat("- EC_filtered_heatmap_gene_availability.csv\n")
cat("- EC_filtered_heatmap_sample_composition.csv\n")

cat("\n=== PHEATMAP HEATMAP CREATION COMPLETE ===\n")
