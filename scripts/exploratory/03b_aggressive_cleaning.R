## More aggressive EC / immune contamination cleaning (exploratory)
## Exploratory variant of stage 03; not used by the main branches.
##
## Inputs:  GSE159354_Clean_EC-Immune_harmony.rds (stage 03)
## Outputs: GSE159354_UltraClean_EC-Immune_harmony.rds, Ultra_Clean_Diagnostics.pdf, UltraClean_Expression_Heatmap.pdf
## 95th / 99th percentile thresholds, a macrophage-specific score, and an
## EC:immune gene-ratio check.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

clean_lung_seurat <- readRDS("GSE159354_Clean_EC-Immune_harmony.rds")

# ===== MORE AGGRESSIVE CLEANING STRATEGIES =====

# Strategy 1: Much more stringent thresholds (99th percentile instead of 90th)
# Strategy 2: Specific macrophage marker contamination detection
# Strategy 3: Multi-step filtering approach

# ===== 1) Identify remaining problematic cells more aggressively =====

# Calculate very specific macrophage signature (the main culprits)
macrophage_genes <- c("MRC1", "CD163", "CD68", "LYZ", "FCGR3A", "MS4A7")
macrophage_genes_available <- macrophage_genes[macrophage_genes %in% rownames(clean_lung_seurat)]

# Add specific macrophage score
clean_lung_seurat <- AddModuleScore(clean_lung_seurat,
                                   features = list(macrophage_genes_available),
                                   name = "Macrophage_Score")

# ===== 2) More stringent contamination detection =====
ec_cells <- WhichCells(clean_lung_seurat, expression = cellset == "EC")
immune_cells <- WhichCells(clean_lung_seurat, expression = cellset == "Immune")

# Use 95th percentile for general immune score, 99th for macrophage score
stringent_immune_threshold <- quantile(clean_lung_seurat$Immune_Score1[ec_cells], 0.95, na.rm = TRUE)
stringent_macro_threshold <- quantile(clean_lung_seurat$Macrophage_Score1[ec_cells], 0.99, na.rm = TRUE)

# Also check for high individual gene expression of key contaminants
individual_gene_contamination <- list()

for (gene in c("MRC1", "CD163", "CD68")) {
  if (gene %in% rownames(clean_lung_seurat)) {
    gene_expr <- FetchData(clean_lung_seurat, vars = gene, cells = ec_cells)
    threshold_99 <- quantile(gene_expr[[gene]], 0.99, na.rm = TRUE)
    individual_gene_contamination[[gene]] <- ec_cells[gene_expr[[gene]] > threshold_99]
    cat("EC cells with high", gene, "expression:", length(individual_gene_contamination[[gene]]), "\n")
  }
}

# ===== 3) Flag highly contaminated EC cells =====
highly_contaminated_EC <- unique(c(
  # High general immune score
  ec_cells[clean_lung_seurat$Immune_Score1[ec_cells] > stringent_immune_threshold],
  # High macrophage score
  ec_cells[clean_lung_seurat$Macrophage_Score1[ec_cells] > stringent_macro_threshold],
  # High individual gene expression
  unlist(individual_gene_contamination)
))

# Also flag immune cells with high EC scores (more stringent)
stringent_ec_threshold <- quantile(clean_lung_seurat$EC_Score1[immune_cells], 0.95, na.rm = TRUE)
highly_contaminated_Immune <- immune_cells[clean_lung_seurat$EC_Score1[immune_cells] > stringent_ec_threshold]

cat("\n=== STRINGENT CONTAMINATION DETECTION ===\n")
cat("Highly contaminated EC cells:", length(highly_contaminated_EC), "\n")
cat("Highly contaminated Immune cells:", length(highly_contaminated_Immune), "\n")

# ===== 4) Visualize the most problematic cells =====
# Mark the cells we're about to remove
clean_lung_seurat$stringent_quality <- "Clean"
clean_lung_seurat$stringent_quality[highly_contaminated_EC] <- "Contaminated_EC"
clean_lung_seurat$stringent_quality[highly_contaminated_Immune] <- "Contaminated_Immune"

stringent_quality_plot <- DimPlot(clean_lung_seurat, group.by = "stringent_quality",
                                 reduction = "umap",
                                 cols = c("Clean" = "lightgray",
                                         "Contaminated_EC" = "#FFA500",
                                         "Contaminated_Immune" = "#9370DB")) +
                                 ggtitle("Stringent Quality Assessment - Cells to Remove")

print(stringent_quality_plot)

# Plot macrophage scores
macro_score_plot <- FeaturePlot(clean_lung_seurat, features = "Macrophage_Score1",
                               reduction = "umap") +
                               ggtitle("Macrophage Gene Signature Score")

print(macro_score_plot)

# ===== 5) Create ultra-clean object =====
cells_to_remove_stringent <- c(highly_contaminated_EC, highly_contaminated_Immune)
cells_to_keep_stringent <- !colnames(clean_lung_seurat) %in% cells_to_remove_stringent

ultra_clean_seurat <- subset(clean_lung_seurat, cells = colnames(clean_lung_seurat)[cells_to_keep_stringent])

cat("\n=== ULTRA-CLEAN RESULTS ===\n")
cat("Clean object cells:", ncol(clean_lung_seurat), "\n")
cat("Ultra-clean cells:", ncol(ultra_clean_seurat), "\n")
cat("Additional cells removed:", length(cells_to_remove_stringent), "\n")
cat("Total retention from original:", round((ncol(ultra_clean_seurat)/ncol(clean_lung_seurat))*100, 1), "%\n")

# ===== 6) Alternative approach: Manual gene ratio filtering =====
# Sometimes the cleanest approach is to look at gene expression ratios

# Calculate EC:Immune gene ratio for each cell
ec_core_genes <- c("PECAM1", "CDH5", "VWF", "KDR")[c("PECAM1", "CDH5", "VWF", "KDR") %in% rownames(ultra_clean_seurat)]
immune_core_genes <- c("PTPRC", "CD68", "MRC1")[c("PTPRC", "CD68", "MRC1") %in% rownames(ultra_clean_seurat)]

cat("EC core genes found:", paste(ec_core_genes, collapse = ", "), "\n")
cat("Immune core genes found:", paste(immune_core_genes, collapse = ", "), "\n")

if (length(ec_core_genes) > 0 & length(immune_core_genes) > 0) {

  # Calculate mean expression for each cell - fixed approach
  ec_data <- FetchData(ultra_clean_seurat, vars = ec_core_genes, layer = "data")
  immune_data <- FetchData(ultra_clean_seurat, vars = immune_core_genes, layer = "data")

  # Calculate means properly
  if (length(ec_core_genes) == 1) {
    ec_mean_expr <- ec_data[[1]]
  } else {
    ec_mean_expr <- rowMeans(ec_data)
  }

  if (length(immune_core_genes) == 1) {
    immune_mean_expr <- immune_data[[1]]
  } else {
    immune_mean_expr <- rowMeans(immune_data)
  }

  # Ensure same length and same cell order
  cell_names <- intersect(names(ec_mean_expr), names(immune_mean_expr))
  ec_mean_expr <- ec_mean_expr[cell_names]
  immune_mean_expr <- immune_mean_expr[cell_names]

  # Calculate ratio (EC genes / Immune genes)
  gene_ratio <- log2((ec_mean_expr + 0.01) / (immune_mean_expr + 0.01))
  names(gene_ratio) <- cell_names

  # Add to metadata properly
  ultra_clean_seurat$EC_to_Immune_ratio <- 0  # Initialize
  ultra_clean_seurat$EC_to_Immune_ratio[cell_names] <- gene_ratio

  # Flag cells with wrong ratios for their assigned type
  ec_cells_ultra <- WhichCells(ultra_clean_seurat, expression = cellset == "EC")
  immune_cells_ultra <- WhichCells(ultra_clean_seurat, expression = cellset == "Immune")

  # EC cells should have positive ratios, immune cells negative ratios
  misassigned_ec <- ec_cells_ultra[ultra_clean_seurat$EC_to_Immune_ratio[ec_cells_ultra] < 0]
  misassigned_immune <- immune_cells_ultra[ultra_clean_seurat$EC_to_Immune_ratio[immune_cells_ultra] > 0]

  cat("\n=== GENE RATIO ANALYSIS ===\n")
  cat("EC cells with negative EC:Immune ratio:", length(misassigned_ec), "\n")
  cat("Immune cells with positive EC:Immune ratio:", length(misassigned_immune), "\n")

  # Plot the ratios
  ratio_plot <- FeaturePlot(ultra_clean_seurat, features = "EC_to_Immune_ratio",
                           reduction = "umap") +
                           ggtitle("EC:Immune Gene Expression Ratio\n(Red=More EC genes, Blue=More Immune genes)") +
                           scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

  print(ratio_plot)

  # Option to remove misassigned cells
  if (length(misassigned_ec) + length(misassigned_immune) < ncol(ultra_clean_seurat) * 0.05) {  # Less than 5%
    purest_seurat <- subset(ultra_clean_seurat,
                           cells = colnames(ultra_clean_seurat)[!colnames(ultra_clean_seurat) %in% c(misassigned_ec, misassigned_immune)])

    cat("Created 'purest_seurat' object with", ncol(purest_seurat), "cells\n")
  } else {
    cat("Too many misassigned cells (", length(misassigned_ec) + length(misassigned_immune),
        "out of", ncol(ultra_clean_seurat), ") - not creating purest object\n")
  }
} else {
  cat("Not enough core genes found for ratio analysis\n")
}

# ===== 7) Final ultra-clean UMAP =====
ultra_clean_plot <- DimPlot(ultra_clean_seurat, group.by = "cellset", reduction = "umap",
                           cols = c("EC" = "#E31A1C", "Immune" = "#1F78B4")) +
                           ggtitle("Ultra-Clean UMAP: Endothelial vs Immune Cells") +
                           theme_minimal()

print(ultra_clean_plot)

# ===== 8) Save results =====
saveRDS(ultra_clean_seurat, file = "GSE159354_UltraClean_EC-Immune_harmony.rds")

# Save diagnostic plots
pdf("Ultra_Clean_Diagnostics.pdf", height = 12, width = 16)
print(stringent_quality_plot)
print(macro_score_plot)
if (exists("ratio_plot")) print(ratio_plot)
print(ultra_clean_plot)
dev.off()

cat("\n=== NEXT STEP ===\n")
cat("Run the heatmap code again with 'ultra_clean_seurat' object\n")
cat("This should show much cleaner separation between EC and immune signatures\n")

# Marker-panel heatmap on the ultra-clean object (see R/marker_panel_heatmap.R)
plot_marker_panel_heatmap(ultra_clean_seurat, label = "UltraClean")
report_cleaning_delta(clean_lung_seurat, ultra_clean_seurat, "CLEAN", "ULTRA-CLEAN")
