## Doublet flagging and EC / immune cross-contamination cleaning
## Stage 03 of the GSE159354 EC / immune pipeline.
##
## Inputs:  GSE159354_FilteredClusters_EC-Immune_harmony.rds (stage 02)
## Outputs: GSE159354_Clean_EC-Immune_harmony.rds, Cell_Quality_Diagnostics_v5.pdf
## Doublets: top 5% of z(nFeature)+z(nCount). Contamination: EC cells above the
## 90th percentile of the immune module score, and immune cells above the 90th
## percentile of the EC module score. The "Clean" object is the starting point
## of the immune branch (stage 30). Two more aggressive variants live in
## scripts/exploratory/.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

# ===== DOUBLET DETECTION AND CELL TYPE PURIFICATION =====
filtered_lung_seurat <- readRDS("GSE159354_FilteredClusters_EC-Immune_harmony.rds")

# Install and load doublet detection packages if needed
# BiocManager::install("scDblFinder")

# ===== SEURAT V5 COMPATIBLE DOUBLET AND QUALITY CONTROL =====

# Alternative approach without SingleCellExperiment conversion
# This uses Seurat-native functions that work better with v5

# ===== 1) Calculate quality metrics =====
# Add mitochondrial gene percentage (if not already present)
if (!"percent.mt" %in% colnames(filtered_lung_seurat@meta.data)) {
  filtered_lung_seurat[["percent.mt"]] <- PercentageFeatureSet(filtered_lung_seurat, pattern = "^MT-")
}

# Add ribosomal gene percentage
filtered_lung_seurat[["percent.ribo"]] <- PercentageFeatureSet(filtered_lung_seurat, pattern = "^RP[SL]")

# Add hemoglobin gene percentage (potential doublet indicator)
filtered_lung_seurat[["percent.hb"]] <- PercentageFeatureSet(filtered_lung_seurat, pattern = "^HB[^(P)]")

# ===== 2) Simple doublet detection based on quality metrics =====
# High gene count + high UMI count often indicates doublets
# Calculate doublet score based on multiple metrics

# Z-score normalize the metrics
filtered_lung_seurat$nFeature_RNA_z <- scale(filtered_lung_seurat$nFeature_RNA)[,1]
filtered_lung_seurat$nCount_RNA_z <- scale(filtered_lung_seurat$nCount_RNA)[,1]

# Simple doublet score (high gene count + high UMI count)
filtered_lung_seurat$doublet_score_simple <- filtered_lung_seurat$nFeature_RNA_z + filtered_lung_seurat$nCount_RNA_z

# Flag potential doublets (top 5% of doublet scores)
doublet_threshold <- quantile(filtered_lung_seurat$doublet_score_simple, 0.95)
filtered_lung_seurat$is_doublet <- filtered_lung_seurat$doublet_score_simple > doublet_threshold

# ===== 3) Calculate gene signature scores =====
# Immune genes that are present in the dataset
immune_genes_all <- c("PTPRC","CD3D","CD8A","CD4","NKG7","MS4A1","MRC1","CD163",
                     "LYZ","FCGR3A","MS4A7","HLA-DRA","HLA-DRB1","HLA-DPA1",
                     "HLA-DPB1","CD68")
immune_genes_available <- immune_genes_all[immune_genes_all %in% rownames(filtered_lung_seurat)]

# EC genes that are present in the dataset
ec_genes_all <- c("PECAM1","CDH5","VWF","KDR","CLDN5","EMCN","ESAM")
ec_genes_available <- ec_genes_all[ec_genes_all %in% rownames(filtered_lung_seurat)]

cat("Immune genes available:", length(immune_genes_available), "out of", length(immune_genes_all), "\n")
cat("EC genes available:", length(ec_genes_available), "out of", length(ec_genes_all), "\n")

# Calculate module scores using Seurat v5 compatible method
filtered_lung_seurat <- AddModuleScore(filtered_lung_seurat,
                                      features = list(immune_genes_available),
                                      name = "Immune_Score")

filtered_lung_seurat <- AddModuleScore(filtered_lung_seurat,
                                      features = list(ec_genes_available),
                                      name = "EC_Score")

# ===== 4) Identify problematic cells =====
# Get cell identifiers
ec_cells <- WhichCells(filtered_lung_seurat, expression = cellset == "EC")
immune_cells <- WhichCells(filtered_lung_seurat, expression = cellset == "Immune")

# Calculate thresholds for contamination
# Use 90th percentile to be more stringent
pure_ec_immune_threshold <- quantile(filtered_lung_seurat$Immune_Score1[ec_cells], 0.90, na.rm = TRUE)
pure_immune_ec_threshold <- quantile(filtered_lung_seurat$EC_Score1[immune_cells], 0.90, na.rm = TRUE)

# Flag problematic cells
filtered_lung_seurat$problematic_EC <- (filtered_lung_seurat$cellset == "EC" &
                                       filtered_lung_seurat$Immune_Score1 > pure_ec_immune_threshold)

filtered_lung_seurat$problematic_Immune <- (filtered_lung_seurat$cellset == "Immune" &
                                           filtered_lung_seurat$EC_Score1 > pure_immune_ec_threshold)

# ===== 5) Summary statistics =====
cat("\n=== CELL QUALITY SUMMARY ===\n")
cat("Total cells:", ncol(filtered_lung_seurat), "\n")
cat("Predicted doublets (top 5% by complexity):", sum(filtered_lung_seurat$is_doublet, na.rm = TRUE), "\n")
cat("Problematic EC cells (high immune score):", sum(filtered_lung_seurat$problematic_EC, na.rm = TRUE), "\n")
cat("Problematic Immune cells (high EC score):", sum(filtered_lung_seurat$problematic_Immune, na.rm = TRUE), "\n")

# Print score thresholds
cat("EC immune contamination threshold:", round(pure_ec_immune_threshold, 3), "\n")
cat("Immune EC contamination threshold:", round(pure_immune_ec_threshold, 3), "\n")

# ===== 6) Create diagnostic plots =====
# Plot doublet scores
doublet_plot <- FeaturePlot(filtered_lung_seurat, features = "doublet_score_simple",
                           reduction = "umap") +
                           ggtitle("Simple Doublet Score (nGenes + nUMIs)")

# Plot quality metrics
nfeature_plot <- FeaturePlot(filtered_lung_seurat, features = "nFeature_RNA",
                            reduction = "umap") +
                            ggtitle("Number of Genes")

ncount_plot <- FeaturePlot(filtered_lung_seurat, features = "nCount_RNA",
                          reduction = "umap") +
                          ggtitle("Number of UMIs")

# Plot module scores
immune_score_plot <- FeaturePlot(filtered_lung_seurat, features = "Immune_Score1",
                                reduction = "umap") +
                                ggtitle("Immune Gene Signature Score")

ec_score_plot <- FeaturePlot(filtered_lung_seurat, features = "EC_Score1",
                            reduction = "umap") +
                            ggtitle("EC Gene Signature Score")

# Create quality classification
filtered_lung_seurat$cell_quality <- "Good"
filtered_lung_seurat$cell_quality[filtered_lung_seurat$is_doublet] <- "Doublet"
filtered_lung_seurat$cell_quality[filtered_lung_seurat$problematic_EC] <- "Contaminated_EC"
filtered_lung_seurat$cell_quality[filtered_lung_seurat$problematic_Immune] <- "Contaminated_Immune"

filtered_lung_seurat$cell_quality <- factor(filtered_lung_seurat$cell_quality,
                                           levels = c("Good", "Doublet", "Contaminated_EC", "Contaminated_Immune"))

quality_plot <- DimPlot(filtered_lung_seurat, group.by = "cell_quality",
                       reduction = "umap",
                       cols = c("Good" = "lightgray", "Doublet" = "black",
                               "Contaminated_EC" = "#FFA500", "Contaminated_Immune" = "#9370DB")) +
                       ggtitle("Cell Quality Assessment")

print(doublet_plot)
print(nfeature_plot)
print(ncount_plot)
print(immune_score_plot)
print(ec_score_plot)
print(quality_plot)

# ===== 7) Create clean object =====
# Remove problematic cells
cells_to_remove <- filtered_lung_seurat$is_doublet |
                   filtered_lung_seurat$problematic_EC |
                   filtered_lung_seurat$problematic_Immune

cells_to_keep <- !cells_to_remove
clean_lung_seurat <- subset(filtered_lung_seurat, cells = colnames(filtered_lung_seurat)[cells_to_keep])

cat("\n=== CLEANING RESULTS ===\n")
cat("Original filtered cells:", ncol(filtered_lung_seurat), "\n")
cat("Clean cells:", ncol(clean_lung_seurat), "\n")
cat("Cells removed:", sum(cells_to_remove, na.rm = TRUE), "\n")
cat("Percentage retained:", round((ncol(clean_lung_seurat)/ncol(filtered_lung_seurat))*100, 1), "%\n")

# ===== 8) Plot clean results =====
clean_celltype_plot <- DimPlot(clean_lung_seurat, group.by = "cellset", reduction = "umap",
                               cols = c("EC" = "#E31A1C", "Immune" = "#1F78B4")) +
                               ggtitle("Clean UMAP: Endothelial vs Immune Cells") +
                               theme_minimal()

print(clean_celltype_plot)

# ===== 9) Save everything =====
saveRDS(clean_lung_seurat, file = "GSE159354_Clean_EC-Immune_harmony.rds")

# Save all diagnostic plots
pdf("Cell_Quality_Diagnostics_v5.pdf", height = 12, width = 16)
print(doublet_plot)
print(nfeature_plot)
print(ncount_plot)
print(immune_score_plot)
print(ec_score_plot)
print(quality_plot)
print(clean_celltype_plot)
dev.off()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Clean object saved: GSE159354_Clean_EC-Immune_harmony.rds\n")
cat("Diagnostic plots saved: Cell_Quality_Diagnostics_v5.pdf\n")
cat("Use 'clean_lung_seurat' for your heatmap analysis!\n")
