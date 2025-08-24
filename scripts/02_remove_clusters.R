## Remove non-EC / non-immune clusters from the integrated object
## Stage 02 of the GSE159354 EC / immune pipeline.
##
## Inputs:  GSE159354_ThreeGroup_EC-Immune_harmony.rds (stage 01)
## Outputs: GSE159354_FilteredClusters_EC-Immune_harmony.rds, Filtered_UMAP_Analysis.pdf
## Clusters 7, 12, 15, 16, 17 and 19 of the joint EC+immune clustering were
## judged to be contaminating / low-quality populations on the marker plots and
## are dropped here.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(tidyr) })

combined_seurat <- readRDS("GSE159354_ThreeGroup_EC-Immune_harmony.rds")

## ===== Remove specific clusters and create new object =====
# Clusters to remove: 17, 12, 16, 15, 7, and 19
clusters_to_remove <- c("17", "12", "16", "15", "7", "19")

# Create a new object by subsetting - remove unwanted clusters
filtered_lung_seurat <- subset(combined_seurat,
                              subset = seurat_clusters %in% setdiff(levels(combined_seurat$seurat_clusters),
                                                                   clusters_to_remove))

# Check how many cells remain
cat("Original object cells:", ncol(combined_seurat), "\n")
cat("Filtered object cells:", ncol(filtered_lung_seurat), "\n")
cat("Cells removed:", ncol(combined_seurat) - ncol(filtered_lung_seurat), "\n")

## ===== UMAP: Endothelial Cells vs Immune Cells (filtered object) =====
filtered_celltype_plot <- DimPlot(filtered_lung_seurat, group.by = "cellset", reduction = "umap",
                                  cols = c("EC" = "#E31A1C", "Immune" = "#1F78B4")) +
                                  ggtitle("Filtered UMAP: Endothelial Cells vs Immune Cells") +
                                  theme_minimal() +
                                  theme(legend.title = element_text(size = 12),
                                        panel.grid = element_blank())

print(filtered_celltype_plot)

## ===== Individual Group UMAPs =====
# Get cell identifiers for each group
ctrl_cells <- WhichCells(filtered_lung_seurat, expression = group == "CTRL")
ipf_cells <- WhichCells(filtered_lung_seurat, expression = group == "IPF")
ssc_cells <- WhichCells(filtered_lung_seurat, expression = group == "SSc-ILD")

# Extract UMAP coordinates for each group
ctrl_coords <- filtered_lung_seurat@reductions$umap@cell.embeddings[ctrl_cells, ]
ipf_coords <- filtered_lung_seurat@reductions$umap@cell.embeddings[ipf_cells, ]
ssc_coords <- filtered_lung_seurat@reductions$umap@cell.embeddings[ssc_cells, ]

# Check what the columns are actually called
print("UMAP coordinate column names:")
print(colnames(ctrl_coords))

# Create individual group plots
ctrl_plot <- ggplot(data.frame(ctrl_coords), aes(x = umap_1, y = umap_2)) +
  geom_point(color = "#2E8B57", size = 0.1, alpha = 0.7) +
  ggtitle("UMAP: CTRL - All Cells (Filtered)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

ipf_plot <- ggplot(data.frame(ipf_coords), aes(x = umap_1, y = umap_2)) +
  geom_point(color = "#9370DB", size = 0.1, alpha = 0.7) +
  ggtitle("UMAP: IPF - All Cells (Filtered)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

ssc_plot <- ggplot(data.frame(ssc_coords), aes(x = umap_1, y = umap_2)) +
  geom_point(color = "#DC143C", size = 0.1, alpha = 0.7) +
  ggtitle("UMAP: SSc-ILD - All Cells (Filtered)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(ctrl_plot)
print(ipf_plot)
print(ssc_plot)

## ===== CTRL vs Disease (IPF + SSc-ILD) UMAP =====
# Create disease status column
filtered_lung_seurat$disease_status <- ifelse(filtered_lung_seurat$group == "CTRL", "CTRL", "Disease")
filtered_lung_seurat$disease_status <- factor(filtered_lung_seurat$disease_status, levels = c("CTRL", "Disease"))

# Create CTRL vs Disease plot
ctrl_vs_disease_plot <- DimPlot(filtered_lung_seurat,
                               group.by = "disease_status",
                               reduction = "umap",
                               pt.size = 0.01,
                               cols = c("CTRL" = "#2E8B57", "Disease" = "#DC143C")) +
                               ggtitle("UMAP: CTRL vs Disease (IPF + SSc-ILD) - Filtered") +
                               theme_minimal() +
                               theme(panel.grid = element_blank(),
                                     plot.title = element_text(hjust = 0.5))

print(ctrl_vs_disease_plot)

## ===== Summary Statistics for Filtered Object =====
cat("\n=== FILTERED OBJECT SUMMARY ===\n")
cat("Clusters removed:", paste(clusters_to_remove, collapse = ", "), "\n")
cat("Remaining clusters:", paste(sort(as.numeric(unique(filtered_lung_seurat$seurat_clusters))), collapse = ", "), "\n")

# Create summary table for filtered object
filtered_cell_summary <- filtered_lung_seurat@meta.data %>%
  group_by(group, cellset) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  pivot_wider(names_from = cellset, values_from = n_cells, values_fill = 0) %>%
  mutate(Total = EC + Immune)

cat("\n=== FILTERED CELL COUNT SUMMARY ===\n")
print(filtered_cell_summary)

# Comparison of original vs filtered
cat("\n=== BEFORE vs AFTER FILTERING ===\n")
original_total <- ncol(combined_seurat)
filtered_total <- ncol(filtered_lung_seurat)
cat("Original total cells:", original_total, "\n")
cat("Filtered total cells:", filtered_total, "\n")
cat("Percentage retained:", round((filtered_total/original_total)*100, 1), "%\n")

## ===== Save the filtered object =====
saveRDS(filtered_lung_seurat, file = "GSE159354_FilteredClusters_EC-Immune_harmony.rds")

## ===== Create a combined plot showing all 4 UMAPs =====
library(gridExtra)

# Arrange all plots in a 2x2 grid
combined_plot <- grid.arrange(
  filtered_celltype_plot,
  ctrl_vs_disease_plot,
  arrangeGrob(ctrl_plot, ipf_plot, ssc_plot, ncol = 3),
  ncol = 1,
  heights = c(1, 1, 1)
)

# Save all plots as PDF
pdf("Filtered_UMAP_Analysis.pdf", height = 15, width = 20)
print(filtered_celltype_plot)
print(ctrl_vs_disease_plot)
print(ctrl_plot)
print(ipf_plot)
print(ssc_plot)
# Print a summary page with all plots
grid.arrange(
  filtered_celltype_plot,
  ctrl_vs_disease_plot,
  arrangeGrob(ctrl_plot, ipf_plot, ssc_plot, ncol = 3),
  ncol = 1,
  heights = c(1, 1, 1),
  top = "Filtered Lung scRNA-seq Analysis Summary"
)
dev.off()

cat("\n=== Analysis Complete ===\n")
cat("New filtered object saved as: GSE159354_FilteredClusters_EC-Immune_harmony.rds\n")
cat("All plots saved as: Filtered_UMAP_Analysis.pdf\n")
