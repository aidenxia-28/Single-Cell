## Immune branch: re-cluster the cleaned CD45+ cells
## Stage 30 of the GSE159354 EC / immune pipeline.
##
## Inputs:  GSE159354_Clean_EC-Immune_harmony.rds (stage 03)
## Outputs: Immune_only_seurat_processed.rds, Immune_UMAP_by_group.pdf, Immune_UMAP_clusters.pdf
## Immune cells are re-normalised (3000 HVGs, 50 PCs, Harmony on donor / group
## / fraction / gsm where they vary) and clustered at resolutions 0.3-1.2;
## resolution 0.8 is used downstream.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ── 1) Load the clean EC+Immune object ─────────────────────────────────────────
# If the object is already in memory, this keeps it. Otherwise it loads from disk.
if (!exists("clean_lung_seurat")) {
  obj_path <- "GSE159354_Clean_EC-Immune_harmony.rds"
  if (!file.exists(obj_path)) stop("Cannot find '", obj_path, "' in the working directory.")
  clean_lung_seurat <- readRDS(obj_path)
}

# ── 2) Subset to IMMUNE only ───────────────────────────────────────────────────
immune_seurat <- subset(clean_lung_seurat, subset = cellset == "Immune")
stopifnot(ncol(immune_seurat) > 0)

# Seurat v5 layered assay? Join layers before PCA/Harmony
if (packageVersion("Seurat") >= "5.0.0" && inherits(immune_seurat[["RNA"]], "Assay5")) {
  immune_seurat[["RNA"]] <- JoinLayers(immune_seurat[["RNA"]])
}

# ── 3) Prepare metadata; keep only covariates that actually vary ───────────────
meta_candidates <- intersect(c("donor","group","fraction","gsm"), colnames(immune_seurat@meta.data))

# Make them clean factors (replace NA with "Unknown", drop unused levels)
for (v in meta_candidates) {
  x <- as.character(immune_seurat[[v]][,1])
  x[is.na(x) | x == ""] <- "Unknown"
  immune_seurat[[v]] <- droplevels(factor(x))
}

# Count levels and pick those with ≥ 2 levels
lev_counts <- sapply(meta_candidates, function(v) nlevels(immune_seurat[[v]]))
vars_for_harmony <- names(lev_counts[lev_counts >= 2])

cat("\n=== HARMONY COVARIATE CHECK (Immune subset) ===\n")
print(lev_counts)
cat("Using for Harmony:", ifelse(length(vars_for_harmony) > 0, paste(vars_for_harmony, collapse = ", "), "NONE (skip Harmony)"), "\n\n")

# ── 4) Standard preprocessing ──────────────────────────────────────────────────
DefaultAssay(immune_seurat) <- "RNA"
immune_seurat <- NormalizeData(immune_seurat, normalization.method = "LogNormalize",
                               scale.factor = 1e4, verbose = FALSE)
immune_seurat <- FindVariableFeatures(immune_seurat, nfeatures = 3000,
                                      selection.method = "vst", verbose = FALSE)
immune_seurat <- ScaleData(immune_seurat, verbose = FALSE)
immune_seurat <- RunPCA(immune_seurat, npcs = 50, verbose = FALSE)

# ── 5) Harmony if possible; otherwise skip ─────────────────────────────────────
if (length(vars_for_harmony) >= 1) {
  message("Running Harmony over: ", paste(vars_for_harmony, collapse = ", "))
  # Note: we rely on default dims for compatibility across Seurat/Harmony versions
  immune_seurat <- RunHarmony(immune_seurat, group.by.vars = vars_for_harmony, verbose = TRUE)
  immune_seurat <- RunUMAP(immune_seurat, reduction = "harmony", dims = 1:30)
  immune_seurat <- FindNeighbors(immune_seurat, reduction = "harmony", dims = 1:30)
} else {
  message("No covariate has ≥2 levels in the Immune subset; skipping Harmony.")
  immune_seurat <- RunUMAP(immune_seurat, reduction = "pca", dims = 1:30)
  immune_seurat <- FindNeighbors(immune_seurat, reduction = "pca", dims = 1:30)
}

immune_seurat <- FindClusters(immune_seurat, resolution = 0.6)

# ── 6) Quick visuals & save ────────────────────────────────────────────────────
p_by_group   <- DimPlot(immune_seurat, reduction = "umap", group.by = "group") +
  ggtitle("Immune cells: UMAP by disease group")

p_clusters   <- DimPlot(immune_seurat, reduction = "umap", label = TRUE) +
  ggtitle("Immune cells: UMAP by clusters")

print(p_by_group)
print(p_clusters)

saveRDS(immune_seurat, file = "Immune_only_seurat_processed.rds")
ggsave("Immune_UMAP_by_group.pdf",   p_by_group, width = 8, height = 6)
ggsave("Immune_UMAP_clusters.pdf",   p_clusters, width = 8, height = 6)

cat("\n=== DONE ===\nSaved object: Immune_only_seurat_processed.rds\n",
    "Plots: Immune_UMAP_by_group.pdf, Immune_UMAP_clusters.pdf\n", sep = "")


resolutions <- c(0.3, 0.5, 0.8, 1.0, 1.2)
for (res in resolutions) {
  immune_seurat <- FindClusters(immune_seurat,
                               resolution = res,
                               verbose = FALSE)
}

# Set default clustering (resolution 0.8 is often good for immune subclustering)
Idents(immune_seurat) <- immune_seurat$RNA_snn_res.0.8

# UMAP colored by clusters (default resolution)
immune_cluster_plot <- DimPlot(immune_seurat,
                              label = TRUE,
                              reduction = "umap",
                              label.size = 4) +
                      ggtitle("Immune Cell Subclusters (res 0.8)") +
                      theme_minimal() +
                      theme(legend.position = "none")  # Remove legend since we have labels

# UMAP colored by disease group
immune_group_plot <- DimPlot(immune_seurat,
                            group.by = "group",
                            reduction = "umap") +
                    ggtitle("Immune Cells: Disease Groups") +
                    theme_minimal()

# UMAP colored by donor
immune_donor_plot <- DimPlot(immune_seurat,
                            group.by = "donor",
                            reduction = "umap") +
                    ggtitle("Immune Cells: Donors") +
                    theme_minimal()

# Create disease status for CTRL vs Disease comparison
immune_seurat$disease_status <- ifelse(immune_seurat$group == "CTRL", "CTRL", "Disease")
immune_seurat$disease_status <- factor(immune_seurat$disease_status, levels = c("CTRL", "Disease"))

# CTRL vs Disease plot
immune_ctrl_disease_plot <- DimPlot(immune_seurat,
                                   group.by = "disease_status",
                                   reduction = "umap",
                                   cols = c("CTRL" = "#2E8B57", "Disease" = "#DC143C")) +
                           ggtitle("Immune Cells: CTRL vs Disease") +
                           theme_minimal()

print(immune_cluster_plot)
print(immune_group_plot)
print(immune_donor_plot)
print(immune_ctrl_disease_plot)

# ===== 7) Resolution comparison plots =====
cat("Creating resolution comparison plots...\n")

res_plots <- list()
for (res in resolutions) {
  cluster_col <- paste0("RNA_snn_res.", res)
  res_plots[[as.character(res)]] <- DimPlot(immune_seurat,
                                           group.by = cluster_col,
                                           label = TRUE,
                                           reduction = "umap") +
                                   ggtitle(paste("Resolution", res)) +
                                   theme_minimal() +
                                   theme(legend.position = "none")
}

# Print all resolution plots
for (plot in res_plots) {
  print(plot)
}

# ===== 8) Summary statistics =====
cat("\n=== IMMUNE SUBCLUSTER SUMMARY ===\n")

# Cluster summary at default resolution
cluster_summary <- immune_seurat@meta.data %>%
  group_by(RNA_snn_res.0.8, group) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = n_cells, values_fill = 0) %>%
  mutate(Total = CTRL + IPF + `SSc-ILD`) %>%
  arrange(as.numeric(RNA_snn_res.0.8))

cat("Cluster distribution by disease group:\n")
print(cluster_summary)

# Summary by resolution
cat("\n=== CLUSTER COUNTS BY RESOLUTION ===\n")
for (res in resolutions) {
  cluster_col <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(immune_seurat@meta.data[[cluster_col]]))
  cat("Resolution", res, ":", n_clusters, "clusters\n")
}

saveRDS(immune_seurat, file = "Immune_only_seurat_processed.rds")
