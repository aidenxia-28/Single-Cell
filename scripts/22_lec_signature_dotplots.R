## Lymphatic branch: LEC signature dot plots
## Stage 22 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: LEC_Signature_DotPlot_*.pdf, LEC_*_cluster_*.csv
## Two versions: all eight subclusters in order, and the reordered / merged
## view (LEC0, LEC6, LEC1, LEC4, LEC2, LEC3/5/7) used in the write-up. The
## notebook had the plotting code twice; it is one function here.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(reshape2) })

sce_LECs <- readRDS("LECs_subclustered.rds")
cat("Clusters in dataset:\n"); print(table(sce_LECs$LEC_subcluster))
cat("Total cells:", ncol(sce_LECs), "\n")

if (packageVersion("Seurat") >= "5.0.0") {
  sce_LECs[["RNA"]] <- JoinLayers(sce_LECs[["RNA"]])
}

#' Scaled-expression / percent-expressing dot plot of `genes` (bottom to top)
#' across the levels of `cluster_col` in `cluster_order`.
lec_signature_dotplot <- function(obj, genes, cluster_col, cluster_order, tag, title,
                                  aspect = c(1.2, 1.5)) {
  available_genes <- intersect(genes, rownames(obj))
  cat("Genes available:", paste(available_genes, collapse = ", "), "\n")
  cat("Missing genes:", paste(setdiff(genes, available_genes), collapse = ", "), "\n")

  Idents(obj) <- obj[[cluster_col]][, 1]
  avg_exp <- AverageExpression(obj, features = available_genes, group.by = cluster_col, slot = "data")$RNA
  for (cl in cluster_order) if (!cl %in% colnames(avg_exp)) avg_exp[, cl] <- 0
  avg_exp <- avg_exp[, cluster_order, drop = FALSE]
  avg_exp_scaled <- t(scale(t(avg_exp)))

  percent_exp <- matrix(0, nrow = length(available_genes), ncol = length(cluster_order),
                        dimnames = list(available_genes, cluster_order))
  data_layer <- if (packageVersion("Seurat") >= "5.0.0") LayerData(obj, assay = "RNA", layer = "data") else
                GetAssayData(obj, slot = "data", assay = "RNA")
  for (cl in cluster_order) {
    cluster_cells <- Cells(obj)[obj[[cluster_col]][, 1] == cl]
    if (length(cluster_cells) > 0) {
      for (gene in available_genes) {
        gene_exp <- data_layer[gene, cluster_cells]
        percent_exp[gene, cl] <- sum(gene_exp > 0) / length(gene_exp) * 100
      }
    }
  }

  plot_data <- merge(
    setNames(melt(avg_exp_scaled), c("Gene", "Cluster", "Scaled_Expression")),
    setNames(melt(percent_exp),    c("Gene", "Cluster", "Percent_Expressing")),
    by = c("Gene", "Cluster"))
  plot_data$Gene    <- factor(plot_data$Gene, levels = available_genes)
  plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)
  plot_data$Size_Bin <- cut(plot_data$Percent_Expressing, breaks = c(-Inf, 20, 40, 60, 80, Inf),
                            labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"), include.lowest = TRUE)
  size_mapping <- c("0-20%" = 1, "20-40%" = 3, "40-60%" = 4, "60-80%" = 5, "80-100%" = 6)

  p <- ggplot(plot_data, aes(x = Cluster, y = Gene)) +
    geom_point(aes(color = Scaled_Expression, size = Size_Bin), alpha = 0.8) +
    scale_color_gradient2(low = "#FFF5EE", mid = "#FF8C00", high = "#8B4513", midpoint = 0,
                          name = "Scaled\nExpression") +
    scale_size_manual(values = size_mapping, name = "% Expressing",
                      guide = guide_legend(override.aes = list(color = "black"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"), legend.title = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 10, 10, "mm")) +
    labs(title = title, x = "LEC Subclusters", y = "Genes")
  p_tight_clean <- p + scale_x_discrete(expand = expansion(mult = 0.02, add = 0.35)) +
    coord_cartesian(clip = "off") + theme(aspect.ratio = aspect[1], plot.margin = margin(12, 24, 12, 24, "pt"))
  p_tighter <- p + scale_x_discrete(expand = expansion(mult = 0.01, add = 0.15)) +
    coord_cartesian(clip = "off") + theme(aspect.ratio = aspect[2], plot.margin = margin(12, 36, 12, 36, "pt"))
  print(p); print(p_tight_clean); print(p_tighter)

  ggsave(paste0("LEC_Signature_DotPlot_", tag, "_BrownOrange.pdf"), p,             width = 12, height = 10, dpi = 300)
  ggsave(paste0("LEC_Signature_DotPlot_", tag, "_TightClean.pdf"),  p_tight_clean, width = 10, height = 10, dpi = 300)
  ggsave(paste0("LEC_Signature_DotPlot_", tag, "_Tighter.pdf"),     p_tighter,     width = 9,  height = 12, dpi = 300)
  write.csv(plot_data,   paste0("LEC_signature_plot_data_", tag, ".csv"), row.names = FALSE)
  write.csv(avg_exp,     paste0("LEC_", tag, "_cluster_average_expression.csv"))
  write.csv(percent_exp, paste0("LEC_", tag, "_cluster_percent_expressing.csv"))

  cluster_group_table <- table(obj[[cluster_col]][, 1], obj$group)
  write.csv(cluster_group_table, paste0("LEC_", tag, "_cluster_composition_by_group.csv"))
  cat("\nGroup composition per cluster:\n"); print(cluster_group_table)
  props <- prop.table(cluster_group_table, margin = 2)
  for (group in colnames(props)) {
    cat("Most enriched in", group, ":", names(which.max(props[, group])), "(", round(max(props[, group]), 3), ")\n")
  }
  for (gene in available_genes) {
    cat(gene, "highest in:", colnames(avg_exp)[which.max(avg_exp[gene, ])], "\n")
  }
  invisible(plot_data)
}

# ---- 1) All eight subclusters, LEC0 ... LEC7 ----------------------------------
genes_all <- c("HIF1A",                                            # bottom
               "ICAM1", "SELP", "TNFSF4", "CCL2", "CXCL1", "CXCL2",  # hypoxia / inflammation
               "CD74",                                             # antigen presentation
               "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DRB1",
               "D2C4", "PROX1", "FLT4", "TEK", "EPAS1", "S1PR1", "EBF1", "MEOX1")  # top
lec_signature_dotplot(sce_LECs, genes_all, "LEC_subcluster", sort(unique(as.character(sce_LECs$LEC_subcluster))),
                      tag = "NewClusters", title = "LEC Signature Expression - New Subclusters (All Groups)",
                      aspect = c(1.2, 1.5))

# ---- 2) Reordered: LEC0, LEC6, LEC1, LEC4, LEC2, LEC3/5/7 merged -----------------
sce_LECs$LEC_reordered <- as.character(sce_LECs$LEC_subcluster)
sce_LECs$LEC_reordered[sce_LECs$LEC_subcluster %in% c("LEC3", "LEC5", "LEC7")] <- "LEC3/5/7"
cat("\nOriginal clusters merged into LEC3/5/7:\n")
print(table(sce_LECs$LEC_subcluster[sce_LECs$LEC_reordered == "LEC3/5/7"]))

genes_reordered <- c("HIF1A", "EPAS1", "ICAM1", "CCL2", "CXCL1", "CXCL2", "CD74",
                     "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "D2C4", "FLT4", "EBF1")
lec_signature_dotplot(sce_LECs, genes_reordered, "LEC_reordered", c("LEC0", "LEC6", "LEC1", "LEC4", "LEC2", "LEC3/5/7"),
                      tag = "Reordered_0-6-1-4-2-357",
                      title = "LEC Signature Expression - Reordered Clusters (0,6,1,4,2,3/5/7)",
                      aspect = c(1.8, 2.0))
