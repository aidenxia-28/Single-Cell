## CellChat ligand-receptor analysis between LEC subclusters and the renumbered
## immune clusters. One function replaces the "complete" (LEC0-LEC7) and
## "targeted" (LEC0/3/5) notebook blocks; the directional heatmap variants that
## followed them are collected in plot_lec_to_immune_heatmaps().
##
## Raw counts of the two objects are combined on their shared genes; groups are
## labelled "LEC_<subcluster>" and "Immune_<cluster>"; CellChatDB.human;
## communications with fewer than 10 cells in a group are filtered out.

## Packages are attached when the first function is called, so sourcing this file
## does not require CellChat to be installed.
.load_cellchat_packages <- function() {
  suppressPackageStartupMessages({
    library(CellChat); library(Seurat); library(dplyr); library(pheatmap)
  })
}

#' @param sce_LECs       LEC Seurat object with Idents() = LEC_subcluster
#' @param immune         immune Seurat object with Idents() = renamed_clusters
#' @param lec_clusters   which LEC subclusters to include (NULL = all)
#' @param prefix         file prefix, e.g. "Complete" or "Targeted"
#' @param n_pathways     how many top pathways to draw as circle plots
#' @param top_n          how many top interactions to keep in the Top_* tables
run_lec_immune_cellchat <- function(sce_LECs, immune, lec_clusters = NULL, prefix = "Complete",
                                    n_pathways = 8, top_n = 25) {
  .load_cellchat_packages()
  if (!is.null(lec_clusters)) sce_LECs <- subset(sce_LECs, idents = lec_clusters)
  cat("LEC clusters included:", paste(levels(droplevels(Idents(sce_LECs))), collapse = ", "), "\n")
  cat("LEC cells:", ncol(sce_LECs), " | Immune cells:", ncol(immune), "\n")

  ## combined raw-count matrix on shared genes
  lec_raw    <- as.matrix(sce_LECs[["RNA"]]$counts)
  immune_raw <- as.matrix(immune[["RNA"]]$counts)
  common_genes <- intersect(rownames(lec_raw), rownames(immune_raw))
  cat("Common genes:", length(common_genes), "\n")
  combined_raw <- cbind(lec_raw[common_genes, ], immune_raw[common_genes, ])
  meta_data <- data.frame(
    cell_types = c(paste0("LEC_", Idents(sce_LECs)), paste0("Immune_", Idents(immune))),
    row.names  = colnames(combined_raw)
  )
  stopifnot(nrow(meta_data) == ncol(combined_raw))

  ## CellChat
  cellchat <- createCellChat(object = combined_raw, meta = meta_data, group.by = "cell_types")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  lec_ids    <- grep("LEC_",    levels(cellchat@idents), value = TRUE)
  immune_ids <- grep("Immune_", levels(cellchat@idents), value = TRUE)
  group_sizes <- as.numeric(table(cellchat@idents)[rownames(cellchat@net$count)])

  ## network circle plots
  pdf(paste0(prefix, "_LEC_Immune_Networks.pdf"), width = 14, height = 12)
  netVisual_circle(cellchat@net$count,  vertex.weight = group_sizes, weight.scale = TRUE,
                   label.edge = FALSE, title.name = paste0("Number of Interactions (", prefix, ")"))
  netVisual_circle(cellchat@net$weight, vertex.weight = group_sizes, weight.scale = TRUE,
                   label.edge = FALSE, title.name = paste0("Interaction Strength (", prefix, ")"))
  dev.off()

  ## overall heatmap + top pathways
  pdf(paste0(prefix, "_LEC_Immune_Pathways.pdf"), width = 16, height = 12)
  netVisual_heatmap(cellchat, signaling = NULL, color.heatmap = "Reds")
  for (pathway in head(cellchat@netP$pathways, n_pathways)) {
    tryCatch(netVisual_aggregate(cellchat, signaling = pathway, layout = "circle"),
             error = function(e) cat("Skipping pathway", pathway, "due to error\n"))
  }
  dev.off()

  ## directional tables
  lec_to_immune <- subsetCommunication(cellchat, sources.use = lec_ids,    targets.use = immune_ids)
  immune_to_lec <- subsetCommunication(cellchat, sources.use = immune_ids, targets.use = lec_ids)
  top_l2i <- lec_to_immune %>% arrange(desc(prob)) %>% head(top_n)
  top_i2l <- immune_to_lec %>% arrange(desc(prob)) %>% head(top_n)
  write.csv(lec_to_immune, paste0(prefix, "_LEC_to_Immune_Communications.csv"), row.names = FALSE)
  write.csv(immune_to_lec, paste0(prefix, "_Immune_to_LEC_Communications.csv"), row.names = FALSE)
  write.csv(top_l2i, paste0(prefix, "_Top", top_n, "_LEC_to_Immune_Interactions.csv"), row.names = FALSE)
  write.csv(top_i2l, paste0(prefix, "_Top", top_n, "_Immune_to_LEC_Interactions.csv"), row.names = FALSE)
  saveRDS(cellchat, paste0(prefix, "_LEC_Immune_CellChat.rds"))

  ## summary
  cat("\n===", prefix, "LEC-IMMUNE CELLCHAT SUMMARY ===\n")
  cat("Cell types:", nlevels(cellchat@idents), " | LEC -> Immune:", nrow(lec_to_immune),
      " | Immune -> LEC:", nrow(immune_to_lec), " | pathways:", length(cellchat@netP$pathways), "\n")
  show_top <- function(tab, title) {
    cat("\nTop 5", title, ":\n")
    for (i in seq_len(min(5, nrow(tab)))) {
      cat(paste0(i, ". ", tab$source[i], " -> ", tab$target[i], " (", tab$ligand[i], "-", tab$receptor[i],
                 ", prob=", round(tab$prob[i], 3), ")\n"))
    }
  }
  show_top(top_l2i, "LEC-to-Immune interactions")
  show_top(top_i2l, "Immune-to-LEC interactions")
  for (lec in lec_ids) {
    spec <- lec_to_immune %>% filter(source == lec) %>% arrange(desc(prob)) %>% head(3)
    if (nrow(spec) > 0) {
      cat("\n", lec, "top targets:\n")
      for (j in seq_len(nrow(spec))) {
        cat("  ", spec$target[j], " (", spec$ligand[j], "-", spec$receptor[j], ", prob=", round(spec$prob[j], 3), ")\n")
      }
    }
  }
  signaling_summary <- lec_to_immune %>% group_by(source) %>%
    summarise(n_targets = n_distinct(target), avg_prob = round(mean(prob), 3), max_prob = round(max(prob), 3),
              top_target = target[which.max(prob)], .groups = "drop") %>% arrange(desc(avg_prob))
  print(signaling_summary)
  invisible(cellchat)
}

#' LEC -> immune signalling heatmaps for a CellChat object produced above.
#' Draws (1) CellChat's own directional heatmap and bubble plot with isolated
#' nodes removed, (2) a pheatmap of the weight matrix restricted to LEC rows
#' and immune columns, with and without clustering of the immune columns.
plot_lec_to_immune_heatmaps <- function(cellchat, prefix = "Targeted", row_order = NULL,
                                        n_pathways = 6) {
  .load_cellchat_packages()
  lec_ids    <- grep("LEC_",    levels(cellchat@idents), value = TRUE)
  immune_ids <- grep("Immune_", levels(cellchat@idents), value = TRUE)
  lec_to_immune <- subsetCommunication(cellchat, sources.use = lec_ids, targets.use = immune_ids)

  pdf(paste0(prefix, "_LEC_Signaling_Heatmaps.pdf"), width = 18, height = 8)
  if (nrow(lec_to_immune) > 0) {
    netVisual_heatmap(cellchat, sources.use = lec_ids, targets.use = immune_ids, color.heatmap = "Reds",
                      title.name = "LEC Signaling to Immune Cells", remove.isolate = TRUE, width = 12, height = 6)
    netVisual_bubble(cellchat, sources.use = lec_ids, targets.use = immune_ids, remove.isolate = TRUE,
                     title.name = "LEC -> Immune Signaling Strength")
  }
  for (pathway in head(cellchat@netP$pathways, n_pathways)) {
    tryCatch(netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", sources.use = lec_ids,
                                 remove.isolate = TRUE),
             error = function(e) cat("Skipping pathway", pathway, "due to error\n"))
  }
  dev.off()

  ## custom matrix: LECs as rows, immune clusters as columns
  comm_matrix <- cellchat@net$weight
  m <- comm_matrix[grep("LEC_", rownames(comm_matrix), value = TRUE),
                   grep("Immune_", colnames(comm_matrix), value = TRUE), drop = FALSE]
  rownames(m) <- gsub("^LEC_LEC", "LEC", rownames(m))
  if (!is.null(row_order)) m <- m[intersect(row_order, rownames(m)), , drop = FALSE]
  cat("LEC -> immune weight matrix:", nrow(m), "LEC rows x", ncol(m), "immune columns\n")

  pdf(paste0(prefix, "_LEC_Signaling_Matrix.pdf"), width = 20, height = 6)
  pheatmap(m, color = colorRampPalette(c("white", "red"))(100), cluster_rows = FALSE, cluster_cols = FALSE,
           show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 12, fontsize_col = 8,
           main = "LEC Signaling Strength to Immune Cells", angle_col = 45, cellwidth = 15, cellheight = 40,
           border_color = "grey90")
  pheatmap(m, color = colorRampPalette(c("#FDE725", "white", "#440154"))(100), cluster_rows = FALSE,
           cluster_cols = TRUE, show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 12, fontsize_col = 8,
           main = "LEC -> Immune (Clustered by Response)", angle_col = 45, cellwidth = 15, cellheight = 40)
  nz <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]
  if (nrow(nz) > 0 && ncol(nz) > 0) {
    pheatmap(nz, color = colorRampPalette(c("white", "red"))(50), cluster_rows = FALSE, cluster_cols = FALSE,
             main = "LEC -> Immune Signaling (Non-zero only)", fontsize = 10, cellwidth = 15, cellheight = 15)
  }
  dev.off()
  invisible(m)
}
