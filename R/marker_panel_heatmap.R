## Single-cell marker-panel heatmap (EC vs immune) used to judge how clean an
## EC / immune split is after each filtering step.
##
## The original notebook repeated this block five times, once per object
## (combined, cluster-filtered, clean, ultra-clean, cluster-filtered-EC). It is
## now one function; the logic is unchanged.

suppressPackageStartupMessages({
  library(Seurat)
  library(pheatmap)
})

## Gene panel: lymphatic, arterial, venous and core EC markers, then immune markers
marker_panel <- list(
  LEC      = c("PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "VEGFC"),
  Arterial = c("FABP4", "GJA5", "EFNB2", "DLL4", "SOX17", "HEY2"),
  Venous   = c("NR2F2", "EPHB4"),
  Core_EC  = c("PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EMCN", "ESAM"),
  Immune   = c("PTPRC", "CD3D", "CD8A", "CD4", "NKG7", "MS4A1", "MRC1", "CD163", "LYZ",
               "FCGR3A", "MS4A7", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD68")
)

marker_panel_category_colors <- c("LEC" = "#2E8B57", "Arterial" = "#FF6347", "Venous" = "#4169E1",
                                  "Core_EC" = "#E31A1C", "Immune" = "#1F78B4")

#' Pull log-normalized expression for a gene set, Seurat v4/v5 safe, without
#' densifying the whole assay.
fetch_normalized_matrix <- function(obj, genes, cells, assay = DefaultAssay(obj)) {
  if (inherits(obj[[assay]], "Assay5")) {
    df <- FetchData(obj, vars = genes, cells = cells, layer = "data")
    t(as.matrix(df[, genes, drop = FALSE]))
  } else {
    mat <- GetAssayData(obj, assay = assay, slot = "data")
    as.matrix(mat[genes, cells, drop = FALSE])
  }
}

#' Row z-score with clipping at +/- clip for contrast.
row_zscore <- function(M, clip = 2.5) {
  row_means <- rowMeans(M, na.rm = TRUE)
  row_sds   <- sqrt(pmax(rowMeans((M - row_means)^2), 0))
  row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1
  z <- (M - row_means) / row_sds
  z[z >  clip] <-  clip
  z[z < -clip] <- -clip
  z[!is.finite(z)] <- 0
  z
}

#' Plot the EC / immune marker-panel heatmap for a Seurat object.
#'
#' @param obj   Seurat object with a `cellset` metadata column ("EC" / "Immune").
#' @param label short label used in titles and file names, e.g. "Clean".
#' @param panel named list of gene vectors (default: marker_panel).
#' @param save_pdf write <label>_Expression_Heatmap.pdf and
#'        Enhanced_<label>_Expression_Heatmap.pdf into the working directory.
#' @return invisibly, a list with the z-score matrix and the ordered cell types.
plot_marker_panel_heatmap <- function(obj, label, panel = marker_panel, save_pdf = TRUE,
                                      title = paste0("Expression (Z-score) - ", label)) {
  if (!"cellset" %in% colnames(obj@meta.data)) {
    stop("obj@meta.data$cellset not found (expected: 'EC' or 'Immune').")
  }
  cells_use <- colnames(obj)[obj$cellset %in% c("EC", "Immune")]
  stopifnot(length(cells_use) > 0)

  genes_panel <- unique(unlist(panel, use.names = FALSE))
  genes_avail <- genes_panel[genes_panel %in% rownames(obj[[DefaultAssay(obj)]])]
  if (length(genes_avail) == 0) stop("None of the panel genes are present in the assay rows.")

  M <- fetch_normalized_matrix(obj, genes_avail, cells_use)
  z <- row_zscore(M)

  ## order columns EC -> Immune with a visual gap between them
  ct  <- factor(obj$cellset[cells_use], levels = c("EC", "Immune"))
  ord <- order(ct, na.last = TRUE)
  z_ord  <- z[, ord, drop = FALSE]
  ct_ord <- droplevels(ct[ord])
  gap_at <- if (length(unique(ct_ord)) == 2) sum(ct_ord == "EC") else NULL

  anno_col <- data.frame(`Cell type` = ct_ord, check.names = FALSE)
  rownames(anno_col) <- colnames(z_ord)
  anno_colors <- list(`Cell type` = cellset_colors)

  ## row annotation: which part of the panel each gene belongs to
  gene_categories <- data.frame(
    Gene_Type = unlist(lapply(names(panel), function(nm) rep(nm, sum(panel[[nm]] %in% genes_avail)))),
    row.names = unlist(lapply(panel, function(g) g[g %in% genes_avail]), use.names = FALSE)
  )
  gene_categories <- gene_categories[genes_avail, , drop = FALSE]
  gene_category_colors <- list(Gene_Type = marker_panel_category_colors[names(panel)])

  breaks <- seq(-2.5, 2.5, length.out = 101)
  cols   <- grDevices::colorRampPalette(c("#FDE725", "white", "#440154"))(length(breaks) - 1)

  cat("\n=== HEATMAP SUMMARY (", label, ") ===\n", sep = "")
  cat("Total cells used:", length(cells_use), "\n")
  cat("EC cells:", sum(ct_ord == "EC"), "\n")
  cat("Immune cells:", sum(ct_ord == "Immune"), "\n")
  cat("Genes available:", length(genes_avail), "of", length(genes_panel), "\n")
  missing_genes <- setdiff(genes_panel, genes_avail)
  if (length(missing_genes) > 0) cat("Missing gene names:", paste(missing_genes, collapse = ", "), "\n")

  draw <- function(with_gene_categories, fontsize_row) {
    pheatmap(
      z_ord,
      color             = cols,
      breaks            = breaks,
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      show_colnames     = FALSE,
      show_rownames     = TRUE,
      fontsize_row      = fontsize_row,
      border_color      = NA,
      annotation_col    = anno_col,
      annotation_row    = if (with_gene_categories) gene_categories else NA,
      annotation_colors = if (with_gene_categories) c(anno_colors, gene_category_colors) else anno_colors,
      gaps_col          = gap_at,
      main              = if (with_gene_categories) paste(title, "with Gene Categories") else title,
      legend_breaks     = c(-2.5, 0, 2.5),
      legend_labels     = c("-2.5", "0", "2.5"),
      use_raster        = TRUE
    )
  }

  draw(FALSE, 9)
  if (save_pdf) {
    pdf(paste0(label, "_Expression_Heatmap.pdf"), height = 10, width = 12)
    draw(FALSE, 9)
    dev.off()
    pdf(paste0("Enhanced_", label, "_Expression_Heatmap.pdf"), height = 12, width = 14)
    draw(TRUE, 8)
    dev.off()
    cat("Saved:", paste0(label, "_Expression_Heatmap.pdf"), "and",
        paste0("Enhanced_", label, "_Expression_Heatmap.pdf"), "\n")
  }
  invisible(list(z = z_ord, celltype = ct_ord, genes = genes_avail))
}

#' Print how many EC / immune cells a cleaning step removed.
report_cleaning_delta <- function(before, after, before_label = "BEFORE", after_label = "AFTER") {
  b_ec <- sum(before$cellset == "EC");  b_im <- sum(before$cellset == "Immune")
  a_ec <- sum(after$cellset  == "EC");  a_im <- sum(after$cellset  == "Immune")
  cat("\n=== COMPARISON SUMMARY ===\n")
  cat(before_label, ":\n  EC cells:", b_ec, "\n  Immune cells:", b_im, "\n  Total:", b_ec + b_im, "\n")
  cat(after_label,  ":\n  EC cells:", a_ec, "\n  Immune cells:", a_im, "\n  Total:", a_ec + a_im, "\n")
  cat("CELLS REMOVED:\n")
  cat("  EC cells removed:", b_ec - a_ec, "(", round(100 * (b_ec - a_ec) / max(b_ec, 1), 1), "%)\n")
  cat("  Immune cells removed:", b_im - a_im, "(", round(100 * (b_im - a_im) / max(b_im, 1), 1), "%)\n")
}
