## Donor-level pseudobulk differential expression between two LEC subclusters
## (DESeq2), with volcano plots. One function replaces the three near-identical
## notebook blocks for LEC1 vs LEC3, LEC1 vs LEC6 and LEC1 vs LEC5.
##
## Design: cells are summed per donor x subcluster (>= 5 cells per pseudobulk
## sample); DESeq2 with design ~ cluster; contrast test vs ref. Each comparison
## is run twice: on all donors ("ALL") and on IPF + SSc-ILD donors only
## ("DISEASED").

## Packages are attached when the first function is called, so sourcing this file
## (which config.R does for every stage) does not require DESeq2 to be installed.
.load_pseudobulk_packages <- function() {
  suppressPackageStartupMessages({
    library(Seurat); library(DESeq2); library(ggplot2); library(ggrepel); library(dplyr)
  })
}

#' Sum raw counts per donor x LEC_subcluster.
create_cluster_pseudobulk <- function(seurat_obj, comparison_name, min_cells = 5,
                                      cluster_col = "LEC_subcluster") {
  .load_pseudobulk_packages()
  expr_matrix <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
  metadata <- seurat_obj@meta.data
  metadata$cell_id   <- rownames(metadata)
  metadata$sample_id <- paste(metadata$donor, metadata[[cluster_col]], sep = "_")

  pseudobulk_list <- list()
  sample_metadata <- data.frame()
  for (sample in unique(metadata$sample_id)) {
    cells_in_sample <- metadata$cell_id[metadata$sample_id == sample]
    if (length(cells_in_sample) >= min_cells) {
      pseudobulk_list[[sample]] <- rowSums(expr_matrix[, cells_in_sample, drop = FALSE])
      sample_info <- metadata[metadata$sample_id == sample, ][1, ]
      sample_metadata <- rbind(sample_metadata, data.frame(
        sample_id  = sample,
        donor      = sample_info$donor,
        cluster    = sample_info[[cluster_col]],
        cell_count = length(cells_in_sample),
        comparison = comparison_name
      ))
    }
  }
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  colnames(pseudobulk_matrix) <- names(pseudobulk_list)
  list(counts = pseudobulk_matrix, metadata = sample_metadata)
}

#' DESeq2 on a pseudobulk matrix, contrast test vs ref. Genes need >= 10 counts
#' in >= 3 samples.
run_cluster_deseq2 <- function(counts, metadata, comparison_name, ref, test) {
  .load_pseudobulk_packages()
  keep_genes <- rowSums(counts >= 10) >= 3
  counts <- counts[keep_genes, ]
  cat("Analyzing", nrow(counts), "genes across", ncol(counts), "samples for", comparison_name, "\n")
  metadata$cluster <- factor(as.character(metadata$cluster), levels = c(ref, test))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ cluster)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("cluster", test, ref), alpha = 0.05)
  as.data.frame(res) %>%
    mutate(gene = rownames(.), comparison = comparison_name) %>%
    arrange(padj, desc(abs(log2FoldChange)))
}

#' Volcano plot with the top genes labelled on both sides.
create_enhanced_volcano <- function(de_results, title, ref, test, top_genes_to_label = 15) {
  plot_data <- de_results %>%
    mutate(
      significance = case_when(
        is.na(padj) ~ "Not significant",
        padj >= 0.05 ~ "Not significant",
        padj < 0.05 & abs(log2FoldChange) < 0.5 ~ "Significant (low FC)",
        padj < 0.05 & abs(log2FoldChange) >= 0.5 & abs(log2FoldChange) < 1 ~ "Significant",
        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Highly significant"
      ),
      neg_log10_padj = pmin(-log10(padj + 1e-300), 50),
      direction = ifelse(log2FoldChange > 0, paste(test, "higher"), paste(ref, "higher"))
    ) %>%
    filter(!is.na(padj))

  top_test <- plot_data %>% filter(padj < 0.05, log2FoldChange > 0) %>% arrange(padj) %>% head(top_genes_to_label / 2)
  top_ref  <- plot_data %>% filter(padj < 0.05, log2FoldChange < 0) %>% arrange(padj) %>% head(top_genes_to_label / 2)
  label_data <- plot_data[plot_data$gene %in% unique(c(top_test$gene, top_ref$gene)), ]

  n_up_test <- sum(plot_data$padj < 0.05 & plot_data$log2FoldChange > 0, na.rm = TRUE)
  n_up_ref  <- sum(plot_data$padj < 0.05 & plot_data$log2FoldChange < 0, na.rm = TRUE)
  subtitle  <- paste0(test, " higher: ", n_up_test, " genes | ", ref, " higher: ", n_up_ref, " genes")
  x_lim <- max(abs(plot_data$log2FoldChange), na.rm = TRUE) * 1.1

  ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Not significant" = "gray80", "Significant (low FC)" = "lightblue",
                                  "Significant" = "orange", "Highly significant" = "red")) +
    geom_text_repel(data = label_data, aes(label = gene), size = 3.5, max.overlaps = 25,
                    box.padding = 0.5, point.padding = 0.2) +
    geom_vline(xintercept = c(-1, -0.5, 0.5, 1), linetype = c("dashed", "dotted", "dotted", "dashed"),
               alpha = 0.6, color = c("red", "blue", "blue", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6, color = "black") +
    annotate("text", x = -Inf, y = Inf, label = paste("Total genes:", nrow(plot_data)),
             hjust = -0.1, vjust = 1.2, size = 3, color = "gray50") +
    labs(title = title, subtitle = subtitle,
         x = paste0("Log2 Fold Change (", test, " vs ", ref, ")"),
         y = "-Log10 (Adjusted P-value)", color = "Significance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 12), legend.position = "bottom",
          panel.grid.minor = element_blank()) +
    xlim(c(-x_lim, x_lim))
}

summarize_de_results <- function(results, comparison_name, ref, test) {
  sig <- results[!is.na(results$padj) & results$padj < 0.05, ]
  cat("\n", comparison_name, ":\n")
  cat("Total genes tested:", nrow(results), "\n")
  cat("Significant genes (padj < 0.05):", nrow(sig), "\n")
  cat(test, "higher:", sum(sig$log2FoldChange > 0), "\n")
  cat(ref,  "higher:", sum(sig$log2FoldChange < 0), "\n")
  cat("Highly significant (|FC| >= 1):", sum(abs(sig$log2FoldChange) >= 1), "\n")
  for (side in list(list(test, sig[sig$log2FoldChange > 0, ]), list(ref, sig[sig$log2FoldChange < 0, ]))) {
    top <- side[[2]] %>% arrange(padj) %>% head(5)
    if (nrow(top) > 0) {
      cat("\nTop genes higher in", side[[1]], ":\n")
      for (i in seq_len(nrow(top))) {
        cat(sprintf("  %s: FC=%.2f, padj=%.2e\n", top$gene[i], abs(top$log2FoldChange[i]), top$padj[i]))
      }
    }
  }
  list(total = nrow(results), significant = nrow(sig),
       test_higher = sum(sig$log2FoldChange > 0), ref_higher = sum(sig$log2FoldChange < 0),
       high_fc = sum(abs(sig$log2FoldChange) >= 1))
}

#' Significant gene lists (padj < 0.05, |log2FC| > fc_threshold), split by direction.
extract_gene_lists <- function(de_results, comparison_name, ref, test,
                               padj_threshold = 0.05, fc_threshold = 0.25, save = TRUE) {
  sig <- de_results[!is.na(de_results$padj) & de_results$padj < padj_threshold &
                      abs(de_results$log2FoldChange) > fc_threshold, ]
  test_up <- sig$gene[sig$log2FoldChange > 0]
  ref_up  <- sig$gene[sig$log2FoldChange < 0]
  cat(comparison_name, ":\n")
  cat("  ", test, "upregulated genes:", length(test_up), "\n")
  cat("  ", ref,  "upregulated genes:", length(ref_up), "\n")
  if (save) {
    tag <- gsub(" ", "_", comparison_name)
    if (length(test_up) > 0) write.csv(data.frame(gene = test_up), paste0(tag, "_", test, "_upregulated_genes.csv"), row.names = FALSE)
    if (length(ref_up)  > 0) write.csv(data.frame(gene = ref_up),  paste0(tag, "_", ref,  "_upregulated_genes.csv"), row.names = FALSE)
  }
  list(test_up = test_up, ref_up = ref_up, all_sig = sig$gene)
}

#' Run the full ref-vs-test pseudobulk comparison (ALL and DISEASED donors),
#' write the CSVs and volcano PDFs, and return both result tables.
run_cluster_vs_cluster_de <- function(sce_LECs, ref = "LEC1", test = "LEC3",
                                      disease_groups = c("IPF", "SSc-ILD"), top_genes_to_label = 20) {
  .load_pseudobulk_packages()
  tag <- paste0(ref, "_vs_", test)
  if (packageVersion("Seurat") >= "5.0.0") sce_LECs[["RNA"]] <- JoinLayers(sce_LECs[["RNA"]])

  ## ---- ALL donors ----
  cat("\n=== COMPARISON 1: ALL", ref, "vs ALL", test, "(All Conditions) ===\n")
  sce_all <- subset(sce_LECs, subset = LEC_subcluster %in% c(ref, test))
  print(table(sce_all$LEC_subcluster, sce_all$group))
  all_pb <- create_cluster_pseudobulk(sce_all, "All_Conditions")
  cat("Pseudobulk samples:", ncol(all_pb$counts), "\n"); print(table(all_pb$metadata$cluster))
  all_results <- run_cluster_deseq2(all_pb$counts, all_pb$metadata, paste0("All_", tag), ref, test)

  ## ---- DISEASED donors ----
  cat("\n=== COMPARISON 2: DISEASED", ref, "vs", test, "(", paste(disease_groups, collapse = " + "), "only ) ===\n")
  sce_dis <- subset(sce_LECs, subset = LEC_subcluster %in% c(ref, test) & group %in% disease_groups)
  print(table(sce_dis$LEC_subcluster, sce_dis$group))
  dis_pb <- create_cluster_pseudobulk(sce_dis, "Diseased_Only")
  cat("Pseudobulk samples:", ncol(dis_pb$counts), "\n"); print(table(dis_pb$metadata$cluster))
  diseased_results <- run_cluster_deseq2(dis_pb$counts, dis_pb$metadata, paste0("Diseased_", tag), ref, test)

  ## ---- volcano plots ----
  volcano_all <- create_enhanced_volcano(all_results, paste("ALL", ref, "vs", test, "(All Conditions)"),
                                         ref, test, top_genes_to_label)
  volcano_dis <- create_enhanced_volcano(diseased_results, paste("DISEASED", ref, "vs", test, "(IPF + SSc-ILD only)"),
                                         ref, test, top_genes_to_label)
  print(volcano_all); print(volcano_dis)
  ggsave(paste0("ALL_", tag, "_Volcano.pdf"),      volcano_all, width = 10, height = 8, dpi = 300)
  ggsave(paste0("DISEASED_", tag, "_Volcano.pdf"), volcano_dis, width = 10, height = 8, dpi = 300)

  ## ---- summaries and tables ----
  all_summary <- summarize_de_results(all_results,      paste("ALL", tag),      ref, test)
  dis_summary <- summarize_de_results(diseased_results, paste("DISEASED", tag), ref, test)
  write.csv(all_results,      paste0("ALL_", tag, "_DE_results.csv"),      row.names = FALSE)
  write.csv(diseased_results, paste0("DISEASED_", tag, "_DE_results.csv"), row.names = FALSE)
  summary_stats <- data.frame(
    Comparison        = c(paste0("ALL_", tag), paste0("DISEASED_", tag)),
    Total_genes       = c(all_summary$total,       dis_summary$total),
    Significant_genes = c(all_summary$significant, dis_summary$significant),
    Test_higher       = c(all_summary$test_higher, dis_summary$test_higher),
    Ref_higher        = c(all_summary$ref_higher,  dis_summary$ref_higher),
    High_FC_genes     = c(all_summary$high_fc,     dis_summary$high_fc)
  )
  names(summary_stats)[4:5] <- c(paste0(test, "_higher"), paste0(ref, "_higher"))
  write.csv(summary_stats, paste0(tag, "_Volcano_Summary_Stats.csv"), row.names = FALSE)
  top_all <- all_results      %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj) %>% head(50)
  top_dis <- diseased_results %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj) %>% head(50)
  write.csv(top_all, paste0("TOP_50_ALL_", tag, ".csv"),      row.names = FALSE)
  write.csv(top_dis, paste0("TOP_50_DISEASED_", tag, ".csv"), row.names = FALSE)

  ## ---- gene lists for GO + overlap between ALL and DISEASED ----
  cat("\n=== EXTRACTING GENE LISTS FOR GO ANALYSIS ===\n")
  all_lists <- extract_gene_lists(all_results,      paste("ALL", tag),      ref, test)
  dis_lists <- extract_gene_lists(diseased_results, paste("DISEASED", tag), ref, test)
  cat("\n=== COMPARATIVE ANALYSIS: ALL vs DISEASED ===\n")
  for (side in c("test_up", "ref_up")) {
    a <- all_lists[[side]]; d <- dis_lists[[side]]
    if (length(a) > 0 && length(d) > 0) {
      cat(if (side == "test_up") test else ref, "upregulated genes:\n")
      cat("  Shared (ALL & DISEASED):", length(intersect(a, d)), "\n")
      cat("  ALL conditions only:", length(setdiff(a, d)), "\n")
      cat("  DISEASED only:", length(setdiff(d, a)), "\n")
    }
  }
  cat("\n===", tag, "PSEUDOBULK ANALYSIS COMPLETE ===\n")
  invisible(list(all = all_results, diseased = diseased_results,
                 gene_lists = list(all = all_lists, diseased = dis_lists)))
}
