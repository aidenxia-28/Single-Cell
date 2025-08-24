## Lymphatic branch: donor-level pseudobulk DE between disease groups within LEC1 and within LEC3
## Stage 25 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: LEC1_DE_*.csv, LEC3_DE_*.csv, LEC1_vs_LEC3_*_Summary.csv, *_Volcano.pdf, LEC1_vs_LEC3_Top_Genes_Heatmap.pdf
## Counts are summed per donor x group x subcluster (>= 5 cells); DESeq2 with
## design ~ group; contrasts IPF vs CTRL, SSc-ILD vs CTRL and IPF vs SSc-ILD,
## run separately inside LEC1 and LEC3.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# =============================================================================
# LEC1 vs LEC3 Pseudobulk Differential Expression Analysis
# =============================================================================


library(Seurat)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(patchwork)

# =============================================================================
# LOAD DATA AND SUBSET CLUSTERS
# =============================================================================

cat("Loading LEC data and subsetting LEC1 and LEC3...\n")
sce_LECs <- readRDS("LECs_subclustered.rds")

# Subset to only LEC1 and LEC3
sce_comparison <- subset(sce_LECs, subset = LEC_subcluster %in% c("LEC1", "LEC3"))

cat("Total cells in LEC1 and LEC3:\n")
print(table(sce_comparison$LEC_subcluster, sce_comparison$group))

# Check we have sufficient donors per group
donor_summary <- sce_comparison@meta.data %>%
  group_by(group, LEC_subcluster, donor) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(group, LEC_subcluster) %>%
  summarise(donor_count = n(), total_cells = sum(cell_count), .groups = 'drop')

cat("\nDonor and cell counts per group and cluster:\n")
print(donor_summary)

# Join layers for Seurat v5 compatibility
if(packageVersion("Seurat") >= "5.0.0") {
  sce_comparison[["RNA"]] <- JoinLayers(sce_comparison[["RNA"]])
}

# =============================================================================
# CREATE PSEUDOBULK PROFILES
# =============================================================================

cat("\n=== CREATING PSEUDOBULK PROFILES ===\n")

# Function to create pseudobulk profiles
create_pseudobulk <- function(seurat_obj, group_by = c("donor", "group", "LEC_subcluster")) {

  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")

  # Get metadata
  metadata <- seurat_obj@meta.data[, group_by, drop = FALSE]
  metadata$cell_id <- rownames(metadata)

  # Create grouping variable
  metadata$sample_id <- paste(metadata$donor, metadata$group, metadata$LEC_subcluster, sep = "_")

  # Aggregate counts by sample
  pseudobulk_list <- list()
  sample_metadata <- data.frame()

  for(sample in unique(metadata$sample_id)) {
    # Get cells for this sample
    cells_in_sample <- metadata$cell_id[metadata$sample_id == sample]

    if(length(cells_in_sample) >= 5) {  # Minimum 5 cells per pseudobulk sample
      # Sum counts across cells
      pseudobulk_counts <- rowSums(expr_matrix[, cells_in_sample, drop = FALSE])
      pseudobulk_list[[sample]] <- pseudobulk_counts

      # Store sample metadata
      sample_info <- metadata[metadata$sample_id == sample, ][1, ]
      sample_metadata <- rbind(sample_metadata, data.frame(
        sample_id = sample,
        donor = sample_info$donor,
        group = sample_info$group,
        cluster = sample_info$LEC_subcluster,
        cell_count = length(cells_in_sample)
      ))
    }
  }

  # Convert to matrix
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  colnames(pseudobulk_matrix) <- names(pseudobulk_list)

  return(list(counts = pseudobulk_matrix, metadata = sample_metadata))
}

# Create pseudobulk profiles
pseudobulk_data <- create_pseudobulk(sce_comparison)
pseudobulk_counts <- pseudobulk_data$counts
pseudobulk_metadata <- pseudobulk_data$metadata

cat("Pseudobulk samples created:", ncol(pseudobulk_counts), "\n")
cat("Genes:", nrow(pseudobulk_counts), "\n")

# Show sample distribution
cat("\nPseudobulk sample distribution:\n")
print(table(pseudobulk_metadata$group, pseudobulk_metadata$cluster))

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ2
# =============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Function to run DESeq2 analysis
run_deseq2_analysis <- function(counts, metadata, cluster_subset, comparison_groups) {

  # Subset to specific cluster
  cluster_samples <- metadata$sample_id[metadata$cluster == cluster_subset]
  cluster_counts <- counts[, cluster_samples, drop = FALSE]
  cluster_metadata <- metadata[metadata$cluster == cluster_subset, ]

  # Filter lowly expressed genes (at least 10 counts in at least 3 samples)
  keep_genes <- rowSums(cluster_counts >= 10) >= 3
  cluster_counts <- cluster_counts[keep_genes, ]

  cat("For", cluster_subset, "- analyzing", nrow(cluster_counts), "genes across", ncol(cluster_counts), "samples\n")

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = cluster_counts,
    colData = cluster_metadata,
    design = ~ group
  )

  # Run DESeq2
  dds <- DESeq(dds)

  # Extract results for each comparison
  results_list <- list()

  for(comp in comparison_groups) {
    contrast <- c("group", comp[1], comp[2])
    res <- results(dds, contrast = contrast, alpha = 0.05)
    res_df <- as.data.frame(res) %>%
      mutate(
        gene = rownames(.),
        cluster = cluster_subset,
        comparison = paste(comp[1], "vs", comp[2])
      ) %>%
      arrange(padj, desc(abs(log2FoldChange)))

    results_list[[paste(comp[1], "vs", comp[2])]] <- res_df
  }

  return(results_list)
}

# Define comparisons
comparisons <- list(
  c("IPF", "CTRL"),
  c("SSc-ILD", "CTRL"),
  c("IPF", "SSc-ILD")
)

# Run analysis for both clusters
lec1_results <- run_deseq2_analysis(pseudobulk_counts, pseudobulk_metadata, "LEC1", comparisons)
lec3_results <- run_deseq2_analysis(pseudobulk_counts, pseudobulk_metadata, "LEC3", comparisons)

# =============================================================================
# EXTRACT TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

cat("\n=== TOP DIFFERENTIALLY EXPRESSED GENES ===\n")

# Function to extract top genes
extract_top_genes <- function(results_list, top_n = 20) {

  top_genes_summary <- data.frame()

  for(comp_name in names(results_list)) {
    res <- results_list[[comp_name]]

    # Get significant genes
    sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]

    if(nrow(sig_genes) > 0) {
      # Top upregulated
      top_up <- sig_genes[sig_genes$log2FoldChange > 0, ][1:min(top_n, sum(sig_genes$log2FoldChange > 0)), ]
      if(nrow(top_up) > 0) {
        top_up$direction <- "Upregulated"
        top_up$rank <- 1:nrow(top_up)
      }

      # Top downregulated
      top_down <- sig_genes[sig_genes$log2FoldChange < 0, ][1:min(top_n, sum(sig_genes$log2FoldChange < 0)), ]
      if(nrow(top_down) > 0) {
        top_down$direction <- "Downregulated"
        top_down$rank <- 1:nrow(top_down)
      }

      # Combine
      if(nrow(top_up) > 0 && nrow(top_down) > 0) {
        top_genes <- rbind(top_up, top_down)
      } else if(nrow(top_up) > 0) {
        top_genes <- top_up
      } else if(nrow(top_down) > 0) {
        top_genes <- top_down
      } else {
        next
      }

      top_genes$comparison <- comp_name
      top_genes_summary <- rbind(top_genes_summary, top_genes)
    }
  }

  return(top_genes_summary)
}

# Extract top genes for both clusters
lec1_top_genes <- extract_top_genes(lec1_results, top_n = 15)
lec3_top_genes <- extract_top_genes(lec3_results, top_n = 15)

# Add cluster information
lec1_top_genes$cluster <- "LEC1"
lec3_top_genes$cluster <- "LEC3"

# Combine results
all_top_genes <- rbind(lec1_top_genes, lec3_top_genes)

# =============================================================================
# COMPARATIVE ANALYSIS BETWEEN LEC1 AND LEC3
# =============================================================================

cat("\n=== COMPARING LEC1 vs LEC3 RESPONSES ===\n")

# Function to compare responses between clusters
compare_cluster_responses <- function(lec1_res, lec3_res, comparison_name) {

  # Get significant genes from each cluster
  lec1_sig <- lec1_res[!is.na(lec1_res$padj) & lec1_res$padj < 0.05, ]
  lec3_sig <- lec3_res[!is.na(lec3_res$padj) & lec3_res$padj < 0.05, ]

  # Find common and unique genes
  lec1_genes <- lec1_sig$gene
  lec3_genes <- lec3_sig$gene

  common_genes <- intersect(lec1_genes, lec3_genes)
  lec1_unique <- setdiff(lec1_genes, lec3_genes)
  lec3_unique <- setdiff(lec3_genes, lec1_genes)

  # For common genes, compare fold changes
  common_comparison <- data.frame()
  if(length(common_genes) > 0) {
    lec1_common <- lec1_sig[lec1_sig$gene %in% common_genes, ]
    lec3_common <- lec3_sig[lec3_sig$gene %in% common_genes, ]

    # Merge by gene
    common_merged <- merge(lec1_common, lec3_common, by = "gene", suffixes = c("_LEC1", "_LEC3"))

    common_comparison <- common_merged %>%
      mutate(
        comparison = comparison_name,
        fold_change_diff = log2FoldChange_LEC1 - log2FoldChange_LEC3,
        response_pattern = case_when(
          sign(log2FoldChange_LEC1) == sign(log2FoldChange_LEC3) ~ "Same_direction",
          TRUE ~ "Opposite_direction"
        )
      ) %>%
      arrange(desc(abs(fold_change_diff)))
  }

  return(list(
    common_genes = common_genes,
    lec1_unique = lec1_unique,
    lec3_unique = lec3_unique,
    common_comparison = common_comparison,
    stats = c(
      LEC1_sig = length(lec1_genes),
      LEC3_sig = length(lec3_genes),
      Common = length(common_genes),
      LEC1_unique = length(lec1_unique),
      LEC3_unique = length(lec3_unique)
    )
  ))
}

# Compare responses for each disease comparison
comparison_results <- list()
for(comp_name in names(lec1_results)) {
  comparison_results[[comp_name]] <- compare_cluster_responses(
    lec1_results[[comp_name]],
    lec3_results[[comp_name]],
    comp_name
  )
}

# Print comparison summaries
for(comp_name in names(comparison_results)) {
  cat("\n", comp_name, "comparison:\n")
  print(comparison_results[[comp_name]]$stats)
}

# =============================================================================
# VISUALIZATION - VOLCANO PLOTS
# =============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Function to create volcano plot
create_volcano_plot <- function(de_results, title, top_genes_to_label = 10) {

  # Prepare data
  plot_data <- de_results %>%
    mutate(
      significance = case_when(
        is.na(padj) ~ "Not significant",
        padj >= 0.05 ~ "Not significant",
        padj < 0.05 & abs(log2FoldChange) < 1 ~ "Significant",
        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Highly significant"
      ),
      neg_log10_padj = -log10(padj + 1e-300)  # Add small constant to avoid -Inf
    )

  # Get top genes to label
  top_genes <- plot_data %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(top_genes_to_label)

  # Create plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c(
      "Not significant" = "gray70",
      "Significant" = "lightblue",
      "Highly significant" = "red"
    )) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 15) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 (Adjusted P-value)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# Create volcano plots for key comparisons
volcano_plots <- list()

for(comp_name in c("IPF vs CTRL", "SSc-ILD vs CTRL")) {
  volcano_plots[[paste("LEC1", comp_name)]] <- create_volcano_plot(
    lec1_results[[comp_name]],
    paste("LEC1:", comp_name)
  )

  volcano_plots[[paste("LEC3", comp_name)]] <- create_volcano_plot(
    lec3_results[[comp_name]],
    paste("LEC3:", comp_name)
  )
}

# Display volcano plots
for(plot_name in names(volcano_plots)) {
  print(volcano_plots[[plot_name]])
}

# =============================================================================
# HEATMAP OF TOP GENES
# =============================================================================

# Get top upregulated genes from each comparison
top_up_genes <- all_top_genes %>%
  filter(direction == "Upregulated", rank <= 10) %>%
  group_by(cluster, comparison) %>%
  slice_head(n = 8) %>%
  ungroup()

if(nrow(top_up_genes) > 0) {

  # Get unique genes for heatmap
  heatmap_genes <- unique(top_up_genes$gene)

  # Create expression matrix for heatmap
  heatmap_expr <- GetAssayData(sce_comparison, slot = "data", assay = "RNA")[heatmap_genes, , drop = FALSE]

  # Calculate average expression by cluster and group
  avg_expr <- AverageExpression(sce_comparison,
                               features = heatmap_genes,
                               group.by = c("LEC_subcluster", "group"),
                               slot = "data")$RNA

  # Scale for heatmap
  avg_expr_scaled <- t(scale(t(avg_expr)))

  # Create annotations
  col_info <- data.frame(
    sample = colnames(avg_expr_scaled),
    cluster = sapply(strsplit(colnames(avg_expr_scaled), "_"), `[`, 1),
    group = sapply(strsplit(colnames(avg_expr_scaled), "_"), `[`, 2)
  )
  rownames(col_info) <- col_info$sample

  annotation_cols <- list(
    cluster = c("LEC1" = "#1f77b4", "LEC3" = "#ff7f0e"),
    group = c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")
  )

  # Create heatmap
  pheatmap(
    avg_expr_scaled,
    annotation_col = col_info[, c("cluster", "group"), drop = FALSE],
    annotation_colors = annotation_cols,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166AC", "white", "#762A83"))(100),
    main = "Top Upregulated Genes - LEC1 vs LEC3",
    fontsize_row = 8,
    fontsize_col = 10,
    filename = "LEC1_vs_LEC3_Top_Genes_Heatmap.pdf",
    width = 8,
    height = 10
  )
}

# =============================================================================
# SUMMARY TABLES AND REPORTS
# =============================================================================

cat("\n=== CREATING SUMMARY REPORTS ===\n")

# Create comprehensive summary table
create_summary_table <- function(cluster_results, cluster_name) {

  summary_table <- data.frame()

  for(comp_name in names(cluster_results)) {
    res <- cluster_results[[comp_name]]
    sig_res <- res[!is.na(res$padj) & res$padj < 0.05, ]

    if(nrow(sig_res) > 0) {
      # Top 10 upregulated
      top_up <- sig_res[sig_res$log2FoldChange > 0, ] %>% head(10)
      if(nrow(top_up) > 0) {
        top_up$direction <- "Upregulated"
        top_up$cluster <- cluster_name
        top_up$comparison <- comp_name
        summary_table <- rbind(summary_table, top_up[, c("gene", "log2FoldChange", "padj", "direction", "cluster", "comparison")])
      }

      # Top 10 downregulated
      top_down <- sig_res[sig_res$log2FoldChange < 0, ] %>% head(10)
      if(nrow(top_down) > 0) {
        top_down$direction <- "Downregulated"
        top_down$cluster <- cluster_name
        top_down$comparison <- comp_name
        summary_table <- rbind(summary_table, top_down[, c("gene", "log2FoldChange", "padj", "direction", "cluster", "comparison")])
      }
    }
  }

  return(summary_table)
}

# Create summary tables
lec1_summary <- create_summary_table(lec1_results, "LEC1")
lec3_summary <- create_summary_table(lec3_results, "LEC3")
combined_summary <- rbind(lec1_summary, lec3_summary)

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save all DE results
for(comp_name in names(lec1_results)) {
  # LEC1 results
  write.csv(lec1_results[[comp_name]],
           paste0("LEC1_DE_", gsub(" ", "_", comp_name), ".csv"),
           row.names = FALSE)

  # LEC3 results
  write.csv(lec3_results[[comp_name]],
           paste0("LEC3_DE_", gsub(" ", "_", comp_name), ".csv"),
           row.names = FALSE)
}

# Save summary tables
write.csv(combined_summary, "LEC1_vs_LEC3_Top_DE_Genes_Summary.csv", row.names = FALSE)
write.csv(all_top_genes, "LEC1_vs_LEC3_All_Top_Genes.csv", row.names = FALSE)

# Save comparison results
comparison_stats <- data.frame()
for(comp_name in names(comparison_results)) {
  stats <- comparison_results[[comp_name]]$stats
  stats_df <- data.frame(
    Comparison = comp_name,
    LEC1_significant = stats["LEC1_sig"],
    LEC3_significant = stats["LEC3_sig"],
    Common_genes = stats["Common"],
    LEC1_unique = stats["LEC1_unique"],
    LEC3_unique = stats["LEC3_unique"]
  )
  comparison_stats <- rbind(comparison_stats, stats_df)
}

write.csv(comparison_stats, "LEC1_vs_LEC3_Comparison_Statistics.csv", row.names = FALSE)

# Save pseudobulk data
write.csv(pseudobulk_metadata, "Pseudobulk_Sample_Metadata.csv", row.names = FALSE)

# Save volcano plots
for(plot_name in names(volcano_plots)) {
  ggsave(paste0(gsub(" ", "_", plot_name), "_Volcano.pdf"),
         volcano_plots[[plot_name]], width = 8, height = 6)
}

cat("Files saved:\n")
cat("- Individual DE results: LEC1_DE_*.csv, LEC3_DE_*.csv\n")
cat("- Summary tables: LEC1_vs_LEC3_Top_DE_Genes_Summary.csv\n")
cat("- Comparison statistics: LEC1_vs_LEC3_Comparison_Statistics.csv\n")
cat("- Volcano plots: *_Volcano.pdf\n")
cat("- Top genes heatmap: LEC1_vs_LEC3_Top_Genes_Heatmap.pdf\n")

# =============================================================================
# KEY FINDINGS SUMMARY
# =============================================================================

cat("\n=== KEY FINDINGS SUMMARY ===\n")

for(comp_name in names(comparison_results)) {
  cat("\n", comp_name, ":\n")
  stats <- comparison_results[[comp_name]]$stats
  cat("- LEC1 significant genes:", stats["LEC1_sig"], "\n")
  cat("- LEC3 significant genes:", stats["LEC3_sig"], "\n")
  cat("- Common responding genes:", stats["Common"], "\n")
  cat("- LEC1-specific responses:", stats["LEC1_unique"], "\n")
  cat("- LEC3-specific responses:", stats["LEC3_unique"], "\n")
}

# Show top genes for deeper investigation
cat("\nTOP UPREGULATED GENES FOR FURTHER INVESTIGATION:\n")

top_candidates <- combined_summary %>%
  filter(direction == "Upregulated") %>%
  arrange(padj) %>%
  head(20)

if(nrow(top_candidates) > 0) {
  for(i in 1:min(10, nrow(top_candidates))) {
    gene_info <- top_candidates[i, ]
    cat(sprintf("%d. %s (%s, %s): FC=%.2f, padj=%.2e\n",
               i, gene_info$gene, gene_info$cluster, gene_info$comparison,
               gene_info$log2FoldChange, gene_info$padj))
  }
}

cat("\n=== PSEUDOBULK ANALYSIS COMPLETE ===\n")
