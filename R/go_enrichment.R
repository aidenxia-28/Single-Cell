## GO enrichment (clusterProfiler, org.Hs.eg.db) for the gene lists produced by
## run_cluster_vs_cluster_de(). One function replaces the two notebook blocks
## for LEC1 vs LEC3 and LEC1 vs LEC6.
##
## For each comparison four gene sets are tested against the universe of all
## tested genes: <test> up / <ref> up, each in ALL donors and in DISEASED donors.
## BP, MF and CC ontologies; BH-adjusted p < 0.05, q < 0.2.

## Packages are attached when the first function is called, so sourcing this file
## does not require clusterProfiler to be installed.
.load_go_packages <- function() {
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(dplyr); library(enrichplot)
  })
}

convert_to_entrez <- function(gene_symbols, comparison_name) {
  .load_go_packages()
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  if (length(gene_symbols) == 0) {
    cat("No genes to convert for", comparison_name, "\n")
    return(character(0))
  }
  conv <- suppressMessages(bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  cat(comparison_name, ": converted", nrow(conv), "out of", length(gene_symbols), "genes\n")
  conv$ENTREZID
}

run_go_enrichment <- function(entrez_ids, comparison_name, universe = NULL) {
  .load_go_packages()
  if (length(entrez_ids) < 5) {
    cat("Skipping", comparison_name, "- too few genes (", length(entrez_ids), ")\n")
    return(NULL)
  }
  cat("Running GO enrichment for", comparison_name, "with", length(entrez_ids), "genes\n")
  go_results <- list()
  for (ont in c("BP", "MF", "CC")) {
    tryCatch({
      ego <- enrichGO(gene = entrez_ids, universe = universe, OrgDb = org.Hs.eg.db, ont = ont,
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
      if (!is.null(ego) && nrow(ego) > 0) {
        go_results[[ont]] <- ego
        cat("  ", ont, ": found", nrow(ego), "enriched terms\n")
      }
    }, error = function(e) cat("  Error in", ont, "enrichment:", e$message, "\n"))
  }
  go_results
}

create_go_plots <- function(go_result_list, title_prefix, save_prefix) {
  .load_go_packages()
  plots <- list()
  for (ont in names(go_result_list)) {
    ego <- go_result_list[[ont]]
    if (nrow(ego) > 0) {
      p_dot <- dotplot(ego, showCategory = 15) + ggtitle(paste(title_prefix, "-", ont)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      p_bar <- barplot(ego, showCategory = 10) + ggtitle(paste(title_prefix, "-", ont, "Top 10")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      plots[[paste(ont, "dot")]] <- p_dot; plots[[paste(ont, "bar")]] <- p_bar
      print(p_dot); print(p_bar)
      ggsave(paste0(save_prefix, "_", ont, "_dotplot.pdf"), p_dot, width = 12, height = 8)
      ggsave(paste0(save_prefix, "_", ont, "_barplot.pdf"), p_bar, width = 10, height = 8)
    }
  }
  plots
}

extract_top_pathways <- function(go_results, ont_type, top_n = 10) {
  if (is.null(go_results) || is.null(go_results[[ont_type]])) return(data.frame())
  ego <- go_results[[ont_type]]
  if (nrow(ego) == 0) return(data.frame())
  ego@result %>% arrange(p.adjust) %>% head(top_n) %>%
    select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
    mutate(Gene_Ratio_Numeric = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
    }))
}

get_top_pathway_names <- function(go_results, ont = "BP", top_n = 5) {
  if (is.null(go_results) || is.null(go_results[[ont]]) || nrow(go_results[[ont]]) == 0) {
    return("No significant pathways found")
  }
  paste(go_results[[ont]]@result %>% arrange(p.adjust) %>% head(top_n) %>% pull(Description), collapse = "; ")
}

save_go_results <- function(go_results, file_prefix) {
  if (is.null(go_results)) return(invisible(NULL))
  for (ont in names(go_results)) {
    if (nrow(go_results[[ont]]) > 0) {
      write.csv(go_results[[ont]]@result, paste0(file_prefix, "_", ont, "_GO_results.csv"), row.names = FALSE)
    }
  }
}

#' GO enrichment for one ref-vs-test comparison. Reads the DE CSVs written by
#' run_cluster_vs_cluster_de() from the working directory.
run_go_for_comparison <- function(ref = "LEC1", test = "LEC3", padj_threshold = 0.05, fc_threshold = 0.25) {
  .load_go_packages()
  tag <- paste0(ref, "_vs_", test)
  cat("Loading", tag, "differential expression results...\n")
  all_res <- read.csv(paste0("ALL_", tag, "_DE_results.csv"))
  dis_res <- read.csv(paste0("DISEASED_", tag, "_DE_results.csv"))
  cat("ALL:", nrow(all_res), "genes; DISEASED:", nrow(dis_res), "genes\n")

  split_lists <- function(de, name) {
    sig <- de[!is.na(de$padj) & de$padj < padj_threshold & abs(de$log2FoldChange) > fc_threshold, ]
    out <- list(test_up = sig$gene[sig$log2FoldChange > 0], ref_up = sig$gene[sig$log2FoldChange < 0])
    cat(name, ": ", test, " up ", length(out$test_up), " | ", ref, " up ", length(out$ref_up), "\n", sep = "")
    out
  }
  all_lists <- split_lists(all_res, paste("ALL", tag))
  dis_lists <- split_lists(dis_res, paste("DISEASED", tag))

  cat("\n=== CONVERTING GENE SYMBOLS TO ENTREZ IDS ===\n")
  sets <- list()
  sets[[paste0("ALL_", test, "_UP")]]      <- list(genes = all_lists$test_up, universe = "all")
  sets[[paste0("ALL_", ref, "_UP")]]       <- list(genes = all_lists$ref_up,  universe = "all")
  sets[[paste0("DISEASED_", test, "_UP")]] <- list(genes = dis_lists$test_up, universe = "dis")
  sets[[paste0("DISEASED_", ref, "_UP")]]  <- list(genes = dis_lists$ref_up,  universe = "dis")
  universes <- list(all = convert_to_entrez(all_res$gene, "ALL universe"),
                    dis = convert_to_entrez(dis_res$gene, "DISEASED universe"))

  cat("\n=== RUNNING GO ENRICHMENT ANALYSIS ===\n")
  go <- list(); entrez <- list()
  for (nm in names(sets)) {
    entrez[[nm]] <- convert_to_entrez(sets[[nm]]$genes, nm)
    go[[nm]]     <- run_go_enrichment(entrez[[nm]], nm, universes[[sets[[nm]]$universe]])
  }

  cat("\n=== CREATING GO VISUALIZATION PLOTS ===\n")
  for (nm in names(go)) if (!is.null(go[[nm]])) create_go_plots(go[[nm]], gsub("_", " ", nm), nm)

  ## top BP pathway comparison across the four gene sets
  bp_combined <- data.frame()
  for (nm in names(go)) {
    paths <- extract_top_pathways(go[[nm]], "BP")
    if (nrow(paths) > 0) { paths$Comparison <- nm; bp_combined <- rbind(bp_combined, paths) }
  }
  if (nrow(bp_combined) > 0) {
    p <- bp_combined %>% group_by(Comparison) %>% slice_head(n = 8) %>%
      ggplot(aes(x = Gene_Ratio_Numeric, y = reorder(Description, Gene_Ratio_Numeric))) +
      geom_point(aes(color = -log10(p.adjust), size = Count)) +
      scale_color_gradient(low = "blue", high = "red", name = "-Log10(padj)") +
      scale_size_continuous(range = c(2, 8), name = "Gene Count") +
      facet_wrap(~ Comparison, scales = "free_y") +
      labs(title = paste("Top GO Biological Process Pathways -", ref, "vs", test),
           x = "Gene Ratio", y = "GO Terms") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 8),
            strip.text = element_text(face = "bold"))
    print(p)
    ggsave(paste0("GO_Pathway_Comparison_", tag, ".pdf"), p, width = 16, height = 12)
  }
  write.csv(bp_combined, paste0("GO_Biological_Process_Comparison_", tag, ".csv"), row.names = FALSE)

  ## pathway overlaps (BP) between ALL and DISEASED, and between the two clusters
  cat("\n=== PATHWAY OVERLAP ANALYSIS ===\n")
  sig_paths <- lapply(go, function(g) {
    if (is.null(g) || is.null(g$BP) || nrow(g$BP) == 0) return(character(0))
    g$BP@result$ID[g$BP@result$p.adjust < 0.05]
  })
  for (nm in names(sig_paths)) cat(nm, ":", length(sig_paths[[nm]]), "significant BP pathways\n")
  ov <- function(a, b) length(intersect(sig_paths[[a]], sig_paths[[b]]))
  cat(test, "shared pathways (ALL vs DISEASED):", ov(paste0("ALL_", test, "_UP"), paste0("DISEASED_", test, "_UP")), "\n")
  cat(ref,  "shared pathways (ALL vs DISEASED):", ov(paste0("ALL_", ref, "_UP"),  paste0("DISEASED_", ref, "_UP")),  "\n")
  cat("Shared pathways between clusters (ALL):",      ov(paste0("ALL_", ref, "_UP"),      paste0("ALL_", test, "_UP")),      "\n")
  cat("Shared pathways between clusters (DISEASED):", ov(paste0("DISEASED_", ref, "_UP"), paste0("DISEASED_", test, "_UP")), "\n")

  ## save everything
  for (nm in names(go)) save_go_results(go[[nm]], nm)
  summary_report <- data.frame(
    Comparison      = names(go),
    Gene_Count      = sapply(entrez, length),
    Top_BP_Pathways = sapply(go, get_top_pathway_names, ont = "BP")
  )
  write.csv(summary_report, paste0("GO_Analysis_Summary_Report_", tag, ".csv"), row.names = FALSE)
  gene_lists_summary <- data.frame(
    Comparison          = names(sets),
    Gene_Count          = sapply(sets, function(s) length(s$genes)),
    Converted_to_Entrez = sapply(entrez, length)
  )
  write.csv(gene_lists_summary, paste0("GO_Analysis_Gene_Lists_Summary_", tag, ".csv"), row.names = FALSE)

  cat("\n=== KEY PATHWAY FINDINGS (", tag, ") ===\n", sep = "")
  for (i in seq_len(nrow(summary_report))) {
    cat("\n", summary_report$Comparison[i], "(", summary_report$Gene_Count[i], "genes):\n")
    cat("  Top pathways:", substr(summary_report$Top_BP_Pathways[i], 1, 200), "\n")
  }
  cat("\n=== GO PATHWAY ANALYSIS FOR", tag, "COMPLETE ===\n")
  invisible(list(go = go, bp_combined = bp_combined, summary = summary_report))
}
