## Lymphatic branch: donor-level pseudobulk DE between LEC subclusters (LEC1 as reference)
## Stage 26 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: ALL_/DISEASED_<ref>_vs_<test>_DE_results.csv, *_Volcano.pdf, TOP_50_*.csv, *_upregulated_genes.csv, <ref>_vs_<test>_Volcano_Summary_Stats.csv
## LEC1 (the CTRL-enriched, homeostatic subcluster) is compared against the
## disease-enriched LEC3 and LEC5 and against LEC6. Each comparison is run on
## all donors and on IPF + SSc-ILD donors only. See R/pseudobulk_de.R.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

sce_LECs <- readRDS("LECs_subclustered.rds")
cat("LEC subclusters:\n"); print(table(sce_LECs$LEC_subcluster, sce_LECs$group))

comparisons <- list(c("LEC1", "LEC3"), c("LEC1", "LEC6"), c("LEC1", "LEC5"))

de_results <- list()
for (cmp in comparisons) {
  tag <- paste0(cmp[1], "_vs_", cmp[2])
  cat("\n\n################", tag, "################\n")
  de_results[[tag]] <- run_cluster_vs_cluster_de(sce_LECs, ref = cmp[1], test = cmp[2])
}
saveRDS(de_results, "LEC_cluster_vs_cluster_DE_results.rds")
cat("\nAll cluster-vs-cluster comparisons complete.\n")
