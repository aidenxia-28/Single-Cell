## Lymphatic branch: GO enrichment of the LEC1-vs-LEC3 and LEC1-vs-LEC6 gene lists
## Stage 27 of the GSE159354 EC / immune pipeline.
##
## Inputs:  ALL_/DISEASED_<ref>_vs_<test>_DE_results.csv (stage 26)
## Outputs: <set>_<BP|MF|CC>_GO_results.csv, <set>_*_dotplot.pdf / _barplot.pdf, GO_Pathway_Comparison_<ref>_vs_<test>.pdf, GO_*_Summary_*.csv
## clusterProfiler::enrichGO (org.Hs.eg.db), universe = all tested genes, BH <
## 0.05, q < 0.2. Add c("LEC1", "LEC5") to the list to run the third
## comparison.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

comparisons <- list(c("LEC1", "LEC3"), c("LEC1", "LEC6"))

go_results <- list()
for (cmp in comparisons) {
  tag <- paste0(cmp[1], "_vs_", cmp[2])
  cat("\n\n################ GO:", tag, "################\n")
  go_results[[tag]] <- run_go_for_comparison(ref = cmp[1], test = cmp[2])
}
saveRDS(go_results, "LEC_GO_enrichment_results.rds")
