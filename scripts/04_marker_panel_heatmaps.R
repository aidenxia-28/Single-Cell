## EC / immune marker-panel heatmaps after each cleaning step
## Stage 04 of the GSE159354 EC / immune pipeline.
##
## Inputs:  GSE159354_ThreeGroup_EC-Immune_harmony.rds (stage 01),
##          GSE159354_FilteredClusters_EC-Immune_harmony.rds (stage 02),
##          GSE159354_Clean_EC-Immune_harmony.rds (stage 03)
## Outputs: <Label>_Expression_Heatmap.pdf and Enhanced_<Label>_Expression_Heatmap.pdf
##          for Combined, Filtered and Clean
## The heatmap shows row-z-scored expression of lymphatic / arterial / venous / core-EC /
## immune markers in every cell, cells ordered EC then immune, so contamination of
## either compartment is visible as a block of the wrong markers. The function is in
## R/marker_panel_heatmap.R; the notebook repeated the code once per object.
##
## Run from the repository root:  Rscript scripts/04_marker_panel_heatmaps.R

source("config.R")   # paths, shared helpers; switches to the results directory

objects <- c(
  Combined = "GSE159354_ThreeGroup_EC-Immune_harmony.rds",
  Filtered = "GSE159354_FilteredClusters_EC-Immune_harmony.rds",
  Clean    = "GSE159354_Clean_EC-Immune_harmony.rds"
)

loaded <- list()
for (label in names(objects)) {
  f <- objects[[label]]
  if (!file.exists(f)) { cat("Skipping", label, "-", f, "not found\n"); next }
  cat("\n#####", label, "#####\n")
  loaded[[label]] <- readRDS(f)
  plot_marker_panel_heatmap(loaded[[label]], label = label)
}

if (all(c("Filtered", "Clean") %in% names(loaded))) {
  report_cleaning_delta(loaded$Filtered, loaded$Clean, "BEFORE CLEANING", "AFTER CLEANING")
}
