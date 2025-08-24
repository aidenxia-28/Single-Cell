## LEC <-> immune ligand-receptor analysis with CellChat
## Stage 40 of the GSE159354 EC / immune pipeline (joins the two branches).
##
## Inputs:  LECs_subclustered.rds (stage 20), Immune_Renamed_0to20_Final.rds (stage 32)
## Outputs: Complete_* and Targeted_* CellChat objects, network / pathway PDFs, communication CSVs, Targeted_LEC_Signaling_Heatmaps.pdf, Targeted_LEC_Signaling_Matrix.pdf
## "Complete" uses all eight LEC subclusters; "Targeted" restricts the LEC side
## to the disease-enriched LEC0, LEC3 and LEC5. See R/cellchat_lec_immune.R.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

sce_LECs       <- readRDS("LECs_subclustered.rds")
immune_filtered <- readRDS("Immune_Renamed_0to20_Final.rds")

Idents(sce_LECs)        <- sce_LECs$LEC_subcluster
Idents(immune_filtered) <- immune_filtered$renamed_clusters

# 1) all LEC subclusters (LEC0-LEC7) vs all 21 immune clusters
complete_cellchat <- run_lec_immune_cellchat(sce_LECs, immune_filtered, lec_clusters = NULL,
                                             prefix = "Complete", n_pathways = 8, top_n = 25)

# 2) disease-enriched LEC0 / LEC3 / LEC5 as senders
targeted_cellchat <- run_lec_immune_cellchat(sce_LECs, immune_filtered, lec_clusters = c("LEC0", "LEC3", "LEC5"),
                                             prefix = "Targeted", n_pathways = 6, top_n = 30)
plot_lec_to_immune_heatmaps(targeted_cellchat, prefix = "Targeted", row_order = c("LEC3", "LEC5", "LEC0"))

cat("\nCellChat analysis complete.\n")
