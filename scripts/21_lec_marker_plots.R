## Lymphatic branch: violin and feature plots of candidate LEC-state genes
## Stage 21 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: LEC_Candidate_Gene_Plots.pdf
## The notebook drew each gene by hand; this is the same set of plots in a
## loop.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

sce_LECs <- readRDS("LECs_subclustered.rds")

# Violin plots are drawn on LEC1-LEC7 (LEC0 excluded, as in the original analysis);
# feature plots use all LECs.
sce_subset <- subset(sce_LECs, subset = LEC_subcluster %in% paste0("LEC", 1:7))

candidate_genes <- c("SPARCL1", "IL3RA", "POSTN", "KLF4", "TBX1", "EBF1", "COL15A1", "STC1", "HIF1A", "EDNRB", "FLT4")
candidate_genes <- candidate_genes[candidate_genes %in% rownames(sce_LECs)]

pdf("LEC_Candidate_Gene_Plots.pdf", width = 10, height = 7)
for (gene in candidate_genes) {
  print(
    VlnPlot(sce_subset, features = gene, group.by = "LEC_subcluster", pt.size = 0.1) +
      labs(title = paste(gene, "Expression Across LEC Clusters 1-7"),
           x = "LEC Clusters", y = "Expression Level") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
  )
  print(
    FeaturePlot(sce_LECs, features = gene, pt.size = 0.5, cols = c("lightgray", "blue"))
  )
}
dev.off()
cat("Saved LEC_Candidate_Gene_Plots.pdf for", length(candidate_genes), "genes\n")
