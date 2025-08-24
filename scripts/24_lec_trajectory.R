## Lymphatic branch: Slingshot (and optional Monocle3) trajectories through the LEC subclusters
## Stage 24 of the GSE159354 EC / immune pipeline.
##
## Inputs:  LECs_subclustered.rds (stage 20)
## Outputs: LEC_Slingshot_Trajectories.rds, LEC_Slingshot_Pseudotime.csv, LEC_*Trajector*.pdf, LEC_Pseudotime_*.pdf, LEC_Gene_*.pdf, LEC_Trajectory_Gene_Analysis.csv (+ LEC_Monocle3_* when monocle3 is installed)
## Slingshot is run twice (automatic, and guided with LEC0 as start and
## LEC6/LEC7 as ends). Subclusters are then labelled Normal (LEC1, LEC6: CTRL-
## enriched), Proinflammatory (LEC0, LEC2, LEC3, LEC5: disease-enriched) or
## Intermediate (LEC4, LEC7) and PROX1 / EPAS1 / HIF1A / SPARCL1 are followed
## along pseudotime.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# =============================================================================
# LEC Trajectory Analysis - LEC0 through LEC7
# =============================================================================

# =============================================================================
# LEC Trajectory Analysis - Slingshot & Monocle3 (with installation handling)
# =============================================================================

library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(patchwork)

# =============================================================================
# MONOCLE3 IS OPTIONAL
# =============================================================================
# Install it with devtools::install_github("cole-trapnell-lab/monocle3") if you
# want the Monocle3 graph; Slingshot runs either way.
monocle3_available <- requireNamespace("monocle3", quietly = TRUE)
if (monocle3_available) {
  suppressPackageStartupMessages(library(monocle3))
  cat("Monocle3 loaded\n")
} else {
  cat("Monocle3 not available; proceeding with Slingshot only\n")
}

# =============================================================================
# LOAD LEC DATA
# =============================================================================

cat("Loading LEC subclustered object for trajectory analysis...\n")
sce_LECs <- readRDS("LECs_subclustered.rds")

cat("LEC clusters for trajectory analysis:\n")
print(table(sce_LECs$LEC_subcluster))

cat("Total cells:", ncol(sce_LECs), "\n")
cat("Groups:", paste(unique(sce_LECs$group), collapse = ", "), "\n")

# =============================================================================
# PREPARE DATA FOR TRAJECTORY ANALYSIS
# =============================================================================

# Set cluster identities
Idents(sce_LECs) <- sce_LECs$LEC_subcluster

# Join layers for Seurat v5 compatibility
if(packageVersion("Seurat") >= "5.0.0") {
  sce_LECs[["RNA"]] <- JoinLayers(sce_LECs[["RNA"]])
}

# Get cluster order
cluster_order <- paste0("LEC", 0:7)
available_clusters <- intersect(cluster_order, unique(sce_LECs$LEC_subcluster))

cat("Available clusters for trajectory:", paste(available_clusters, collapse = ", "), "\n")

# =============================================================================
# SLINGSHOT TRAJECTORY ANALYSIS
# =============================================================================

cat("\n=== SLINGSHOT TRAJECTORY ANALYSIS ===\n")

# Convert to SingleCellExperiment for slingshot
cat("Converting to SingleCellExperiment...\n")
sce_for_slingshot <- as.SingleCellExperiment(sce_LECs)

# Get UMAP coordinates
umap_coords <- Embeddings(sce_LECs, "umap")
colnames(umap_coords) <- c("UMAP1", "UMAP2")

# Add UMAP to SCE object
reducedDim(sce_for_slingshot, "UMAP") <- umap_coords

# Add cluster information
sce_for_slingshot$clusters <- sce_LECs$LEC_subcluster

# Run Slingshot trajectory analysis
cat("Running Slingshot trajectory inference...\n")

# Option 1: Let Slingshot automatically determine start/end clusters
sds_auto <- slingshot(sce_for_slingshot,
                      clusterLabels = "clusters",
                      reducedDim = "UMAP",
                      approx_points = 150)

cat("Slingshot found", length(slingCurves(sds_auto)), "trajectory/trajectories\n")

# Option 2: Specify potential start and end clusters based on biology
# You can modify these based on your biological knowledge
# For example, if LEC0 represents progenitor-like cells and LEC7 represents mature cells
sds_guided <- slingshot(sce_for_slingshot,
                        clusterLabels = "clusters",
                        reducedDim = "UMAP",
                        start.clus = "LEC0",  # Modify based on your biology
                        end.clus = c("LEC6", "LEC7"),  # Modify based on your biology
                        approx_points = 150)

cat("Guided Slingshot found", length(slingCurves(sds_guided)), "trajectory/trajectories\n")

# =============================================================================
# VISUALIZE SLINGSHOT TRAJECTORIES
# =============================================================================

cat("Creating Slingshot trajectory visualizations...\n")

# Define colors for clusters
cluster_colors <- brewer.pal(length(available_clusters), "Set2")
names(cluster_colors) <- available_clusters

# Create trajectory plot - automatic version
trajectory_data_auto <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Cluster = sce_LECs$LEC_subcluster,
  Group = sce_LECs$group
)

p_slingshot_auto <- ggplot(trajectory_data_auto, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Cluster), size = 0.8, alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal() +
  labs(title = "Slingshot Trajectories - Automatic Detection",
       subtitle = paste("Found", length(slingCurves(sds_auto)), "trajectories")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Add trajectory curves to automatic plot
for(i in seq_along(slingCurves(sds_auto))) {
  curve_data <- data.frame(slingCurves(sds_auto)[[i]]$s[slingCurves(sds_auto)[[i]]$ord, ])
  colnames(curve_data) <- c("UMAP1", "UMAP2")
  p_slingshot_auto <- p_slingshot_auto +
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              color = "black", size = 2, alpha = 0.8) +
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              color = rainbow(length(slingCurves(sds_auto)))[i], size = 1.5)
}

print(p_slingshot_auto)

# Create trajectory plot - guided version
p_slingshot_guided <- ggplot(trajectory_data_auto, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Cluster), size = 0.8, alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal() +
  labs(title = "Slingshot Trajectories - Guided (LEC0 start)",
       subtitle = paste("Found", length(slingCurves(sds_guided)), "trajectories")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Add trajectory curves to guided plot
for(i in seq_along(slingCurves(sds_guided))) {
  curve_data <- data.frame(slingCurves(sds_guided)[[i]]$s[slingCurves(sds_guided)[[i]]$ord, ])
  colnames(curve_data) <- c("UMAP1", "UMAP2")
  p_slingshot_guided <- p_slingshot_guided +
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              color = "black", size = 2, alpha = 0.8) +
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              color = rainbow(length(slingCurves(sds_guided)))[i], size = 1.5)
}

print(p_slingshot_guided)

# =============================================================================
# PSEUDOTIME ANALYSIS WITH SLINGSHOT
# =============================================================================

cat("Analyzing pseudotime with Slingshot...\n")

# Get pseudotime values (using guided version)
pseudotime_matrix <- slingPseudotime(sds_guided)

# Create pseudotime plots for each trajectory
if(ncol(pseudotime_matrix) > 0) {

  pseudotime_plots <- list()

  for(i in 1:ncol(pseudotime_matrix)) {
    trajectory_name <- paste("Trajectory", i)

    # Create data frame with pseudotime
    pt_data <- data.frame(
      UMAP1 = umap_coords[, 1],
      UMAP2 = umap_coords[, 2],
      Pseudotime = pseudotime_matrix[, i],
      Cluster = sce_LECs$LEC_subcluster,
      Group = sce_LECs$group
    )

    # Remove cells not assigned to this trajectory (NA pseudotime)
    pt_data <- pt_data[!is.na(pt_data$Pseudotime), ]

    if(nrow(pt_data) > 0) {
      p_pseudotime <- ggplot(pt_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = Pseudotime), size = 1.0, alpha = 0.8) +
        scale_color_viridis_c(name = "Pseudotime") +
        theme_minimal() +
        labs(title = paste("Pseudotime -", trajectory_name),
             subtitle = paste(nrow(pt_data), "cells assigned")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

      pseudotime_plots[[i]] <- p_pseudotime
      print(p_pseudotime)
    }
  }
}

# =============================================================================
# MONOCLE3 TRAJECTORY ANALYSIS (if available)
# =============================================================================

if(monocle3_available) {
  cat("\n=== MONOCLE3 TRAJECTORY ANALYSIS ===\n")

  # Convert Seurat to Monocle3 CDS object
  cat("Converting to Monocle3 CDS object...\n")

  # Extract necessary data
  expression_matrix <- GetAssayData(sce_LECs, slot = "counts", assay = "RNA")
  cell_metadata <- sce_LECs@meta.data
  gene_metadata <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names = rownames(expression_matrix)
  )

  # Create CDS object
  cds <- new_cell_data_set(expression_data = expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)

  # Add existing dimensionality reductions
  cds <- preprocess_cds(cds, num_dim = 50, verbose = FALSE)

  # Add UMAP coordinates from Seurat
  reducedDims(cds)$UMAP <- umap_coords

  # Add cluster information
  cds$clusters <- sce_LECs$LEC_subcluster
  cds$group <- sce_LECs$group

  # Learn trajectory graph
  cat("Learning trajectory graph with Monocle3...\n")
  cds <- learn_graph(cds, verbose = FALSE)

  # Order cells along trajectory
  # You can specify root cells based on biological knowledge
  # For example, if LEC0 represents the most progenitor-like state:
  root_cells <- colnames(cds)[cds$clusters == "LEC0"]

  if(length(root_cells) > 0) {
    # Use a subset of root cells if too many
    if(length(root_cells) > 50) {
      root_cells <- sample(root_cells, 50)
    }

    cds <- order_cells(cds, root_cells = root_cells, verbose = FALSE)

    cat("Trajectory ordering complete with", length(root_cells), "root cells from LEC0\n")
  } else {
    # If no LEC0 cells, let Monocle3 choose
    cds <- order_cells(cds, verbose = FALSE)
    cat("Trajectory ordering complete with automatic root selection\n")
  }

  # =============================================================================
  # VISUALIZE MONOCLE3 TRAJECTORIES
  # =============================================================================

  cat("Creating Monocle3 trajectory visualizations...\n")

  # Plot trajectory colored by clusters
  p_monocle_clusters <- plot_cells(cds,
                                  color_cells_by = "clusters",
                                  label_groups_by_cluster = FALSE,
                                  group_label_size = 4) +
    labs(title = "Monocle3 Trajectory - Colored by Clusters") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  print(p_monocle_clusters)

  # Plot trajectory colored by pseudotime
  p_monocle_pseudotime <- plot_cells(cds,
                                    color_cells_by = "pseudotime",
                                    label_groups_by_cluster = FALSE,
                                    group_label_size = 4) +
    labs(title = "Monocle3 Trajectory - Colored by Pseudotime") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  print(p_monocle_pseudotime)

  # Plot trajectory colored by disease group
  p_monocle_group <- plot_cells(cds,
                               color_cells_by = "group",
                               label_groups_by_cluster = FALSE,
                               group_label_size = 4) +
    labs(title = "Monocle3 Trajectory - Colored by Disease Group") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  print(p_monocle_group)

} else {
  cat("\n=== MONOCLE3 NOT AVAILABLE - SKIPPING MONOCLE3 ANALYSIS ===\n")
  cat("Continuing with Slingshot analysis only...\n")
}

# =============================================================================
# GENE EXPRESSION ALONG TRAJECTORIES
# =============================================================================

cat("Analyzing gene expression along trajectories...\n")

# Select key LEC genes for trajectory analysis
trajectory_genes <- c("PROX1", "LYVE1", "FLT4", "PDPN", "CCL21", "VEGFC",
                      "HIF1A", "EPAS1", "ICAM1", "CD74", "CCL2", "TEK")

if(monocle3_available && exists("cds")) {
  available_traj_genes <- intersect(trajectory_genes, rownames(cds))

  cat("Available genes for Monocle3 trajectory analysis:", paste(available_traj_genes, collapse = ", "), "\n")

  if(length(available_traj_genes) >= 4) {
    # Plot gene expression along pseudotime
    p_genes_pseudotime <- plot_cells(cds,
                                     genes = available_traj_genes[1:min(4, length(available_traj_genes))],
                                     show_trajectory_graph = FALSE,
                                     label_cell_groups = FALSE,
                                     label_leaves = FALSE)

    print(p_genes_pseudotime)

    # Plot genes in pseudotime order
    cds_subset <- cds[available_traj_genes, ]
    p_genes_time <- plot_genes_in_pseudotime(cds_subset[available_traj_genes[1:min(6, length(available_traj_genes))], ],
                                             color_cells_by = "clusters",
                                             min_expr = 0.1,
                                             ncol = 2)
    print(p_genes_time)
  }
} else {
  cat("Monocle3 not available - skipping gene expression trajectory analysis\n")
  cat("You can still analyze gene expression along Slingshot pseudotime manually\n")

  # Alternative: Simple correlation analysis with Slingshot pseudotime
  if(exists("pseudotime_matrix") && ncol(pseudotime_matrix) > 0) {
    available_traj_genes <- intersect(trajectory_genes, rownames(sce_LECs))

    if(length(available_traj_genes) > 0) {
      cat("Performing simple correlation analysis with Slingshot pseudotime...\n")

      # Get expression data
      expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[available_traj_genes, , drop = FALSE]

      # Calculate correlations with first trajectory pseudotime
      pt_values <- pseudotime_matrix[, 1]
      valid_cells <- !is.na(pt_values)

      if(sum(valid_cells) > 10) {
        gene_correlations <- apply(expr_data[, valid_cells, drop = FALSE], 1, function(gene_expr) {
          cor(gene_expr, pt_values[valid_cells], method = "spearman")
        })

        cat("Gene correlations with pseudotime (first trajectory):\n")
        print(sort(gene_correlations, decreasing = TRUE))

        # Create simple scatter plots for top correlated genes
        top_genes <- names(sort(abs(gene_correlations), decreasing = TRUE))[1:min(4, length(gene_correlations))]

        for(gene in top_genes) {
          gene_expr <- as.numeric(expr_data[gene, valid_cells])
          plot_data <- data.frame(
            Pseudotime = pt_values[valid_cells],
            Expression = gene_expr,
            Cluster = sce_LECs$LEC_subcluster[valid_cells]
          )

          p_gene_scatter <- ggplot(plot_data, aes(x = Pseudotime, y = Expression)) +
            geom_point(aes(color = Cluster), alpha = 0.6) +
            geom_smooth(method = "loess", color = "black", se = TRUE) +
            labs(title = paste("Gene Expression along Pseudotime:", gene),
                 x = "Slingshot Pseudotime", y = "Expression") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))

          print(p_gene_scatter)
        }
      }
    }
  }
}

# =============================================================================
# TRAJECTORY COMPARISON BY DISEASE GROUP
# =============================================================================

cat("Comparing trajectories across disease groups...\n")

# Split by disease group and compare
disease_groups <- unique(sce_LECs$group)

# Create pseudotime distribution plots by disease group
if(exists("pseudotime_matrix") && ncol(pseudotime_matrix) > 0) {

  pseudotime_comparison <- data.frame(
    Pseudotime = as.vector(pseudotime_matrix[, 1]),  # Use first trajectory
    Group = rep(sce_LECs$group, ncol(pseudotime_matrix))[1:length(as.vector(pseudotime_matrix[, 1]))],
    Cluster = rep(sce_LECs$LEC_subcluster, ncol(pseudotime_matrix))[1:length(as.vector(pseudotime_matrix[, 1]))]
  )

  # Remove NA values
  pseudotime_comparison <- pseudotime_comparison[!is.na(pseudotime_comparison$Pseudotime), ]

  if(nrow(pseudotime_comparison) > 0) {

    # Pseudotime distribution by disease group
    p_pseudotime_dist <- ggplot(pseudotime_comparison, aes(x = Pseudotime, fill = Group)) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(values = c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")) +
      labs(title = "Pseudotime Distribution by Disease Group",
           x = "Pseudotime", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    print(p_pseudotime_dist)

    # Boxplot comparison
    p_pseudotime_box <- ggplot(pseudotime_comparison, aes(x = Group, y = Pseudotime, fill = Group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
      scale_fill_manual(values = c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")) +
      labs(title = "Pseudotime by Disease Group",
           x = "Disease Group", y = "Pseudotime") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))

    print(p_pseudotime_box)
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== TRAJECTORY ANALYSIS SUMMARY ===\n")

cat("Total LECs analyzed:", ncol(sce_LECs), "\n")
cat("Clusters included:", paste(available_clusters, collapse = ", "), "\n")

if(exists("sds_guided")) {
  cat("Slingshot trajectories found:", length(slingCurves(sds_guided)), "\n")

  # Show trajectory assignments
  for(i in 1:length(slingCurves(sds_guided))) {
    assigned_cells <- sum(!is.na(pseudotime_matrix[, i]))
    cat("  Trajectory", i, ":", assigned_cells, "cells assigned\n")
  }
}

if(exists("cds") && monocle3_available) {
  # Monocle3 summary
  pseudotime_range <- range(pseudotime(cds), na.rm = TRUE)
  assigned_cells_monocle <- sum(!is.na(pseudotime(cds)))
  cat("Monocle3 pseudotime range:", round(pseudotime_range[1], 2), "-", round(pseudotime_range[2], 2), "\n")
  cat("Monocle3 cells with pseudotime:", assigned_cells_monocle, "\n")
} else {
  cat("Monocle3 analysis not available\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING TRAJECTORY ANALYSIS RESULTS ===\n")

# Save plots
if(exists("p_slingshot_auto")) {
  ggsave("LEC_Slingshot_Trajectories_Auto.pdf", p_slingshot_auto, width = 10, height = 8)
}
if(exists("p_slingshot_guided")) {
  ggsave("LEC_Slingshot_Trajectories_Guided.pdf", p_slingshot_guided, width = 10, height = 8)
}
if(exists("p_monocle_clusters")) {
  ggsave("LEC_Monocle3_Clusters.pdf", p_monocle_clusters, width = 10, height = 8)
}
if(exists("p_monocle_pseudotime")) {
  ggsave("LEC_Monocle3_Pseudotime.pdf", p_monocle_pseudotime, width = 10, height = 8)
}
if(exists("p_monocle_group")) {
  ggsave("LEC_Monocle3_DiseaseGroup.pdf", p_monocle_group, width = 10, height = 8)
}
if(exists("p_pseudotime_dist")) {
  ggsave("LEC_Pseudotime_Distribution_by_Group.pdf", p_pseudotime_dist, width = 10, height = 6)
}
if(exists("p_pseudotime_box")) {
  ggsave("LEC_Pseudotime_Boxplot_by_Group.pdf", p_pseudotime_box, width = 8, height = 6)
}

# Save trajectory objects
if(exists("sds_guided")) {
  saveRDS(sds_guided, "LEC_Slingshot_Trajectories.rds")
}
if(exists("cds") && monocle3_available) {
  saveRDS(cds, "LEC_Monocle3_CDS.rds")
}

# Save pseudotime data
if(exists("pseudotime_matrix")) {
  slingshot_pseudotime <- data.frame(
    Cell_ID = rownames(pseudotime_matrix),
    pseudotime_matrix,
    Cluster = sce_LECs$LEC_subcluster,
    Group = sce_LECs$group
  )
  write.csv(slingshot_pseudotime, "LEC_Slingshot_Pseudotime.csv", row.names = FALSE)
}

if(exists("cds") && monocle3_available) {
  monocle_pseudotime <- data.frame(
    Cell_ID = colnames(cds),
    Pseudotime = pseudotime(cds),
    Cluster = cds$clusters,
    Group = cds$group
  )
  write.csv(monocle_pseudotime, "LEC_Monocle3_Pseudotime.csv", row.names = FALSE)
}

cat("Files saved:\n")
if(exists("p_slingshot_auto")) {
  cat("- Slingshot trajectory plots (PDF)\n")
}
if(exists("p_monocle_clusters")) {
  cat("- Monocle3 trajectory plots (PDF)\n")
}
cat("- LEC_Slingshot_Trajectories.rds (Slingshot object)\n")
if(exists("cds") && monocle3_available) {
  cat("- LEC_Monocle3_CDS.rds (Monocle3 CDS object)\n")
  cat("- LEC_Monocle3_Pseudotime.csv\n")
}
cat("- LEC_Slingshot_Pseudotime.csv\n")

cat("\n=== LEC TRAJECTORY ANALYSIS COMPLETE ===\n")

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("TRAJECTORY ANALYSIS RESULTS:\n\n")

cat("1. SLINGSHOT RESULTS:\n")
cat("   - Identifies branching trajectories through UMAP space\n")
cat("   - Each trajectory represents a potential differentiation path\n")
cat("   - Pseudotime values show progression along each path\n\n")

cat("2. MONOCLE3 RESULTS:\n")
cat("   - Learns a trajectory graph connecting cell states\n")
cat("   - Pseudotime represents inferred developmental progression\n")
cat("   - Can identify branch points and terminal states\n\n")

cat("3. BIOLOGICAL INTERPRETATION FOR LECS:\n")
cat("   - Early pseudotime: Progenitor-like or quiescent LECs\n")
cat("   - Late pseudotime: Activated, specialized, or mature LECs\n")
cat("   - Branch points: Decision points between different LEC fates\n")
cat("   - Disease differences: How pathology affects LEC trajectories\n\n")

cat("4. NEXT STEPS:\n")
cat("   - Identify genes that change along pseudotime\n")
cat("   - Test for differential trajectory usage between disease groups\n")
cat("   - Validate trajectory predictions with known LEC biology\n")
cat("   - Consider spatial context if location data available\n\n")


# =============================================================================
# LOESS TREND OF SELECTED GENES ALONG THE FIRST SLINGSHOT TRAJECTORY
# =============================================================================

# Define genes to plot
genes_to_plot <- c("PDPN", "EPAS1", "HIF1A", "SPARCL1", "PROX1", "FLT4")

# Get pseudotime matrix from slingshot object
pseudotime_matrix <- slingPseudotime(sds_guided)

# Get pseudotime values for first trajectory
pt_values <- pseudotime_matrix[, 1]
valid_cells <- !is.na(pt_values)

# Create plots for each gene
plot_list <- list()

for(gene in genes_to_plot) {

  # Check if gene exists in the data
  if(gene %in% rownames(sce_LECs)) {

    # Get gene expression data
    gene_expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[gene, , drop = FALSE]
    gene_expr <- as.numeric(gene_expr_data[1, valid_cells])

    # Create plot data
    plot_data <- data.frame(
      Pseudotime = pt_values[valid_cells],
      Expression = gene_expr
    )

    # Create the plot with trajectory line only
    p_gene <- ggplot(plot_data, aes(x = Pseudotime, y = Expression)) +
      geom_smooth(method = "loess", color = "black", se = TRUE, size = 1.2) +
      labs(title = paste("Gene Expression along Pseudotime:", gene),
           x = "Slingshot Pseudotime", y = "Expression") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    plot_list[[gene]] <- p_gene
    print(p_gene)

  } else {
    cat("Gene", gene, "not found in dataset\n")
  }
}

pdf("LEC_Gene_Loess_Along_Pseudotime.pdf", width = 8, height = 6)
for (p in plot_list) print(p)
dev.off()


# =============================================================================
# SPECIFIC GENE TRAJECTORY ANALYSIS - PROX1, EPAS1, HIF1A, SPARCL1
# Normal vs Proinflammatory LEC Differentiation
# =============================================================================

library(reshape2)

# =============================================================================
# DEFINE TARGET GENES AND SETUP
# =============================================================================

# Target genes for trajectory analysis
target_genes <- c("PROX1", "EPAS1", "HIF1A", "SPARCL1")

# Check which genes are available in the dataset
available_genes <- intersect(target_genes, rownames(sce_LECs))
missing_genes <- setdiff(target_genes, available_genes)

cat("Target genes analysis:\n")
cat("Available genes:", paste(available_genes, collapse = ", "), "\n")
if(length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")

  # Check for alternative gene names
  potential_alternatives <- c("PROX-1" = "PROX1",
                             "SPARCl1" = "SPARCL1",
                             "SPARC-like-1" = "SPARCL1")

  for(alt_name in names(potential_alternatives)) {
    if(alt_name %in% rownames(sce_LECs)) {
      cat("Found alternative:", alt_name, "-> using", potential_alternatives[alt_name], "\n")
      available_genes <- c(available_genes, alt_name)
    }
  }
}

if(length(available_genes) == 0) {
  stop("None of the target genes found in the dataset!")
}

# =============================================================================
# CLASSIFY TRAJECTORIES AS NORMAL VS PROINFLAMMATORY
# =============================================================================

# Define proinflammatory vs normal LEC characteristics
# Based on cluster distribution across disease conditions:
# CTRL-enriched: LEC1 (527 vs 104+83), LEC6 (100 vs 29+44)
# Disease-enriched: LEC0 (289 vs 516+222), LEC2 (106 vs 322+216),
#                   LEC3 (21 vs 269+139), LEC5 (9 vs 151+50)

classify_lec_state <- function(cluster_name, disease_group) {
  # Classification based on cluster distribution patterns

  if(cluster_name == "LEC1") {
    # LEC1: Highly enriched in CTRL (527) vs disease (104 IPF + 83 SSc-ILD)
    return("Normal")
  } else if(cluster_name == "LEC6") {
    # LEC6: More abundant in CTRL (100) vs disease (29 IPF + 44 SSc-ILD)
    return("Normal")
  } else if(cluster_name %in% c("LEC0", "LEC2", "LEC3", "LEC5")) {
    # Disease-enriched clusters:
    # LEC0: 289 CTRL vs 516 IPF + 222 SSc-ILD
    # LEC2: 106 CTRL vs 322 IPF + 216 SSc-ILD
    # LEC3: 21 CTRL vs 269 IPF + 139 SSc-ILD (most disease-specific)
    # LEC5: 9 CTRL vs 151 IPF + 50 SSc-ILD (very disease-specific)
    return("Proinflammatory")
  } else if(cluster_name %in% c("LEC4", "LEC7")) {
    # LEC4: Mixed distribution, LEC7: Very few cells
    return("Intermediate")
  } else {
    return("Intermediate")
  }
}

# Add LEC state classification to metadata
sce_LECs$LEC_state <- mapply(classify_lec_state,
                            sce_LECs$LEC_subcluster,
                            sce_LECs$group)

cat("LEC state classification based on cluster distribution:\n")
cat("Normal LECs (CTRL-enriched): LEC1, LEC6\n")
cat("Proinflammatory LECs (Disease-enriched): LEC0, LEC2, LEC3, LEC5\n")
cat("Intermediate LECs: LEC4, LEC7\n\n")

cat("LEC state classification summary:\n")
print(table(sce_LECs$LEC_state, sce_LECs$group))

cat("\nCluster-specific distribution validation:\n")
cluster_state_table <- table(sce_LECs$LEC_subcluster, sce_LECs$LEC_state)
print(cluster_state_table)

# =============================================================================
# GENE EXPRESSION ALONG SLINGSHOT TRAJECTORIES
# =============================================================================

if(exists("pseudotime_matrix") && ncol(pseudotime_matrix) > 0) {

  cat("Creating gene trajectory plots with Slingshot pseudotime...\n")

  # Get expression data for target genes
  expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[available_genes, , drop = FALSE]

  # Create plots for each trajectory
  slingshot_gene_plots <- list()

  for(traj_idx in 1:ncol(pseudotime_matrix)) {

    trajectory_name <- paste("Trajectory", traj_idx)
    pt_values <- pseudotime_matrix[, traj_idx]
    valid_cells <- !is.na(pt_values)

    if(sum(valid_cells) > 10) {

      # Create data frame for this trajectory
      traj_data <- data.frame(
        Cell_ID = colnames(sce_LECs)[valid_cells],
        Pseudotime = pt_values[valid_cells],
        Cluster = sce_LECs$LEC_subcluster[valid_cells],
        Group = sce_LECs$group[valid_cells],
        LEC_state = sce_LECs$LEC_state[valid_cells],
        stringsAsFactors = FALSE
      )

      # Add gene expression data
      for(gene in available_genes) {
        traj_data[[gene]] <- as.numeric(expr_data[gene, valid_cells])
      }

      # Melt data for plotting
      traj_data_long <- melt(traj_data,
                            id.vars = c("Cell_ID", "Pseudotime", "Cluster", "Group", "LEC_state"),
                            measure.vars = available_genes,
                            variable.name = "Gene",
                            value.name = "Expression")

      # Create trajectory plot colored by LEC state
      p_trajectory_state <- ggplot(traj_data_long, aes(x = Pseudotime, y = Expression)) +
        geom_point(aes(color = LEC_state), alpha = 0.6, size = 0.8) +
        geom_smooth(aes(color = LEC_state), method = "loess", se = TRUE, alpha = 0.3) +
        facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
        scale_color_manual(values = c("Normal" = "#2CA02C",
                                     "Proinflammatory" = "#FF4500",
                                     "Intermediate" = "#FFB347")) +
        labs(title = paste("Gene Expression Along", trajectory_name),
             subtitle = "Normal vs Proinflammatory LEC Differentiation",
             x = "Pseudotime",
             y = "Expression Level",
             color = "LEC State") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text = element_text(face = "bold"),
              panel.grid.minor = element_blank())

      print(p_trajectory_state)
      slingshot_gene_plots[[traj_idx]] <- p_trajectory_state

      # Create disease group comparison
      p_trajectory_disease <- ggplot(traj_data_long, aes(x = Pseudotime, y = Expression)) +
        geom_point(aes(color = Group), alpha = 0.6, size = 0.8) +
        geom_smooth(aes(color = Group), method = "loess", se = TRUE, alpha = 0.3) +
        facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
        scale_color_manual(values = c("CTRL" = "#2CA02C",
                                     "IPF" = "#FF7F00",
                                     "SSc-ILD" = "#9467BD")) +
        labs(title = paste("Gene Expression Along", trajectory_name),
             subtitle = "Comparison by Disease Group",
             x = "Pseudotime",
             y = "Expression Level",
             color = "Disease Group") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text = element_text(face = "bold"),
              panel.grid.minor = element_blank())

      print(p_trajectory_disease)

      # Statistical analysis of gene expression differences
      cat(paste("\n=== Statistical Analysis for", trajectory_name, "===\n"))

      for(gene in available_genes) {
        gene_data <- traj_data[, c("Pseudotime", "LEC_state", gene)]
        colnames(gene_data)[3] <- "Expression"

        # Correlation analysis
        normal_cells <- gene_data$LEC_state == "Normal"
        proinflam_cells <- gene_data$LEC_state == "Proinflammatory"

        if(sum(normal_cells) > 5 && sum(proinflam_cells) > 5) {
          cor_normal <- cor(gene_data$Pseudotime[normal_cells],
                           gene_data$Expression[normal_cells],
                           method = "spearman")
          cor_proinflam <- cor(gene_data$Pseudotime[proinflam_cells],
                              gene_data$Expression[proinflam_cells],
                              method = "spearman")

          cat(paste(gene, "- Normal LEC correlation with pseudotime:", round(cor_normal, 3), "\n"))
          cat(paste(gene, "- Proinflammatory LEC correlation with pseudotime:", round(cor_proinflam, 3), "\n"))

          # T-test for expression differences
          t_test_result <- t.test(gene_data$Expression[normal_cells],
                                 gene_data$Expression[proinflam_cells])
          cat(paste(gene, "- Mean expression difference p-value:", format(t_test_result$p.value, scientific = TRUE), "\n"))
        }
      }
    }
  }
}

# =============================================================================
# GENE EXPRESSION ALONG MONOCLE3 TRAJECTORIES
# =============================================================================

if(exists("cds") && monocle3_available) {

  cat("Creating gene trajectory plots with Monocle3 pseudotime...\n")

  # Check which genes are available in the CDS object
  available_genes_monocle <- intersect(available_genes, rownames(cds))

  if(length(available_genes_monocle) > 0) {

    # Add LEC state to CDS object
    cds$LEC_state <- sce_LECs$LEC_state[match(colnames(cds), colnames(sce_LECs))]

    # Plot genes along pseudotime colored by LEC state
    p_monocle_genes_state <- plot_cells(cds,
                                       genes = available_genes_monocle,
                                       color_cells_by = "LEC_state",
                                       show_trajectory_graph = TRUE,
                                       label_cell_groups = FALSE,
                                       label_leaves = FALSE,
                                       cell_size = 0.8) +
      scale_color_manual(values = c("Normal" = "#2CA02C",
                                   "Proinflammatory" = "#FF4500",
                                   "Intermediate" = "#FFB347")) +
      labs(title = "Gene Expression - Monocle3 Trajectory",
           subtitle = "Colored by LEC State (Normal vs Proinflammatory)",
           color = "LEC State") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))

    print(p_monocle_genes_state)

    # Plot genes in pseudotime order
    if(length(available_genes_monocle) >= 2) {
      cds_subset <- cds[available_genes_monocle, ]

      p_genes_pseudotime_state <- plot_genes_in_pseudotime(cds_subset,
                                                          color_cells_by = "LEC_state",
                                                          min_expr = 0.1,
                                                          ncol = 2) +
        scale_color_manual(values = c("Normal" = "#2CA02C",
                                     "Proinflammatory" = "#FF4500",
                                     "Intermediate" = "#FFB347")) +
        labs(title = "Gene Expression Along Pseudotime - Monocle3",
             subtitle = "Normal vs Proinflammatory LEC States",
             color = "LEC State") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

      print(p_genes_pseudotime_state)
    }

    # Create detailed expression analysis
    monocle_data <- data.frame(
      Cell_ID = colnames(cds),
      Pseudotime = pseudotime(cds),
      Cluster = cds$clusters,
      Group = cds$group,
      LEC_state = cds$LEC_state
    )

    # Add gene expression
    expr_matrix <- exprs(cds)[available_genes_monocle, , drop = FALSE]
    for(gene in available_genes_monocle) {
      monocle_data[[gene]] <- as.numeric(expr_matrix[gene, ])
    }

    # Remove cells with NA pseudotime
    monocle_data <- monocle_data[!is.na(monocle_data$Pseudotime), ]

    if(nrow(monocle_data) > 0) {
      # Melt for plotting
      monocle_data_long <- melt(monocle_data,
                               id.vars = c("Cell_ID", "Pseudotime", "Cluster", "Group", "LEC_state"),
                               measure.vars = available_genes_monocle,
                               variable.name = "Gene",
                               value.name = "Expression")

      # Create comprehensive trajectory plot
      p_monocle_comprehensive <- ggplot(monocle_data_long, aes(x = Pseudotime, y = Expression)) +
        geom_point(aes(color = LEC_state), alpha = 0.6, size = 0.8) +
        geom_smooth(aes(color = LEC_state), method = "loess", se = TRUE, alpha = 0.3) +
        facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
        scale_color_manual(values = c("Normal" = "#2CA02C",
                                     "Proinflammatory" = "#FF4500",
                                     "Intermediate" = "#FFB347")) +
        labs(title = "Gene Expression Along Monocle3 Pseudotime",
             subtitle = "Normal vs Proinflammatory LEC Differentiation",
             x = "Monocle3 Pseudotime",
             y = "Expression Level",
             color = "LEC State") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text = element_text(face = "bold"),
              panel.grid.minor = element_blank())

      print(p_monocle_comprehensive)
    }
  }
}

# =============================================================================
# TRAJECTORY BRANCHING ANALYSIS
# =============================================================================

cat("Analyzing trajectory branching patterns...\n")

# Analyze which clusters contribute to normal vs proinflammatory trajectories
if(exists("pseudotime_matrix")) {

  trajectory_analysis <- data.frame(
    Cell_ID = rownames(pseudotime_matrix),
    Cluster = sce_LECs$LEC_subcluster,
    Group = sce_LECs$group,
    LEC_state = sce_LECs$LEC_state
  )

  # Add pseudotime for each trajectory
  for(i in 1:ncol(pseudotime_matrix)) {
    trajectory_analysis[[paste0("Trajectory_", i, "_PT")]] <- pseudotime_matrix[, i]
  }

  # Analyze trajectory usage by LEC state
  cat("Trajectory usage by LEC state:\n")

  for(i in 1:ncol(pseudotime_matrix)) {
    traj_col <- paste0("Trajectory_", i, "_PT")
    assigned_cells <- !is.na(trajectory_analysis[[traj_col]])

    if(sum(assigned_cells) > 0) {
      traj_states <- table(trajectory_analysis$LEC_state[assigned_cells])
      cat(paste("Trajectory", i, ":\n"))
      print(traj_states)
      cat("\n")
    }
  }
}

# =============================================================================
# COMPARATIVE HEATMAP
# =============================================================================

cat("Creating comparative gene expression heatmap...\n")

if(exists("pseudotime_matrix") && ncol(pseudotime_matrix) > 0) {

  # Use first trajectory for heatmap
  pt_values <- pseudotime_matrix[, 1]
  valid_cells <- !is.na(pt_values)

  if(sum(valid_cells) > 20) {
    # Bin pseudotime for heatmap
    n_bins <- 20
    pt_bins <- cut(pt_values[valid_cells], breaks = n_bins, labels = FALSE)

    # Create average expression by pseudotime bin and LEC state
    heatmap_data <- data.frame(
      Pseudotime_bin = pt_bins,
      LEC_state = sce_LECs$LEC_state[valid_cells]
    )

    # Add gene expression
    expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[available_genes, valid_cells, drop = FALSE]
    for(gene in available_genes) {
      heatmap_data[[gene]] <- as.numeric(expr_data[gene, ])
    }

    # Calculate mean expression by bin and state
    heatmap_summary <- heatmap_data %>%
      group_by(Pseudotime_bin, LEC_state) %>%
      summarise_at(available_genes, mean, na.rm = TRUE) %>%
      ungroup()

    # Create heatmap data
    heatmap_long <- melt(heatmap_summary,
                        id.vars = c("Pseudotime_bin", "LEC_state"),
                        measure.vars = available_genes,
                        variable.name = "Gene",
                        value.name = "Mean_Expression")

    # Filter out combinations with too few cells
    heatmap_long <- heatmap_long[!is.na(heatmap_long$Mean_Expression), ]

    if(nrow(heatmap_long) > 0) {
      p_heatmap <- ggplot(heatmap_long, aes(x = Pseudotime_bin, y = Gene)) +
        geom_tile(aes(fill = Mean_Expression)) +
        facet_wrap(~ LEC_state, ncol = 3) +
        scale_fill_viridis_c(name = "Mean\nExpression") +
        labs(title = "Gene Expression Along Pseudotime",
             subtitle = "Normal vs Proinflammatory LEC States",
             x = "Pseudotime Bin",
             y = "Gene") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(face = "bold"))

      print(p_heatmap)
    }
  }
}

# =============================================================================
# SAVE GENE-SPECIFIC TRAJECTORY RESULTS
# =============================================================================

cat("\n=== SAVING GENE-SPECIFIC TRAJECTORY RESULTS ===\n")

# Save gene trajectory plots
if(exists("slingshot_gene_plots")) {
  for(i in seq_along(slingshot_gene_plots)) {
    filename <- paste0("LEC_Gene_Trajectory_", i, "_Slingshot.pdf")
    ggsave(filename, slingshot_gene_plots[[i]], width = 12, height = 10)
    cat("Saved:", filename, "\n")
  }
}

if(exists("p_monocle_comprehensive")) {
  ggsave("LEC_Gene_Trajectory_Monocle3_Comprehensive.pdf", p_monocle_comprehensive, width = 12, height = 10)
  cat("Saved: LEC_Gene_Trajectory_Monocle3_Comprehensive.pdf\n")
}

if(exists("p_heatmap")) {
  ggsave("LEC_Gene_Expression_Heatmap_by_State.pdf", p_heatmap, width = 12, height = 8)
  cat("Saved: LEC_Gene_Expression_Heatmap_by_State.pdf\n")
}

# Save trajectory analysis data with gene expression
if(exists("trajectory_analysis")) {
  trajectory_gene_data <- trajectory_analysis

  # Add gene expression data
  expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[available_genes, , drop = FALSE]
  for(gene in available_genes) {
    trajectory_gene_data[[paste0(gene, "_Expression")]] <- as.numeric(expr_data[gene, ])
  }

  write.csv(trajectory_gene_data, "LEC_Trajectory_Gene_Analysis.csv", row.names = FALSE)
  cat("Saved: LEC_Trajectory_Gene_Analysis.csv\n")
}

cat("\n=== GENE-SPECIFIC TRAJECTORY ANALYSIS COMPLETE ===\n")

# =============================================================================
# INTERPRETATION SUMMARY
# =============================================================================

cat("\n=== INTERPRETATION GUIDE FOR GENE TRAJECTORIES ===\n")

cat("GENE FUNCTIONS IN LEC BIOLOGY:\n\n")

cat("1. PROX1 (Prospero Homeobox 1):\n")
cat("   - Master regulator of lymphatic development\n")
cat("   - Essential for lymphatic identity\n")
cat("   - Expected: High in normal LECs, potentially altered in disease\n\n")

cat("2. EPAS1 (HIF-2α):\n")
cat("   - Hypoxia-inducible factor 2 alpha\n")
cat("   - Regulates angiogenesis and metabolism\n")
cat("   - Expected: Upregulated in hypoxic/inflammatory conditions\n\n")

cat("3. HIF1A (HIF-1α):\n")
cat("   - Hypoxia-inducible factor 1 alpha\n")
cat("   - Key hypoxia response regulator\n")
cat("   - Expected: Elevated in proinflammatory states\n\n")

cat("4. SPARCL1 (SPARC-like 1/Hevin):\n")
cat("   - Extracellular matrix protein\n")
cat("   - Involved in cell adhesion and ECM remodeling\n")
cat("   - Expected: Altered in fibrotic conditions\n\n")

cat("EXPECTED TRAJECTORY PATTERNS:\n")
cat("- Normal LECs (LEC1, LEC6): Stable PROX1, low HIF1A/EPAS1, homeostatic SPARCL1\n")
cat("- Proinflammatory LECs (LEC0, LEC2, LEC3, LEC5): Variable PROX1, high HIF1A/EPAS1, altered SPARCL1\n")
cat("- LEC3 & LEC5: Most disease-specific clusters, expect strongest inflammatory signature\n")
cat("- Disease progression: Shift from LEC1/LEC6 dominance to LEC0/LEC2/LEC3/LEC5 patterns\n\n")

cat("CLUSTER-SPECIFIC INSIGHTS:\n")
cat("- LEC1: Homeostatic LECs (527 CTRL vs 187 total disease)\n")
cat("- LEC6: Another normal state (100 CTRL vs 73 total disease)\n")
cat("- LEC3: Highly disease-specific (21 CTRL vs 408 total disease)\n")
cat("- LEC5: Very disease-specific (9 CTRL vs 201 total disease)\n")
cat("- LEC0: Disease-associated but also present in CTRL\n")
cat("- LEC2: Disease-enriched with some CTRL presence\n\n")


# =============================================================================
# TRAJECTORY PATH VISUALIZATION WITH GENE EXPRESSION
# Shows actual trajectory curves colored by gene expression levels
# =============================================================================


# =============================================================================
# SETUP - ENSURE PREVIOUS ANALYSIS IS LOADED
# =============================================================================

# Target genes for trajectory path analysis
target_genes <- c("PROX1", "EPAS1", "HIF1A", "SPARCL1")

# Check which genes are available in the dataset
available_genes <- intersect(target_genes, rownames(sce_LECs))

cat("Creating trajectory path plots for genes:", paste(available_genes, collapse = ", "), "\n")

# Get UMAP coordinates
umap_coords <- Embeddings(sce_LECs, "umap")
colnames(umap_coords) <- c("UMAP1", "UMAP2")

# =============================================================================
# TRAJECTORY PATH PLOTS - GENE EXPRESSION OVERLAY
# =============================================================================

if(exists("sds_guided") && length(available_genes) > 0) {

  cat("Creating trajectory path plots with gene expression overlay...\n")

  # Get expression data
  expr_data <- GetAssayData(sce_LECs, slot = "data", assay = "RNA")[available_genes, , drop = FALSE]

  # Create base trajectory data
  trajectory_base_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Cluster = sce_LECs$LEC_subcluster,
    Group = sce_LECs$group,
    LEC_state = sce_LECs$LEC_state
  )

  # Add gene expression to base data
  for(gene in available_genes) {
    trajectory_base_data[[gene]] <- as.numeric(expr_data[gene, ])
  }

  # Create plots for each gene
  gene_trajectory_plots <- list()

  for(gene in available_genes) {

    cat(paste("Creating trajectory plot for", gene, "...\n"))

    # Create base plot with gene expression
    p_base <- ggplot(trajectory_base_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes_string(color = gene), size = 1.2, alpha = 0.8) +
      scale_color_viridis_c(name = paste(gene, "\nExpression"), option = "plasma") +
      labs(title = paste("LEC Trajectory Paths -", gene, "Expression"),
           subtitle = "Trajectory curves overlaid on gene expression",
           x = "UMAP1", y = "UMAP2") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            legend.title = element_text(size = 11),
            panel.grid.minor = element_blank())

    # Add trajectory curves
    for(i in seq_along(slingCurves(sds_guided))) {
      curve_data <- data.frame(slingCurves(sds_guided)[[i]]$s[slingCurves(sds_guided)[[i]]$ord, ])
      colnames(curve_data) <- c("UMAP1", "UMAP2")

      # Add thick black outline for visibility
      p_base <- p_base +
        geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                  color = "black", size = 3, alpha = 0.8) +
        geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                  color = "white", size = 2, alpha = 0.9)
    }

    gene_trajectory_plots[[gene]] <- p_base
    print(p_base)
  }

  # Create combined plot showing all genes
  if(length(available_genes) >= 2) {
    combined_gene_plot <- wrap_plots(gene_trajectory_plots, ncol = 2)
    print(combined_gene_plot)

    # Save combined plot
    ggsave("LEC_All_Genes_Trajectory_Paths.pdf", combined_gene_plot,
           width = 16, height = 12)
  }
}

# =============================================================================
# TRAJECTORY PATH WITH LEC STATE OVERLAY
# =============================================================================

if(exists("sds_guided")) {

  cat("Creating trajectory paths colored by LEC state...\n")

  # Define colors for LEC states
  state_colors <- c("Normal" = "#2CA02C",
                   "Proinflammatory" = "#FF4500",
                   "Intermediate" = "#FFB347")

  # Create LEC state trajectory plot
  p_lec_state <- ggplot(trajectory_base_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = LEC_state), size = 1.2, alpha = 0.8) +
    scale_color_manual(values = state_colors, name = "LEC State") +
    labs(title = "LEC Trajectory Paths - Normal vs Proinflammatory States",
         subtitle = "White curves show inferred developmental trajectories",
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.grid.minor = element_blank())

  # Add trajectory curves
  for(i in seq_along(slingCurves(sds_guided))) {
    curve_data <- data.frame(slingCurves(sds_guided)[[i]]$s[slingCurves(sds_guided)[[i]]$ord, ])
    colnames(curve_data) <- c("UMAP1", "UMAP2")

    p_lec_state <- p_lec_state +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = "black", size = 3, alpha = 0.8) +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = "white", size = 2, alpha = 0.9)
  }

  print(p_lec_state)

  # Save LEC state trajectory plot
  ggsave("LEC_Trajectory_Paths_by_State.pdf", p_lec_state,
         width = 10, height = 8)
}

# =============================================================================
# TRAJECTORY PATH WITH CLUSTER LABELS
# =============================================================================

if(exists("sds_guided")) {

  cat("Creating trajectory paths with cluster labels...\n")

  # Define colors for clusters
  cluster_colors <- brewer.pal(length(unique(sce_LECs$LEC_subcluster)), "Set2")
  names(cluster_colors) <- sort(unique(sce_LECs$LEC_subcluster))

  # Create cluster trajectory plot
  p_clusters <- ggplot(trajectory_base_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Cluster), size = 1.2, alpha = 0.8) +
    scale_color_manual(values = cluster_colors, name = "LEC Cluster") +
    labs(title = "LEC Trajectory Paths - Cluster Identity",
         subtitle = "Trajectories connect LEC developmental states (LEC0-LEC7)",
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.grid.minor = element_blank())

  # Add trajectory curves with different colors for each trajectory
  trajectory_curve_colors <- rainbow(length(slingCurves(sds_guided)))

  for(i in seq_along(slingCurves(sds_guided))) {
    curve_data <- data.frame(slingCurves(sds_guided)[[i]]$s[slingCurves(sds_guided)[[i]]$ord, ])
    colnames(curve_data) <- c("UMAP1", "UMAP2")

    p_clusters <- p_clusters +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = "black", size = 3, alpha = 0.8) +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = trajectory_curve_colors[i], size = 2, alpha = 0.9)
  }

  print(p_clusters)

  # Save cluster trajectory plot
  ggsave("LEC_Trajectory_Paths_by_Cluster.pdf", p_clusters,
         width = 10, height = 8)
}

# =============================================================================
# GENE EXPRESSION ALONG TRAJECTORY CURVES
# =============================================================================

if(exists("sds_guided") && exists("pseudotime_matrix") && length(available_genes) > 0) {

  cat("Creating gene expression plots along trajectory curves...\n")

  # For each trajectory, create detailed gene expression plots
  trajectory_gene_curves <- list()

  for(traj_idx in 1:ncol(pseudotime_matrix)) {

    pt_values <- pseudotime_matrix[, traj_idx]
    valid_cells <- !is.na(pt_values)

    if(sum(valid_cells) > 20) {

      cat(paste("Processing trajectory", traj_idx, "with", sum(valid_cells), "cells...\n"))

      # Get curve data for this trajectory
      curve_data <- data.frame(slingCurves(sds_guided)[[traj_idx]]$s[slingCurves(sds_guided)[[traj_idx]]$ord, ])
      colnames(curve_data) <- c("UMAP1", "UMAP2")

      # Create plots for each gene along this trajectory
      traj_gene_plots <- list()

      for(gene in available_genes) {

        # Get cells assigned to this trajectory
        traj_cells_data <- data.frame(
          UMAP1 = umap_coords[valid_cells, 1],
          UMAP2 = umap_coords[valid_cells, 2],
          Pseudotime = pt_values[valid_cells],
          Gene_Expression = as.numeric(expr_data[gene, valid_cells]),
          Cluster = sce_LECs$LEC_subcluster[valid_cells],
          LEC_state = sce_LECs$LEC_state[valid_cells]
        )

        # Create trajectory-specific gene plot
        p_traj_gene <- ggplot(traj_cells_data, aes(x = UMAP1, y = UMAP2)) +
          geom_point(aes(color = Gene_Expression), size = 1.5, alpha = 0.8) +
          scale_color_viridis_c(name = paste(gene, "\nExpression"), option = "plasma") +
          geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                    color = "black", size = 3, alpha = 0.8) +
          geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                    color = "white", size = 2, alpha = 0.9) +
          labs(title = paste("Trajectory", traj_idx, "-", gene, "Expression"),
               subtitle = paste("n =", sum(valid_cells), "cells assigned to this trajectory"),
               x = "UMAP1", y = "UMAP2") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5),
                panel.grid.minor = element_blank())

        traj_gene_plots[[gene]] <- p_traj_gene
        print(p_traj_gene)
      }

      # Combine genes for this trajectory
      if(length(traj_gene_plots) >= 2) {
        combined_traj_genes <- wrap_plots(traj_gene_plots, ncol = 2)
        trajectory_gene_curves[[traj_idx]] <- combined_traj_genes

        # Save trajectory-specific gene plots
        filename <- paste0("LEC_Trajectory_", traj_idx, "_Gene_Expression_Paths.pdf")
        ggsave(filename, combined_traj_genes, width = 14, height = 10)
        cat("Saved:", filename, "\n")
      }
    }
  }
}

# =============================================================================
# DISEASE GROUP COMPARISON ON TRAJECTORY PATHS
# =============================================================================

if(exists("sds_guided")) {

  cat("Creating disease group comparison on trajectory paths...\n")

  # Define disease colors
  disease_colors <- c("CTRL" = "#2CA02C", "IPF" = "#FF7F00", "SSc-ILD" = "#9467BD")

  # Create disease group trajectory plot
  p_disease <- ggplot(trajectory_base_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = Group), size = 1.2, alpha = 0.8) +
    scale_color_manual(values = disease_colors, name = "Disease Group") +
    labs(title = "LEC Trajectory Paths - Disease Group Distribution",
         subtitle = "How different diseases populate the trajectory paths",
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.title = element_text(size = 12),
          panel.grid.minor = element_blank())

  # Add trajectory curves
  for(i in seq_along(slingCurves(sds_guided))) {
    curve_data <- data.frame(slingCurves(sds_guided)[[i]]$s[slingCurves(sds_guided)[[i]]$ord, ])
    colnames(curve_data) <- c("UMAP1", "UMAP2")

    p_disease <- p_disease +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = "black", size = 3, alpha = 0.8) +
      geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
                color = "white", size = 2, alpha = 0.9)
  }

  print(p_disease)

  # Save disease trajectory plot
  ggsave("LEC_Trajectory_Paths_by_Disease.pdf", p_disease,
         width = 10, height = 8)
}

# =============================================================================
# SUMMARY AND FILE OUTPUTS
# =============================================================================

cat("\n=== TRAJECTORY PATH VISUALIZATION COMPLETE ===\n")
cat("Generated trajectory path plots showing:\n")
cat("1. Gene expression overlay on trajectory curves for each target gene\n")
cat("2. LEC state (Normal vs Proinflammatory) distribution along paths\n")
cat("3. Cluster identity and trajectory connections\n")
cat("4. Disease group distribution on trajectory paths\n")
cat("5. Trajectory-specific gene expression patterns\n\n")

cat("Files saved:\n")
if(length(available_genes) >= 2) {
  cat("- LEC_All_Genes_Trajectory_Paths.pdf (combined gene view)\n")
}
cat("- LEC_Trajectory_Paths_by_State.pdf (Normal vs Proinflammatory)\n")
cat("- LEC_Trajectory_Paths_by_Cluster.pdf (Cluster identity)\n")
cat("- LEC_Trajectory_Paths_by_Disease.pdf (Disease group distribution)\n")
if(exists("trajectory_gene_curves")) {
  cat("- LEC_Trajectory_[N]_Gene_Expression_Paths.pdf (individual trajectories)\n")
}

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("TRAJECTORY PATH PLOTS SHOW:\n\n")

cat("1. WHITE CURVES: Inferred developmental trajectories connecting LEC states\n")
cat("2. COLOR GRADIENTS: Gene expression levels along the trajectory paths\n")
cat("3. TRAJECTORY DIRECTION: From progenitor-like to mature/specialized LECs\n\n")

cat("WHAT TO LOOK FOR:\n")
cat("- PROX1: Should be high in normal LECs, may decrease in inflammatory states\n")
cat("- HIF1A/EPAS1: Expected to increase along proinflammatory trajectories\n")
cat("- SPARCL1: May show different patterns in fibrotic vs normal trajectories\n")
cat("- Path differences: Normal vs proinflammatory LECs may follow different routes\n\n")

cat("BIOLOGICAL SIGNIFICANCE:\n")
cat("- Trajectories represent potential differentiation/activation paths\n")
cat("- Gene expression changes along paths show functional transitions\n")
cat("- Disease-specific trajectory usage indicates altered LEC biology\n\n")
