## Immune branch: drop contaminating clusters and renumber the rest 0-20
## Stage 32 of the GSE159354 EC / immune pipeline.
##
## Inputs:  Immune_only_seurat_processed.rds (stage 30/31)
## Outputs: Immune_Renamed_0to20_Final.rds, Final_*.pdf, Clean_Immune_Cluster_Disease_Composition.pdf/.png
## Resolution-0.8 clusters 13 (epithelial contamination), 20 (unknown), 22
## (proliferating B cells) and 24 (mast / basophils) are removed; the remaining
## 21 clusters are renumbered 0-20 and stored in metadata column
## renamed_clusters.
##
## Run from the repository root:  Rscript scripts/<this file>

source("config.R")   # paths, shared helpers; switches to the results directory

# ===== USE SEURAT'S BUILT-IN RENAMING FUNCTIONS =====
library(Seurat)
library(dplyr)
library(ggplot2)

# Load original immune object
if (file.exists("Immune_only_seurat_processed.rds")) {
  immune_seurat_original <- readRDS("Immune_only_seurat_processed.rds")
} else {
  clean_lung_seurat <- readRDS("GSE159354_Clean_EC-Immune_harmony.rds")
  immune_seurat_original <- subset(clean_lung_seurat, subset = cellset == "Immune")
}

cat("Original clusters:", paste(sort(unique(immune_seurat_original$RNA_snn_res.0.8)), collapse = ", "), "\n")

# Remove specified clusters
clusters_to_remove <- c("13", "20", "22", "24")
immune_filtered <- subset(immune_seurat_original,
                         subset = RNA_snn_res.0.8 %in% setdiff(levels(immune_seurat_original$RNA_snn_res.0.8),
                                                              clusters_to_remove))

cat("After filtering clusters:", paste(sort(unique(immune_filtered$RNA_snn_res.0.8)), collapse = ", "), "\n")

# ===== USE RenameIdents() TO RENAME CLUSTERS =====
# Set the current clustering as active identity
Idents(immune_filtered) <- immune_filtered$RNA_snn_res.0.8

# Get remaining clusters in order
remaining_clusters <- sort(as.numeric(levels(immune_filtered$RNA_snn_res.0.8)))
cat("Remaining clusters:", paste(remaining_clusters, collapse = ", "), "\n")

# Create renaming vector: old_name = new_name
new_names <- 0:(length(remaining_clusters) - 1)
rename_vector <- setNames(as.character(new_names), as.character(remaining_clusters))

cat("Renaming mapping:\n")
for (i in 1:length(rename_vector)) {
  cat("  Cluster", names(rename_vector)[i], "-> Cluster", rename_vector[i], "\n")
}

# Apply the renaming using RenameIdents
immune_filtered <- RenameIdents(immune_filtered, rename_vector)

# Add the new cluster names to metadata
immune_filtered$renamed_clusters <- Idents(immune_filtered)

cat("Renaming complete! New clusters:", paste(sort(unique(Idents(immune_filtered))), collapse = ", "), "\n")

# ===== CREATE DISEASE STATUS =====
immune_filtered$disease_status <- ifelse(immune_filtered$group == "CTRL", "CTRL", "Disease")
immune_filtered$disease_status <- factor(immune_filtered$disease_status, levels = c("CTRL", "Disease"))

# ===== CREATE PLOTS =====
clusters_plot <- DimPlot(immune_filtered,
                        reduction = "umap",
                        label = TRUE,
                        label.size = 4) +
                 ggtitle("Renumbered Immune Clusters (0-20)") +
                 theme_minimal() +
                 theme(legend.position = "none")

disease_groups_plot <- DimPlot(immune_filtered,
                              group.by = "group",
                              reduction = "umap",
                              cols = c("CTRL" = "#2E8B57", "IPF" = "#9370DB", "SSc-ILD" = "#DC143C")) +
                      ggtitle("Immune Cells: Disease Groups") +
                      theme_minimal()

ctrl_disease_plot <- DimPlot(immune_filtered,
                            group.by = "disease_status",
                            reduction = "umap",
                            cols = c("CTRL" = "#2E8B57", "Disease" = "#DC143C")) +
                    ggtitle("Immune Cells: CTRL vs Disease") +
                    theme_minimal()

print(clusters_plot)
print(disease_groups_plot)
print(ctrl_disease_plot)

# ===== SUMMARY TABLE =====
cat("\n=== CLUSTER RENAMING SUMMARY ===\n")
cluster_counts <- table(Idents(immune_filtered))
print(cluster_counts)

disease_summary <- immune_filtered@meta.data %>%
  group_by(group) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  mutate(percentage = round((n_cells / sum(n_cells)) * 100, 1))

cat("\nDisease group distribution:\n")
print(disease_summary)

# ===== SAVE RESULTS =====
saveRDS(immune_filtered, file = "Immune_Renamed_0to20_Final.rds")

ggsave("Final_Renumbered_Clusters.pdf", clusters_plot, width = 10, height = 8)
ggsave("Final_Disease_Groups.pdf", disease_groups_plot, width = 10, height = 8)
ggsave("Final_CTRL_vs_Disease.pdf", ctrl_disease_plot, width = 10, height = 8)

# Create combined PDF
pdf("Final_Immune_Analysis_Renamed.pdf", height = 12, width = 16)
print(clusters_plot)
print(disease_groups_plot)
print(ctrl_disease_plot)

library(gridExtra)
grid.arrange(disease_groups_plot, ctrl_disease_plot, ncol = 2,
            top = "Immune Cells: Disease Analysis (Renumbered Clusters)")
dev.off()

# ===== FINAL SUCCESS MESSAGE =====
cat("\n=== SUCCESS! CLUSTERS RENUMBERED ===\n")
cat("Clusters successfully renamed from 0-", length(remaining_clusters)-1, "\n")
cat("Total cells:", ncol(immune_filtered), "\n")
cat("Active identity now contains renumbered clusters\n")
cat("Files saved with 'Final_' prefix\n")

# Update your main object
immune_seurat <- immune_filtered
cat("\n'immune_seurat' now has properly renumbered clusters!\n")


clusters_plot2 <- DimPlot(immune_seurat,
                        reduction = "umap",
                        label = FALSE) +  # Remove labels on plot
                 ggtitle("Immune Clusters") +
                 theme_minimal() +
                 theme(legend.position = "right")  # Show legend

print(clusters_plot2)


# ===== CLEAN IMMUNE CLUSTER DISEASE COMPOSITION PLOT =====
library(tidyr)

# ===== 1) Use existing composition data =====
# (assuming composition_data already exists from previous code)

# If composition_data doesn't exist, recreate it
if (!exists("composition_data")) {
  composition_data <- immune_seurat@meta.data %>%
    mutate(cluster = as.character(renamed_clusters)) %>%
    group_by(cluster, group) %>%
    summarise(n_cells = n(), .groups = 'drop') %>%
    group_by(cluster) %>%
    mutate(
      total_in_cluster = sum(n_cells),
      percentage = (n_cells / total_in_cluster) * 100
    ) %>%
    ungroup() %>%
    complete(cluster, group, fill = list(n_cells = 0, percentage = 0)) %>%
    group_by(cluster) %>%
    mutate(total_in_cluster = sum(n_cells[n_cells > 0])) %>%
    mutate(percentage = ifelse(total_in_cluster > 0, (n_cells / total_in_cluster) * 100, 0)) %>%
    ungroup()

  composition_data$cluster <- factor(composition_data$cluster,
                                    levels = as.character(sort(as.numeric(unique(composition_data$cluster)))))
}

# ===== 2) Create clean stacked bar plot =====
# Define colors
disease_colors <- c("CTRL" = "#2E8B57",      # Green
                   "IPF" = "#FF8C00",        # Orange
                   "SSc-ILD" = "#8A2BE2")    # Purple

clean_composition_plot <- ggplot(composition_data, aes(x = cluster, y = percentage, fill = group)) +
  geom_bar(stat = "identity",
           color = "white",
           size = 0.3,
           width = 0.8) +  # Make bars wider
  scale_fill_manual(values = disease_colors,
                   name = "Disease Group") +
  scale_y_continuous(limits = c(0, 100),
                    breaks = seq(0, 100, 25),
                    expand = c(0, 0)) +
  labs(
    x = "Immune Cell Cluster",
    y = "Percentage",
    title = "Disease Group Composition by Immune Cell Cluster"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black", margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.ticks.length = unit(0.2, "cm")
  )

print(clean_composition_plot)

# ===== 3) Save the clean plot =====
ggsave("Clean_Immune_Cluster_Disease_Composition.pdf", clean_composition_plot,
       width = 14, height = 7, dpi = 300)

ggsave("Clean_Immune_Cluster_Disease_Composition.png", clean_composition_plot,
       width = 14, height = 7, dpi = 300)

# ===== 4) Alternative version with even wider bars =====
extra_wide_plot <- ggplot(composition_data, aes(x = cluster, y = percentage, fill = group)) +
  geom_bar(stat = "identity",
           color = "white",
           size = 0.2,
           width = 0.9) +  # Even wider bars
  scale_fill_manual(values = disease_colors,
                   name = "Disease Group") +
  scale_y_continuous(limits = c(0, 100),
                    breaks = seq(0, 100, 25),
                    expand = c(0, 0)) +
  labs(
    x = "Immune Cell Cluster",
    y = "Percentage",
    title = "Disease Group Composition by Immune Cell Cluster"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.4)
  )

ggsave("Extra_Wide_Immune_Composition.pdf", extra_wide_plot,
       width = 16, height = 8, dpi = 300)

# ===== 5) Create summary without cluttering numbers =====
cat("\n=== CLEAN COMPOSITION PLOT CREATED ===\n")
cat("Improvements made:\n")
cat("- Wider bars (width = 0.8 and 0.9)\n")
cat("- Removed all percentage labels from bars\n")
cat("- Cleaner appearance with better spacing\n")
cat("- Larger plot dimensions for better visibility\n")

cat("\nFiles generated:\n")
cat("1. Clean_Immune_Cluster_Disease_Composition.pdf - Clean version (14x7)\n")
cat("2. Clean_Immune_Cluster_Disease_Composition.png - PNG version\n")
cat("3. Extra_Wide_Immune_Composition.pdf - Extra wide version (16x8)\n")

cat("\nColors:\n")
cat("- Green (#2E8B57): CTRL\n")
cat("- Orange (#FF8C00): IPF\n")
cat("- Purple (#8A2BE2): SSc-ILD\n")

# ===== 6) Print some key insights without cluttering the plot =====
# Most disease-enriched clusters

disease_enriched <- composition_data %>%
  filter(group != "CTRL") %>%
  group_by(cluster) %>%
  summarise(disease_pct = sum(percentage), .groups = 'drop') %>%
  arrange(desc(disease_pct)) %>%
  head(5)

cat("\nTop 5 disease-enriched clusters (IPF + SSc-ILD combined):\n")
for(i in 1:nrow(disease_enriched)) {
  cat("  Cluster", disease_enriched$cluster[i], ":", round(disease_enriched$disease_pct[i], 1), "% disease\n")
}
