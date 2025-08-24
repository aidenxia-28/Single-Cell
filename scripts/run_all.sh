#!/usr/bin/env bash
# Run the whole GSE159354 EC / immune pipeline in order.
#
#   bash scripts/run_all.sh                 # everything
#   bash scripts/run_all.sh ec              # stages 01-04 + endothelial / lymphatic branch
#   bash scripts/run_all.sh immune          # stages 01-04 + immune branch
#   GSE159354_RAW=/path/to/GSE159354_RAW SC_RESULTS=/path/to/results bash scripts/run_all.sh
#
# Each stage reads the .rds of the previous one from $SC_RESULTS (default ./results),
# so a stage can be re-run on its own with Rscript once its inputs exist.
set -euo pipefail
cd "$(dirname "$0")/.."

core=(01_load_integrate 02_remove_clusters 03_doublet_contamination_qc 04_marker_panel_heatmaps)
ec=(10_ec_subtypes 11_ec_cluster_filter 12_ec_subtype_heatmap
    20_lec_subcluster 21_lec_marker_plots 22_lec_signature_dotplots 23_lec_heatmap 24_lec_trajectory
    25_lec_within_cluster_disease_de 26_lec_cluster_vs_cluster_de 27_lec_go_enrichment)
immune=(30_immune_subcluster 31_immune_markers 32_immune_filter_renumber 33_immune_renamed_markers_heatmaps)
joint=(40_cellchat_lec_immune)

case "${1:-all}" in
  ec)     stages=("${core[@]}" "${ec[@]}") ;;
  immune) stages=("${core[@]}" "${immune[@]}") ;;
  all)    stages=("${core[@]}" "${ec[@]}" "${immune[@]}" "${joint[@]}") ;;
  *)      echo "usage: $0 [all|ec|immune]" >&2; exit 1 ;;
esac

for s in "${stages[@]}"; do
  echo "================ $s ================"
  Rscript "scripts/$s.R"
done
echo "Pipeline finished. Results in ${SC_RESULTS:-$(pwd)/results}"
