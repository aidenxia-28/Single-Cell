# What changed between the working notebook and this repository

The analysis was developed interactively in one R notebook
(`notebooks/original_analysis.Rmd`, ~13,500 lines, 64 chunks). This repository is the
same analysis, re-organised so it can be run and read. The computations were not changed;
the list below is complete.

## Structure

- One numbered script per stage in `scripts/`; the stage order and the file names each
  stage reads and writes are in `docs/pipeline.md`. Scripts are meant to be run from the
  repository root (`Rscript scripts/01_load_integrate.R`) or all together with
  `scripts/run_all.sh`.
- `config.R` holds the two configurable paths (raw data, results) and sources `R/`.
  The notebook saved everything into its working directory by bare file name; instead of
  rewriting several hundred file names, `config.R` switches the working directory to the
  results folder, so all those names are unchanged.
- Code that the notebook repeated with only an object or cluster name changed is now one
  function in `R/`:
  - `R/marker_panel_heatmap.R` – the EC/immune marker-panel heatmap (was five copies).
  - `R/pseudobulk_de.R` – donor-level pseudobulk DESeq2 between two LEC subclusters with
    volcano plots (was three copies: LEC1 vs LEC3, LEC1 vs LEC6, LEC1 vs LEC5).
  - `R/go_enrichment.R` – GO enrichment of those gene lists (was two copies).
  - `R/cellchat_lec_immune.R` – LEC/immune CellChat (was two copies plus three heatmap
    variants).
  - `R/manifest.R` – the GEO sample manifest and the per-sample loader.
  - `R/colors.R` – the palettes.
- The LEC signature dot plot (two near-identical chunks) is one function inside
  `scripts/22_lec_signature_dotplots.R`; the per-gene violin/feature plots are a loop in
  `scripts/21_lec_marker_plots.R`.
- Branches that were tried and not used downstream are kept under `scripts/exploratory/`
  (ultra-clean and cluster-based EC filtering; immune re-clustering straight from the
  integrated object).
- Chunks dropped as exact duplicates or print-only repeats: the second copy of the
  "remove clusters 13/20/22/24" block, the "no labels" re-plots of the EC-subtype UMAPs,
  the first (superseded) LEC1-vs-LEC3 volcano block, the broken first CellChat attempt that
  referenced an undefined `cell_labels` vector, and the empty chunks.
- `docs/exploratory/research.nb.html` is the rendered output of an earlier, different
  version of this analysis ("Single Cell 3": adaptive QC, Harmony on donor + fraction +
  group, EC subtypes with a 0.3-SD margin, disease-only LEC subclustering). It is kept
  because it is the only rendered output in the repository; its numbers do not correspond
  to the scripts here.

## Mechanical edits

- Machine-specific `setwd()` / `getwd()` and the hard-coded `/Users/.../Single Cell 4`
  data path were replaced by `config.R`.
- All `install.packages()`, `BiocManager::install()` and `devtools::install_github()`
  calls were removed; package requirements are listed in the README. The Monocle3
  install-on-the-fly block became a `requireNamespace()` check (Slingshot always runs).
- `system("open ...")` calls that opened PDFs in Preview were removed.

## Bugs fixed

- Stage 01: a `readRDS()` of the integrated object placed before the object was saved
  (a leftover from re-running the chunk) and a stray `1` at the start of a comment line.
- Stage 21: `VlnPlot(features = "")` for the IL3RA panel now reads `"IL3RA"`.
- Stage 24: `genes_to_plot` listed `PDPN` twice.
- Stage 31: `cat("\n" + paste(...))` – R has no `+` for strings; now `paste0`.
- Stage 33: `saveRDS(renamed_markers, ...)` ran before `renamed_markers` existed; moved
  after `FindAllMarkers()`. The M1/M2 feature-plot chunk had an empty gene vector
  (`c("", "", "IL-10", "CD209")`) and referred to an object from another chunk; it now
  loops over IL10 / CD209 / TLR4 / CD4 on the renumbered object.
- Stage 40: the summary printed `ncol(combined_matrix)`, which did not exist
  (`combined_raw`).
