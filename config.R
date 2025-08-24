## config.R -- paths and shared helpers for the GSE159354 EC / immune pipeline
##
## Every script in scripts/ starts with `source("config.R")` and is meant to be
## run from the repository root, e.g.
##
##     Rscript scripts/01_load_integrate.R
##
## Two locations are configurable, either by editing the defaults below or by
## setting environment variables before starting R:
##
##   GSE159354_RAW   folder that contains the GEO *_filtered_feature_bc_matrix.h5
##                   files (any depth; files are discovered recursively).
##                   Default: <repo>/data/GSE159354_RAW
##   SC_RESULTS      folder that receives every .rds / .pdf / .csv the pipeline
##                   writes. Default: <repo>/results
##
## The original analysis was written as one R notebook that saved everything
## into its working directory by bare file name. Rather than rewrite several
## hundred file names, config.R switches the working directory to SC_RESULTS
## after loading the helpers, so every relative path in scripts/ lands there.

if (!exists("cfg", inherits = FALSE)) {
  root <- normalizePath(getwd(), mustWork = TRUE)
  if (!file.exists(file.path(root, "config.R"))) {
    stop("Run the scripts from the repository root (the folder that contains config.R).")
  }
  cfg <- list(
    root        = root,
    data_dir    = Sys.getenv("GSE159354_RAW", unset = file.path(root, "data", "GSE159354_RAW")),
    results_dir = Sys.getenv("SC_RESULTS",    unset = file.path(root, "results"))
  )
  dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)

  ## shared helper functions (colors, manifest, heatmaps, DE, GO, CellChat)
  for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) {
    source(f)
  }

  setwd(cfg$results_dir)
  message("GSE159354 pipeline: data in ", cfg$data_dir, "; results in ", cfg$results_dir)
}
