## GEO sample manifest for GSE159354 and the per-sample loader.
##
## The manifest lists all 48 GSM samples in the series (10 donors x sorted
## fractions). This pipeline uses only the CD31+ (endothelial) and CD45+
## (immune) fractions, which exist for every donor: 20 samples in total.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

sample_manifest <- tribble(
  ~gsm,        ~donor,  ~group,    ~fraction,
  "GSM4827158","F19609","IPF",     "EPCAMp_SAbGALp",
  "GSM4827159","F19609","IPF",     "EPCAMp_SAbGALneg",
  "GSM4827160","F19609","IPF",     "CD90p",
  "GSM4827161","F19609","IPF",     "CD31p",
  "GSM4827162","F19609","IPF",     "CD45p",
  "GSM4827163","N1902", "CTRL",    "CD90p",
  "GSM4827164","N1902", "CTRL",    "TripleNeg",
  "GSM4827165","N1902", "CTRL",    "EPCAMp",
  "GSM4827166","N1902", "CTRL",    "CD31p",
  "GSM4827167","N1902", "CTRL",    "CD45p",
  "GSM4827168","S19600","SSc-ILD", "EPCAMp_HT2280neg",
  "GSM4827169","S19600","SSc-ILD", "EPCAMp_HT2280pos",
  "GSM4827170","S19600","SSc-ILD", "TripleNeg",
  "GSM4827171","S19600","SSc-ILD", "CD31p",
  "GSM4827172","S19600","SSc-ILD", "CD45p",
  "GSM4827173","S19584","SSc-ILD", "TripleNeg",
  "GSM4827174","S19584","SSc-ILD", "EPCAMp",
  "GSM4827175","S19584","SSc-ILD", "CD31p",
  "GSM4827176","S19584","SSc-ILD", "CD45p",
  "GSM4827177","S19578","SSc-ILD", "CD90p",
  "GSM4827178","S19578","SSc-ILD", "EPCAMp",
  "GSM4827179","S19578","SSc-ILD", "CD31p",
  "GSM4827180","S19578","SSc-ILD", "CD45p",
  "GSM4827181","N1857", "CTRL",    "TripleNeg",
  "GSM4827182","N1857", "CTRL",    "EPCAMp",
  "GSM4827183","N1857", "CTRL",    "CD31p",
  "GSM4827184","N1857", "CTRL",    "CD45p",
  "GSM4827185","N1857", "CTRL",    "Unsorted",
  "GSM4827186","F18555","IPF",     "TripleNeg",
  "GSM4827187","F18555","IPF",     "EPCAMp",
  "GSM4827188","F18555","IPF",     "CD31p",
  "GSM4827189","F18555","IPF",     "CD45p",
  "GSM4827190","F18555","IPF",     "Unsorted",
  "GSM4827191","F18545","IPF",     "TripleNeg",
  "GSM4827192","F18545","IPF",     "EPCAMp",
  "GSM4827193","F18545","IPF",     "CD31p",
  "GSM4827194","F18545","IPF",     "CD45p",
  "GSM4827195","F18545","IPF",     "Unsorted",
  "GSM4827196","N1865", "CTRL",    "TripleNeg",
  "GSM4827197","N1865", "CTRL",    "EPCAMp",
  "GSM4827198","N1865", "CTRL",    "CD31p",
  "GSM4827199","N1865", "CTRL",    "CD45p",
  "GSM4827200","N1865", "CTRL",    "Unsorted",
  "GSM4827201","F18536","IPF",     "TripleNeg",
  "GSM4827202","F18536","IPF",     "EPCAMp",
  "GSM4827203","F18536","IPF",     "CD31p",
  "GSM4827204","F18536","IPF",     "CD45p",
  "GSM4827205","F18536","IPF",     "Unsorted"
) %>% mutate(group = factor(group, levels = c("CTRL", "IPF", "SSc-ILD")))

#' Discover the filtered_feature_bc_matrix.h5 files under data_dir and join
#' them to the manifest, keeping only the requested fractions.
build_sample_table <- function(data_dir, fractions = c("CD31p", "CD45p"),
                               celltype_mapping = c("CD31p" = "EC", "CD45p" = "Immune")) {
  h5_file_paths <- list.files(data_dir, pattern = "filtered_feature_bc_matrix\\.h5$",
                              full.names = TRUE, recursive = TRUE)
  if (length(h5_file_paths) == 0) {
    stop("No *filtered_feature_bc_matrix.h5 files found under ", data_dir,
         ". Download GSE159354_RAW from GEO and point GSE159354_RAW at it (see config.R).")
  }
  file_discovery_tbl <- tibble(
    file  = normalizePath(h5_file_paths),
    fname = basename(h5_file_paths),
    gsm   = sub("^(GSM\\d+).*", "\\1", basename(h5_file_paths))
  )
  sample_manifest %>%
    left_join(file_discovery_tbl, by = "gsm") %>%
    filter(fraction %in% fractions, !is.na(file), file.exists(file)) %>%
    mutate(cellset   = recode(fraction, !!!celltype_mapping),
           group_chr = as.character(group),
           sample_id = paste(group_chr, cellset, donor, gsm, sep = "_"))
}

#' Read one 10x h5 file into a QC-filtered Seurat object with sample metadata.
#' QC: genes detected in >= 3 cells, 200 < nFeature_RNA < 8000, percent.mt < 20.
load_sample <- function(row, min_features = 200, max_features = 8000, max_mt = 20) {
  count_matrix <- Read10X_h5(row$file)
  rownames(count_matrix) <- make.unique(rownames(count_matrix))
  seurat_obj <- CreateSeuratObject(counts = count_matrix, project = row$sample_id,
                                   min.cells = 3, min.features = min_features)
  seurat_obj$gsm      <- row$gsm
  seurat_obj$donor    <- row$donor
  seurat_obj$group    <- factor(row$group, levels = c("CTRL", "IPF", "SSc-ILD"))
  seurat_obj$fraction <- row$fraction
  seurat_obj$cellset  <- row$cellset
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features &
                                            nFeature_RNA < max_features &
                                            percent.mt    < max_mt)
  ## prefix cell barcodes with sample_id so merged barcodes stay unique and readable
  seurat_obj <- RenameCells(seurat_obj, new.names = paste(row$sample_id, Cells(seurat_obj), sep = ":"))
  seurat_obj
}
