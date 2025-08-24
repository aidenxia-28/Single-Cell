# data/

Put the GEO supplementary archive here:

```
data/GSE159354_RAW/
    GSM4827161_..._filtered_feature_bc_matrix.h5
    GSM4827162_..._filtered_feature_bc_matrix.h5
    ...
```

Download `GSE159354_RAW.tar` from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159354,
extract it into this folder (subfolders are fine; files are discovered recursively), and
`scripts/01_load_integrate.R` will pick up the 20 CD31+ / CD45+ samples listed in
`R/manifest.R`. Any other location works by setting `GSE159354_RAW=/path/to/GSE159354_RAW`.

The raw files are git-ignored.
