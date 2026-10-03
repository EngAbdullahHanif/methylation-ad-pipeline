#!/usr/bin/env bash
set -euo pipefail

Rscript scripts/00_download_data.R
Rscript scripts/01_import_and_qc.R
Rscript scripts/02_normalisation.R
Rscript scripts/03_probe_filtering_and_imputation.R
Rscript scripts/04_exploration_and_batch.R
Rscript scripts/05_differential_methylation.R
Rscript scripts/06_enrichment.R
Rscript -e 'rmarkdown::render("report/analysis_report.Rmd")'
