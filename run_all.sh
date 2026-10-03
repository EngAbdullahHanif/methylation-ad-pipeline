#!/usr/bin/env bash
set -euo pipefail

Rscript scripts/00_download_data.R
Rscript scripts/01_import_and_qc.R
Rscript scripts/02_normalisation.R
Rscript scripts/03_exploration_and_differential.R
printf 'Completed milestones M0-M3.\\n'
