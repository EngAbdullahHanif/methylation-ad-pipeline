#!/usr/bin/env Rscript

# Build the analysis sample sheet, import EPIC IDATs, and perform basic QC.

options(stringsAsFactors = FALSE)
set.seed(20261003)

if (!requireNamespace("GEOquery", quietly = TRUE) ||
    !requireNamespace("minfi", quietly = TRUE)) {
  stop("GEOquery and minfi must be installed before running this script.")
}

project_dir <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "raw", "GSE212682", "idats")
metadata_dir <- file.path(project_dir, "data", "metadata")
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

accession <- "GSE212682"
region_to_keep <- "middle frontal gyrus"
max_per_group <- 5L

gse <- GEOquery::getGEO(accession, GSEMatrix = TRUE, destdir = metadata_dir)[[1]]
geo_metadata <- Biobase::pData(gse)
characteristic_columns <- grep("^characteristics_ch1", names(geo_metadata), value = TRUE)
characteristic_table <- as.data.frame(geo_metadata[, characteristic_columns, drop = FALSE])

find_characteristic <- function(prefix) {
  for (column in names(characteristic_table)) {
    values <- characteristic_table[[column]]
    matching <- grepl(paste0("^", prefix, ":"), values, ignore.case = TRUE)
    if (any(matching)) {
      result <- rep(NA_character_, nrow(geo_metadata))
      result[matching] <- sub(paste0("^", prefix, ":\\s*"), "", values[matching],
                              ignore.case = TRUE)
      return(result)
    }
  }
  rep(NA_character_, nrow(geo_metadata))
}

sample_sheet <- data.frame(
  sample_id = rownames(geo_metadata),
  title = geo_metadata$title,
  region = find_characteristic("brain_region"),
  region_abbreviation = find_characteristic("brain_region_abbreviation"),
  sex = find_characteristic("Sex"),
  age_raw = find_characteristic("age"),
  conf_diag = find_characteristic("conf_diag"),
  adseverityscore = find_characteristic("adseverityscore"),
  bisulfite_batch = find_characteristic("bisulfite_batch"),
  batch = find_characteristic("batch"),
  stringsAsFactors = FALSE
)

sample_sheet$diagnosis <- ifelse(
  sample_sheet$conf_diag == "0", "Control",
  ifelse(sample_sheet$conf_diag == "2", "AD", NA_character_)
)
sample_sheet$age_numeric <- suppressWarnings(as.numeric(sample_sheet$age_raw))

idat_files <- list.files(raw_dir, pattern = "\\.idat\\.gz$", full.names = TRUE,
                         ignore.case = TRUE)
idat_names <- basename(idat_files)
idat_map <- data.frame(
  sample_id = sub("_.*$", "", idat_names),
  idat_file = idat_files,
  stringsAsFactors = FALSE
)
green_files <- idat_map[grepl("_Grn\\.idat\\.gz$", idat_map$idat_file, ignore.case = TRUE), ]
green_files$basename <- sub("_Grn\\.idat\\.gz$", "", basename(green_files$idat_file),
                            ignore.case = TRUE)
green_files$sentrix_id <- sub("^[^_]+_([^_]+)_.*$", "\\1", green_files$basename)
green_files$sentrix_position <- sub("^.*_(R[0-9]+C[0-9]+)$", "\\1", green_files$basename)
green_files$Basename <- file.path(dirname(green_files$idat_file), green_files$basename)

sample_sheet <- merge(sample_sheet, green_files[, c(
  "sample_id", "sentrix_id", "sentrix_position", "Basename"
)], by = "sample_id", all.x = TRUE, sort = FALSE)
sample_sheet <- sample_sheet[
  sample_sheet$region == region_to_keep &
    sample_sheet$diagnosis %in% c("Control", "AD") &
    !is.na(sample_sheet$Basename),
]

set.seed(20261003)
sample_sheet <- do.call(rbind, lapply(split(sample_sheet, sample_sheet$diagnosis), function(group) {
  group[sample(seq_len(nrow(group)), min(nrow(group), max_per_group)), , drop = FALSE]
}))
sample_sheet <- sample_sheet[order(sample_sheet$diagnosis, sample_sheet$sample_id), ]
rownames(sample_sheet) <- NULL

if (nrow(sample_sheet) < 4L || length(unique(sample_sheet$diagnosis)) != 2L) {
  stop("The selected region does not contain both diagnosis groups after filtering.")
}

utils::write.table(
  sample_sheet,
  file.path(metadata_dir, "sample_sheet_clean.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Selected samples: ", nrow(sample_sheet))
print(table(sample_sheet$diagnosis))

rg_set <- minfi::read.metharray.exp(
  targets = sample_sheet,
  extended = TRUE,
  force = TRUE
)
minfi::pData(rg_set)$diagnosis <- sample_sheet$diagnosis
minfi::pData(rg_set)$sample_id <- sample_sheet$sample_id

normalized_set <- minfi::preprocessRaw(rg_set)
genomic_set <- minfi::mapToGenome(normalized_set)
sex_prediction <- minfi::getSex(genomic_set)
rm(normalized_set, genomic_set)
gc()

detection_p <- minfi::detectionP(rg_set)
mean_detection_p <- colMeans(detection_p, na.rm = TRUE)
qc_threshold <- 0.01
qc_status <- ifelse(mean_detection_p <= qc_threshold, "pass", "fail")
qc_summary <- data.frame(
  sample_id = sample_sheet$sample_id,
  diagnosis = sample_sheet$diagnosis,
  mean_detection_p = mean_detection_p,
  qc_threshold = qc_threshold,
  qc_status = qc_status,
  stringsAsFactors = FALSE
)

normalized_set <- minfi::preprocessRaw(rg_set)
genomic_set <- minfi::mapToGenome(normalized_set)
sex_check <- data.frame(
  sample_id = sample_sheet$sample_id,
  metadata_sex = sample_sheet$sex,
  predicted_sex = as.character(sex_prediction$predictedSex),
  stringsAsFactors = FALSE
)
sex_check$sex_mismatch <- !is.na(sex_check$metadata_sex) &
  sex_check$metadata_sex != "" &
  tolower(substr(sex_check$metadata_sex, 1L, 1L)) !=
    tolower(substr(sex_check$predicted_sex, 1L, 1L))

pdf(file.path(results_dir, "qc_mean_detection_p.pdf"), width = 9, height = 5)
barplot(
  qc_summary$mean_detection_p,
  names.arg = qc_summary$sample_id,
  las = 2,
  cex.names = 0.45,
  ylab = "Mean detection p-value",
  main = "Sample detection p-values"
)
abline(h = qc_threshold, col = "red", lty = 2)
dev.off()

pdf(file.path(results_dir, "qc_intensity.pdf"), width = 7, height = 7)
minfi::plotQC(rg_set)
dev.off()

saveRDS(rg_set, file.path(results_dir, "rg_set_raw.rds"))
saveRDS(detection_p, file.path(results_dir, "detection_p.rds"))
saveRDS(qc_summary, file.path(results_dir, "qc_summary.rds"))
saveRDS(sex_check, file.path(results_dir, "sex_check.rds"))
utils::write.table(qc_summary, file.path(results_dir, "qc_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(sex_check, file.path(results_dir, "sex_check.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

message("QC complete. Samples failing the mean detection p-value threshold: ",
        sum(qc_summary$qc_status == "fail"))
message("Sex mismatches flagged: ", sum(sex_check$sex_mismatch, na.rm = TRUE))
