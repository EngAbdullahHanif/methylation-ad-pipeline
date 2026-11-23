#!/usr/bin/env Rscript

# Normalize EPIC data with Noob, compare beta distributions, and filter probes.

options(stringsAsFactors = FALSE)
set.seed(20261003)

if (!requireNamespace("minfi", quietly = TRUE) ||
    !requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  stop("minfi and the EPIC hg19 annotation package are required.")
}

project_dir <- normalizePath(".", mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

rg_set <- readRDS(file.path(results_dir, "rg_set_raw.rds"))
detection_p <- readRDS(file.path(results_dir, "detection_p.rds"))
imported_probe_count <- nrow(rg_set)

# Noob corrects background signal and dye imbalance on the raw red/green channels.
noob_set <- minfi::preprocessNoob(rg_set)
raw_set <- minfi::preprocessRaw(rg_set)
raw_beta <- minfi::getBeta(raw_set)
noob_beta <- minfi::getBeta(noob_set)

pdf(file.path(results_dir, "normalisation_beta_density.pdf"), width = 8, height = 6)
plot(
  density(as.vector(raw_beta), na.rm = TRUE),
  col = "grey40",
  lwd = 2,
  xlab = "Beta value",
  main = "Raw versus Noob beta-value distributions"
)
lines(density(as.vector(noob_beta), na.rm = TRUE), col = "steelblue", lwd = 2)
legend("topright", legend = c("Raw", "Noob"), col = c("grey40", "steelblue"), lwd = 2)
dev.off()
normalisation_summary <- data.frame(
  method = c("raw", "Noob"),
  beta_min = c(min(raw_beta, na.rm = TRUE), min(noob_beta, na.rm = TRUE)),
  beta_median = c(median(raw_beta, na.rm = TRUE), median(noob_beta, na.rm = TRUE)),
  beta_max = c(max(raw_beta, na.rm = TRUE), max(noob_beta, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
rm(raw_set, raw_beta, noob_beta, rg_set)
gc()

annotation_object <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
candidate_probes <- intersect(rownames(noob_set), rownames(detection_p))
annotation <- minfi::getAnnotation(
  annotation_object,
  lociNames = candidate_probes
)
common_probes <- Reduce(intersect, list(
  rownames(noob_set),
  rownames(detection_p),
  rownames(annotation)
))
noob_set <- noob_set[common_probes, ]
detection_p <- detection_p[common_probes, ]
annotation <- annotation[common_probes, ]

# Keep probes detected in at least 80% of the selected samples.
detection_threshold <- 0.01
minimum_detection_fraction <- 0.80
detection_keep <- rowMeans(detection_p <= detection_threshold, na.rm = TRUE) >=
  minimum_detection_fraction

sex_chromosome_keep <- !annotation$chr %in% c("chrX", "chrY")
# These EPIC annotation fields are the CpG and single-base-extension SNP fields
# used by minfi's SNP filter; this avoids creating another full assay object.
snp_keep <- is.na(annotation$CpG_rs) & is.na(annotation$SBE_rs)

mask_version <- "Zhou EPIC v1 hg19 MASK_general (20180909)"
mask_file <- file.path(project_dir, "data", "reference", "EPIC.hg19.MASK_general.20180909.tsv.gz")
masked_probes <- read.delim(mask_file, colClasses = "character")$CpG
if (!length(masked_probes) || anyNA(masked_probes) || anyDuplicated(masked_probes)) {
  stop("Expected a nonempty EPIC v1 mask with unique probe IDs.")
}
before_mask <- detection_keep & sex_chromosome_keep & snp_keep
keep <- before_mask & !common_probes %in% masked_probes
if (!any(keep) || !any(before_mask & !keep)) {
  stop("The EPIC mask did not return a plausible overlap with the QC-passing probes.")
}
filtered_set <- noob_set[keep, ]
filtered_beta <- minfi::getBeta(filtered_set)
missing_before_imputation <- sum(is.na(filtered_beta))

imputation_method <- "none: no missing beta values were observed after filtering"
if (missing_before_imputation > 0L) {
  stop("Missing beta values were observed. Add and document an approved imputation method before continuing.")
}

filter_counts <- data.frame(
  step = c(
    "Imported probes",
    "Common to normalized data, detection p-values, and annotation",
    "Pass detection p-value in at least 80% of samples",
    "Autosomal probes",
    "No CpG or SBE SNP annotation",
    mask_version
  ),
  probes_remaining = c(
    imported_probe_count,
    length(common_probes),
    sum(detection_keep),
    sum(detection_keep & sex_chromosome_keep),
    sum(before_mask),
    sum(keep)
  ),
  probes_removed_at_step = c(
    NA_integer_,
    imported_probe_count - length(common_probes),
    length(common_probes) - sum(detection_keep),
    sum(detection_keep) - sum(detection_keep & sex_chromosome_keep),
    sum(detection_keep & sex_chromosome_keep) - sum(before_mask),
    sum(before_mask & !keep)
  ),
  stringsAsFactors = FALSE
)

utils::write.table(filter_counts, file.path(results_dir, "probe_filter_counts.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(normalisation_summary, file.path(results_dir, "normalisation_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

m2_summary <- data.frame(
  detection_threshold = detection_threshold,
  minimum_detection_fraction = minimum_detection_fraction,
  missing_beta_values = missing_before_imputation,
  imputation_method = imputation_method,
  probe_mask = mask_version,
  mask_excluded_probes = sum(before_mask & !keep),
  stringsAsFactors = FALSE
)
utils::write.table(m2_summary, file.path(results_dir, "m2_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(noob_set, file.path(results_dir, "methyl_set_noob.rds"))
saveRDS(common_probes[before_mask], file.path(results_dir, "probe_universe_before_mask.rds"))
saveRDS(filtered_set, file.path(results_dir, "methyl_set_filtered.rds"))
saveRDS(filtered_beta, file.path(results_dir, "beta_filtered.rds"))

message("Noob normalization and probe filtering complete.")
message("Probes retained: ", nrow(filtered_set))
message("Missing beta values after filtering: ", missing_before_imputation)
