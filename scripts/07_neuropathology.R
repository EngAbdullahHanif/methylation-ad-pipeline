#!/usr/bin/env Rscript

source("scripts/extension_models.R")
cohort <- read_extension_assays("pathology")
metadata <- cohort$metadata
if (anyNA(metadata[, c("batch", "sex", "proportion_neurons")])) stop("Incomplete neuropathology covariates.")
included <- !is.na(metadata$niaa_b)
metadata <- metadata[included, ]

# The prespecified measure is the NIA-AA B score for neurofibrillary tangles.
if (!setequal(unique(metadata$niaa_b), c(2L, 3L))) {
  stop("Expected NIA-AA B scores of 2 and 3 for the prespecified one-level contrast.")
}
formula <- ~ niaa_b + batch + sex + proportion_neurons
design <- model.matrix(formula, metadata)
result <- fit_extension(cohort$m[, included, drop = FALSE], cohort$beta[, included, drop = FALSE], design, "niaa_b")

qc <- read.delim("results/extensions/pathology/qc_summary.tsv")
score_by_batch <- as.data.frame.matrix(with(metadata, table(batch, niaa_b)))
write_extension_table(data.frame(batch = rownames(score_by_batch), n_b2 = score_by_batch[["2"]],
                                 n_b3 = score_by_batch[["3"]]), "neuropathology_score_by_batch")
summary <- cbind(data.frame(region = "middle frontal gyrus", contrast = "NIA-AA B3 minus B2",
                            model = deparse1(formula), n_eligible_arrays = nrow(qc),
                            n_qc_pass = sum(qc$qc_status == "pass"), n_samples = nrow(metadata),
                            n_b2 = sum(metadata$niaa_b == 2), n_b3 = sum(metadata$niaa_b == 3),
                            n_batches_with_both_scores = sum(score_by_batch[["2"]] > 0 & score_by_batch[["3"]] > 0),
                            n_intermediate_diagnoses = sum(metadata$conf_diag == 1)),
                 extension_statistics(result, design))
write_extension_table(summary, "neuropathology_summary")
write_extension_table(head(result, 100), "neuropathology_top")
saveRDS(result, "results/extensions/pathology/site_statistics.rds")
plot_extension_effects(result, "neuropathology.pdf", "Adjusted delta beta (NIA-AA B3 minus B2), 95% CI")
