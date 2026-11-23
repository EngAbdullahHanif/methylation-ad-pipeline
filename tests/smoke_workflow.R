# Run from the repository root with Rscript tests/smoke_workflow.R.
metadata <- read.delim("data/metadata/sample_sheet_clean.tsv")
qc <- read.delim("results/qc_summary.tsv")
pairs <- read.delim("results/batch_overlap.tsv")
design <- read.delim("results/model_design.tsv")
pathology <- read.delim("results/neuropathology_summary.tsv")
pathology_batches <- read.delim("results/neuropathology_score_by_batch.tsv")
regional <- read.delim("results/regional_specificity_summary.tsv")
donor_pairs <- read.delim("results/regional_specificity_pairs.tsv")
stopifnot(pathology$n_samples == sum(pathology_batches$n_b2 + pathology_batches$n_b3),
          pathology$n_b2 == sum(pathology_batches$n_b2), pathology$n_b3 == sum(pathology_batches$n_b3),
          pathology$n_batches_with_both_scores == sum(pathology_batches$n_b2 > 0 & pathology_batches$n_b3 > 0),
          pathology$n_samples <= pathology$n_qc_pass, pathology$n_qc_pass <= pathology$n_eligible_arrays,
          pathology$n_intermediate_diagnoses > 0L, pathology$residual_df == pathology$n_samples - pathology$design_rank,
          !anyDuplicated(donor_pairs$donor_id), !anyDuplicated(donor_pairs$sample_id_mfg),
          !anyDuplicated(donor_pairs$sample_id_dg), !any(donor_pairs$sample_id_dg %in% donor_pairs$sample_id_mfg),
          all(donor_pairs$diagnosis_mfg == donor_pairs$diagnosis_dg), all(donor_pairs$sex_mfg == donor_pairs$sex_dg),
          regional$n_pairs == nrow(donor_pairs), regional$n_control_pairs == sum(donor_pairs$diagnosis_mfg == "Control"),
          regional$n_ad_pairs == sum(donor_pairs$diagnosis_mfg == "AD"),
          regional$residual_df == regional$n_pairs - regional$design_rank,
          isTRUE(all.equal(sum(donor_pairs$diagnosis_information_share), 1)),
          isTRUE(all.equal(regional$max_donor_diagnosis_information_share, max(donor_pairs$diagnosis_information_share))))
extension_summaries <- list(neuropathology = pathology, regional_specificity = regional)
for (name in names(extension_summaries)) {
  info <- extension_summaries[[name]]
  candidates <- read.delim(file.path("results", paste0(name, "_top.tsv")))
  stopifnot(info$design_rank == info$design_columns, info$residual_df > 0,
            nrow(candidates) == min(100L, info$n_probes), !anyDuplicated(candidates$CpG),
            !is.unsorted(candidates$p_value), all(candidates$ci_95_low <= candidates$delta_beta),
            all(candidates$delta_beta <= candidates$ci_95_high),
            file.info(file.path("results", paste0(name, ".pdf")))$size > 0)
}

# Refit the leading beta-scale estimates independently when local assays are available.
pathology_dir <- "results/extensions/pathology"
dentate_dir <- "results/extensions/dentate"
if (file.exists(file.path(pathology_dir, "beta_filtered.rds")) &&
    file.exists(file.path(dentate_dir, "beta_filtered.rds"))) {
  pathology_metadata <- read.delim(file.path(pathology_dir, "sample_sheet_clean.tsv"))
  pathology_beta <- readRDS(file.path(pathology_dir, "beta_filtered.rds"))
  dentate_beta <- readRDS(file.path(dentate_dir, "beta_filtered.rds"))
  stopifnot(regional$n_probes == length(intersect(rownames(pathology_beta), rownames(dentate_beta))))
  sample_ids <- function(beta) sub("_.*$", "", colnames(beta))
  pathology_metadata <- subset(pathology_metadata, !is.na(niaa_b))
  top_pathology <- read.delim("results/neuropathology_top.tsv")[1L, ]
  pathology_metadata$beta <- as.numeric(pathology_beta[top_pathology$CpG,
                                                       match(pathology_metadata$sample_id, sample_ids(pathology_beta))])
  expected <- coef(lm(beta ~ niaa_b + batch + sex + proportion_neurons, pathology_metadata))["niaa_b"]
  stopifnot(isTRUE(all.equal(unname(expected), top_pathology$delta_beta)))

  top_regional <- read.delim("results/regional_specificity_top.tsv")[1L, ]
  batches <- sort(unique(c(donor_pairs$batch_dg, donor_pairs$batch_mfg)))
  batch_differences <- outer(donor_pairs$batch_dg, batches, "==") - outer(donor_pairs$batch_mfg, batches, "==")
  covariates <- model.matrix(~ sex_mfg + I(proportion_neurons_dg - proportion_neurons_mfg), donor_pairs)
  paired_design <- cbind(covariates, batch_differences, diagnosisAD = donor_pairs$diagnosis_mfg == "AD")
  observed <- dentate_beta[top_regional$CpG, match(donor_pairs$sample_id_dg, sample_ids(dentate_beta))] -
    pathology_beta[top_regional$CpG, match(donor_pairs$sample_id_mfg, sample_ids(pathology_beta))]
  expected <- lm.fit(paired_design, observed)$coefficients["diagnosisAD"]
  stopifnot(qr(paired_design)$rank == regional$design_rank,
            isTRUE(all.equal(unname(expected), top_regional$delta_beta)))
  rm(pathology_beta, dentate_beta)
}

stopifnot(nrow(qc) == 29L, sum(qc$diagnosis == "Control") == 18L, sum(qc$diagnosis == "AD") == 11L,
          setequal(metadata$sample_id, qc$sample_id[qc$qc_status == "pass"]),
          all(metadata$age_raw == "90+"), all(is.na(metadata$age_numeric)),
          all(abs(rowSums(metadata[startsWith(names(metadata), "proportion_")]) - 1) < 1e-4))

batch_counts <- with(metadata, table(batch, diagnosis))
stopifnot(identical(pairs$batch, rownames(batch_counts)),
          identical(pairs$n_control, as.integer(batch_counts[, "Control"])),
          identical(pairs$n_ad, as.integer(batch_counts[, "AD"])),
          identical(pairs$shared, unname(rowSums(batch_counts > 0) == 2L)),
          identical(design$analysis, c("primary", "shared_batches")))
cohorts <- list(metadata, metadata[metadata$batch %in% pairs$batch[pairs$shared], ])
design_matrices <- lapply(seq_along(cohorts), function(i) {
  cohort <- cohorts[[i]]
  cohort$diagnosis <- factor(cohort$diagnosis, levels = c("Control", "AD"))
  formula <- as.formula(design$model[i])
  expected_terms <- c("diagnosis", "sex", "proportion_neurons")
  if (length(unique(cohort$batch)) > 1L) expected_terms <- c(expected_terms, "batch")
  matrix <- model.matrix(formula, data = cohort)
  stopifnot(setequal(attr(terms(formula), "term.labels"), expected_terms),
            design$n_samples[i] == nrow(cohort), design$n_control[i] == sum(cohort$diagnosis == "Control"),
            design$n_ad[i] == sum(cohort$diagnosis == "AD"),
            design$n_batches[i] == length(unique(cohort$batch)),
            design$design_columns[i] == ncol(matrix), design$design_rank[i] == qr(matrix)$rank,
            design$design_rank[i] == ncol(matrix),
            design$residual_df[i] == nrow(matrix) - qr(matrix)$rank, design$residual_df[i] > 0,
            "diagnosisAD" %in% colnames(matrix))
  matrix
})

pca <- read.delim("results/pca_scores.tsv")
pca_metadata <- metadata[match(pca$sample_id, metadata$sample_id), ]
pca_fields <- c("sample_id", "diagnosis", "sex", "batch", "bisulfite_batch")
stopifnot(!anyDuplicated(pca$sample_id), nrow(pca) == nrow(metadata),
          !anyNA(pca_metadata$sample_id),
          all(vapply(pca_fields, function(field) {
            identical(pca[[field]], pca_metadata[[field]])
          }, logical(1))))

top <- read.delim("results/differential_top.tsv")
summary <- read.delim("results/differential_summary.tsv")
stopifnot(nrow(top) == min(100L, summary$n_probes), !is.unsorted(top$P.Value))
holdout <- read.delim("results/leave_one_batch_out_summary.tsv")
holdout_effects <- read.delim("results/leave_one_batch_out_effects.tsv")
stopifnot(identical(holdout$omitted_batch, c("none", pairs$batch[pairs$shared])),
          !anyDuplicated(holdout_effects[c("omitted_batch", "CpG")]),
          file.info("results/leave_one_batch_out.pdf")$size > 0L,
          isTRUE(all.equal(holdout_effects$delta_beta[holdout_effects$omitted_batch == "none"], top$delta_beta)))
for (i in seq_len(nrow(holdout))) {
  row <- holdout[i, ]
  cohort <- metadata[metadata$batch != row$omitted_batch, ]
  cohort$diagnosis <- factor(cohort$diagnosis, levels = c("Control", "AD"))
  design_matrix <- model.matrix(as.formula(design$model[1L]), cohort)
  effects <- holdout_effects[holdout_effects$omitted_batch == row$omitted_batch, ]
  stopifnot(row$status == "estimable", row$n_samples == nrow(cohort), row$n_removed == nrow(metadata) - nrow(cohort),
            row$design_rank == qr(design_matrix)$rank, row$design_rank == row$design_columns,
            row$residual_df == nrow(cohort) - row$design_rank, row$residual_df > 0,
            identical(effects$CpG, top$CpG), all(effects$ci_95_low <= effects$delta_beta),
            all(effects$delta_beta <= effects$ci_95_high),
            row$n_leading_same_direction == sum(sign(effects$delta_beta) == sign(top$delta_beta)),
            row$n_primary_hits_passing_both_cutoffs == sum(effects$primary_passed_both_cutoffs &
                                                            effects$fdr < 0.05 & abs(effects$delta_beta) >= 0.05))
}
mask <- read.delim("data/reference/EPIC.hg19.MASK_general.20180909.tsv.gz", colClasses = "character")
mask_summary <- read.delim("results/probe_mask_summary.tsv")
mask_comparison <- read.delim("results/probe_mask_comparison.tsv")
before_top <- which(mask_comparison$rank_before_mask <= 100L)
after_top <- which(mask_comparison$rank_after_mask <= 100L)
stopifnot(!any(top$CpG %in% mask$CpG), mask_summary$probes_after_mask == summary$n_probes,
          mask_summary$probes_before_mask - summary$n_probes == mask_summary$probes_excluded,
          mask_summary$probes_excluded > 0L,
          mask_summary$top100_overlap == length(intersect(before_top, after_top)),
          mask_summary$pre_mask_top100_excluded == sum(mask_comparison$excluded_by_mask[before_top]),
          identical(top$CpG, mask_comparison$CpG[after_top][order(mask_comparison$rank_after_mask[after_top])]))
comparison <- read.delim("results/differential_sensitivity_top.tsv")
sensitivity <- read.delim("results/differential_sensitivity_summary.tsv")
stopifnot(identical(comparison$CpG, top$CpG), identical(comparison$primary_logFC, top$logFC),
          identical(comparison$primary_delta_beta, top$delta_beta),
          identical(comparison$primary_fdr, top$adj.P.Val),
          sensitivity$n_probes == summary$n_probes,
          sensitivity$n_primary_top_probes == nrow(comparison),
          sensitivity$primary_passing_both_cutoffs == summary$significant_fdr_005_delta_beta_005,
          isTRUE(all.equal(sensitivity$top_m_direction_agreement,
                           mean(sign(comparison$primary_logFC) == sign(comparison$shared_logFC)))))
if (file.exists("results/differential_sensitivity.rds")) {
  paired <- readRDS("results/differential_sensitivity.rds")
  stopifnot(nrow(paired) == summary$n_probes, !anyDuplicated(paired$CpG),
            isTRUE(all.equal(sensitivity$m_effect_pearson_r, cor(paired$primary_logFC, paired$shared_logFC))),
            sensitivity$shared_passing_both_cutoffs ==
              sum(paired$shared_fdr < 0.05 & abs(paired$shared_delta_beta) >= 0.05))
  rm(paired)
}
if (file.exists("results/differential_methylation.rds")) {
  full <- readRDS("results/differential_methylation.rds")
  expected <- head(order(full$P.Value, -abs(full$delta_beta), na.last = NA), 100L)
  expected_lambda <- median(qchisq(1 - full$P.Value, df = 1), na.rm = TRUE) /
    qchisq(0.5, df = 1)
  stopifnot(identical(top$CpG, full$CpG[expected]),
            isTRUE(all.equal(summary$genomic_inflation_factor, expected_lambda)))
}

dmr_summary <- read.delim("results/dmr_summary.tsv")
dmr_ranges <- read.delim("results/dmr_ranges.tsv")
if (!is.na(dmr_summary$n_regions)) stopifnot(nrow(dmr_ranges) == dmr_summary$n_regions)
if (identical(dmr_summary$status, "insufficient significant CpGs")) {
  stopifnot(dmr_summary$n_fdr_significant_cpgs < 2L, nrow(dmr_ranges) == 0L)
}

enrichment <- read.delim("results/enrichment_results.tsv")
stopifnot(setequal(unique(enrichment$collection), c("GO", "KEGG")),
          all(nzchar(enrichment$Description, keepNA = TRUE)),
          all(nzchar(enrichment$pathway_id, keepNA = TRUE)))

app <- shiny::shinyAppDir("app")
has_beta_matrix <- file.exists("results/beta_filtered.rds")
if (has_beta_matrix) {
  beta <- readRDS("results/beta_filtered.rds")
  stopifnot(!any(rownames(beta) %in% mask$CpG), nrow(beta) == summary$n_probes)
  ids <- sub("_.*$", "", colnames(beta))
  probes <- head(top$CpG, 5L)
  for (i in seq_along(cohorts)) {
    y <- t(beta[probes, match(cohorts[[i]]$sample_id, ids), drop = FALSE])
    expected <- lm.fit(design_matrices[[i]], y)$coefficients["diagnosisAD", ]
    field <- c("primary_delta_beta", "shared_delta_beta")[i]
    stopifnot(isTRUE(all.equal(unname(expected), head(comparison[[field]], length(probes)))))
  }
  for (omitted in holdout$omitted_batch) {
    cohort <- metadata[metadata$batch != omitted, ]
    cohort$diagnosis <- factor(cohort$diagnosis, levels = c("Control", "AD"))
    cohort$observed_beta <- as.numeric(beta[top$CpG[1L], match(cohort$sample_id, ids)])
    fit <- lm(update(as.formula(design$model[1L]), observed_beta ~ .), cohort)
    expected <- c(coef(fit)["diagnosisAD"], confint(fit, "diagnosisAD"))
    row <- holdout_effects[holdout_effects$omitted_batch == omitted & holdout_effects$primary_rank == 1L, ]
    stopifnot(isTRUE(all.equal(unname(expected), unlist(row[c("delta_beta", "ci_95_low", "ci_95_high")], use.names = FALSE))))
  }
  rm(beta)
}
shiny::testServer(app, {
  session$setInputs(pca_colour = "batch", fdr_cutoff = 1, delta_cutoff = 0)
  stopifnot(grepl("diagnosis", output$overview_counts, fixed = TRUE),
            !is.null(output$overview_qc), !is.null(output$overview_design), !is.null(output$pca_plot),
            nrow(filtered_diff()) == nrow(top), !is.null(output$volcano_plot),
            !is.null(output$diff_table), !is.null(output$sensitivity_summary), !is.null(output$enrichment_plot),
            !is.null(output$enrichment_table))
  if (has_beta_matrix) {
    session$setInputs(probe_query = top$CpG[1L])
    stopifnot(!is.null(output$probe_plot))
  }
})

report_dir <- tempfile("workflow-report-")
dir.create(report_dir)
report <- rmarkdown::render("report/analysis_report.Rmd", output_dir = report_dir,
                            envir = new.env(), quiet = TRUE)
stopifnot(file.info(report)$size > 0L)
unlink(report_dir, recursive = TRUE)
if ("--replication" %in% commandArgs(trailingOnly = TRUE)) {
  replication <- read.delim("results/replication_GSE59685_summary.tsv")
  effects <- read.delim("results/replication_GSE59685_effects.tsv")
  tested <- effects[effects$status == "tested", ]
  stopifnot(identical(top$CpG, effects$CpG), identical(top$delta_beta, effects$discovery_delta_beta),
            replication$n_ad == 60L, replication$n_control == 24L, replication$residual_df > 0L,
            replication$probes_compared == nrow(tested),
            replication$discovery_passing_cutoffs == sum(effects$discovery_passed_cutoffs),
            replication$discovery_passes_tested == sum(tested$discovery_passed_cutoffs),
            isTRUE(all.equal(tested$replication_fdr, p.adjust(tested$replication_p, "BH"))),
            isTRUE(all.equal(replication$direction_agreement, mean(tested$same_direction))),
            isTRUE(all.equal(replication$beta_effect_pearson_r,
                             cor(tested$discovery_delta_beta, tested$replication_delta_beta))),
            file.info("results/replication_GSE59685.pdf")$size > 0L)
}
message("Workflow smoke test passed.")
