#!/usr/bin/env Rscript

# Explore PCs, assess batch structure, and run exploratory differential methylation.

options(stringsAsFactors = FALSE)
set.seed(20261003)

required <- c("limma", "minfi", "matrixStats")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("limma, minfi, and matrixStats are required.")
}

project_dir <- normalizePath(".", mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
metadata_dir <- file.path(project_dir, "data", "metadata")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(file.path(metadata_dir, "sample_sheet_clean.tsv"),
                       check.names = FALSE)
methyl_set <- readRDS(file.path(results_dir, "methyl_set_filtered.rds"))
beta <- readRDS(file.path(results_dir, "beta_filtered.rds"))

beta_sample_ids <- sub("_.*$", "", colnames(beta))
if (anyDuplicated(beta_sample_ids) || anyDuplicated(metadata$sample_id) ||
    !setequal(beta_sample_ids, metadata$sample_id) ||
    !identical(dimnames(beta), dimnames(methyl_set))) {
  stop("Expected one metadata row per beta sample and matching beta and methylation columns.")
}
metadata <- metadata[match(beta_sample_ids, metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AD"))
metadata$sex <- factor(metadata$sex)
metadata$batch <- factor(metadata$batch)
metadata$bisulfite_batch <- factor(metadata$bisulfite_batch)

cell_fields <- c("proportion_neurons", "proportion_olig_opc", "proportion_astrocytes",
                 "proportion_endothelial", "proportion_microglia")
model_fields <- c("diagnosis", "sex", "batch", "bisulfite_batch", cell_fields)
if (!all(model_fields %in% names(metadata)) || anyNA(metadata[, model_fields])) {
  stop("QC-passing samples must have complete diagnosis, sex, batch, and cell proportion data.")
}
cell_proportions <- as.matrix(metadata[, cell_fields])
if (!is.numeric(cell_proportions) || any(!is.finite(cell_proportions)) ||
    any(cell_proportions < 0 | cell_proportions > 1) ||
    any(abs(rowSums(cell_proportions) - 1) > 1e-4)) {
  stop("Expected published cell proportions between 0 and 1 that sum to 1 per sample.")
}

batch_table <- with(metadata, table(batch, diagnosis))
shared_batches <- rownames(batch_table)[rowSums(batch_table > 0L) == 2L]
if (!length(shared_batches)) {
  stop("The QC-passing cohort has no array batch with both diagnosis groups.")
}
batch_overlap <- data.frame(
  batch = rownames(batch_table), n_control = as.integer(batch_table[, "Control"]),
  n_ad = as.integer(batch_table[, "AD"]), shared = rownames(batch_table) %in% shared_batches,
  bisulfite_batches = vapply(split(as.character(metadata$bisulfite_batch), metadata$batch),
                            function(values) paste(sort(unique(values)), collapse = ", "), character(1)),
  row.names = NULL
)
utils::write.table(batch_overlap, file.path(results_dir, "batch_overlap.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

coef_name <- "diagnosisAD"
build_design <- function(included, label) {
  cohort <- droplevels(metadata[included, , drop = FALSE])
  if (nlevels(cohort$diagnosis) != 2L || nlevels(cohort$sex) < 2L) {
    stop(label, " requires both diagnosis groups and variation in recorded sex.")
  }
  terms <- c("diagnosis", "sex", "proportion_neurons")
  if (nlevels(cohort$batch) > 1L) {
    terms <- append(terms, "batch", after = 1L)
  }
  formula <- reformulate(terms)
  design <- model.matrix(formula, data = cohort)
  decomposition <- qr(design)
  residual_df <- nrow(design) - decomposition$rank
  if (decomposition$rank != ncol(design) || residual_df <= 0L ||
      !coef_name %in% colnames(design)) {
    stop(label, " design is not estimable. Rank ", decomposition$rank, " of ", ncol(design),
         "; residual degrees of freedom ", residual_df, ".")
  }
  list(matrix = design, diagnostics = data.frame(
    analysis = label, model = deparse1(formula), coefficient = coef_name,
    n_samples = nrow(cohort), n_control = sum(cohort$diagnosis == "Control"),
    n_ad = sum(cohort$diagnosis == "AD"), n_batches = nlevels(cohort$batch),
    design_columns = ncol(design), design_rank = decomposition$rank, residual_df = residual_df
  ))
}
in_shared_batch <- metadata$batch %in% shared_batches
primary <- build_design(rep(TRUE, nrow(metadata)), "primary")
sensitivity <- build_design(in_shared_batch, "shared_batches")
utils::write.table(rbind(primary$diagnostics, sensitivity$diagnostics),
                   file.path(results_dir, "model_design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

variable_summary <- data.frame(
  variable = c("diagnosis", "sex", "batch", "bisulfite_batch", "age_raw"),
  unique_values = c(
    length(unique(metadata$diagnosis)), length(unique(metadata$sex)),
    length(unique(metadata$batch)), length(unique(metadata$bisulfite_batch)),
    length(unique(metadata$age_raw))
  ),
  stringsAsFactors = FALSE
)
utils::write.table(variable_summary, file.path(results_dir, "m3_variable_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

model_note <- paste("The primary fit uses all QC-passing samples with adjustment for array batch,",
                    "recorded sex and estimated neuronal proportion. The sensitivity fit is restricted",
                    "to array batches containing both diagnoses. Age is excluded because GEO reports",
                    "90+ throughout. Bisulfite batch is recorded separately from the fitted array batch.")
writeLines(model_note, file.path(results_dir, "m3_model_note.txt"))

# PCA uses the most variable probes to make the exploratory plot computationally light.
probe_variance <- matrixStats::rowVars(beta)
variable_probes <- which(is.finite(probe_variance) & probe_variance > 0)
top_probes <- head(variable_probes[order(probe_variance[variable_probes], decreasing = TRUE)], 10000L)
pca <- prcomp(t(beta[top_probes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_scores <- data.frame(
  sample_id = metadata$sample_id,
  PC1 = pca$x[, "PC1"],
  PC2 = pca$x[, "PC2"],
  PC3 = pca$x[, "PC3"],
  metadata[, c("diagnosis", "sex", "batch", "bisulfite_batch"), drop = FALSE],
  row.names = NULL,
  stringsAsFactors = FALSE
)
explained <- (pca$sdev^2) / sum(pca$sdev^2)
pca_variance <- data.frame(
  component = paste0("PC", seq_along(explained)),
  variance_fraction = explained,
  stringsAsFactors = FALSE
)
utils::write.table(pca_scores, file.path(results_dir, "pca_scores.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(pca_variance, file.path(results_dir, "pca_variance.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(results_dir, "pca_by_diagnosis.pdf"), width = 7, height = 6)
plot(pca_scores$PC1, pca_scores$PC2, col = as.integer(pca_scores$diagnosis),
     pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA colored by diagnosis")
text(pca_scores$PC1, pca_scores$PC2, labels = pca_scores$sample_id, pos = 3, cex = 0.6)
legend("topright", legend = levels(metadata$diagnosis), col = seq_along(levels(metadata$diagnosis)), pch = 19)
dev.off()

pdf(file.path(results_dir, "pca_by_batch.pdf"), width = 8, height = 6)
par(mar = c(5, 4, 4, 8) + 0.1)
batch_colors <- grDevices::hcl.colors(nlevels(metadata$batch), "Dark 3")
plot(pca_scores$PC1, pca_scores$PC2, col = batch_colors[as.integer(pca_scores$batch)],
     pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA colored by batch")
legend(par("usr")[2], par("usr")[4], legend = levels(metadata$batch), col = batch_colors,
       pch = 19, cex = 0.8, xpd = TRUE, bty = "n")
dev.off()

# Test continuous PC associations only where the metadata varies.
pc_tests <- do.call(rbind, lapply(seq_len(min(5L, ncol(pca$x))), function(index) {
  pc_name <- paste0("PC", index)
  pc_values <- pca$x[, index]
  variables <- c("diagnosis", "sex", "batch", "bisulfite_batch")
  do.call(rbind, lapply(variables, function(variable) {
    values <- pca_scores[[variable]]
    if (length(unique(values)) < 2L) {
      return(data.frame(pc = pc_name, variable = variable, p_value = NA_real_, note = "constant", stringsAsFactors = FALSE))
    }
    fit <- lm(pc_values ~ values)
    data.frame(pc = pc_name, variable = variable,
               p_value = anova(fit)[1, "Pr(>F)"], note = "descriptive association", stringsAsFactors = FALSE)
  }))
}))
pc_tests$FDR <- p.adjust(pc_tests$p_value, method = "BH")
utils::write.table(pc_tests, file.path(results_dir, "pc_variable_associations.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# Fit the same contrast on each scale. Limma tests the M-value coefficient;
# the beta-scale coefficient uses the same covariate adjustment for interpretation.
fit_diagnosis <- function(m, beta, design) {
  fit <- limma::eBayes(limma::lmFit(m, design), robust = TRUE)
  result <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  weights <- qr.solve(design, diag(nrow(design)))[match(coef_name, colnames(design)), ]
  result$CpG <- rownames(result)
  result$delta_beta <- drop(beta %*% weights)[match(result$CpG, rownames(beta))]
  result
}

m_values <- minfi::getM(methyl_set)
rm(methyl_set)
if (any(!is.finite(m_values))) {
  stop("The filtered methylation data contain non-finite M-values.")
}
results <- fit_diagnosis(m_values, beta, primary$matrix)
shared_results <- fit_diagnosis(m_values[, in_shared_batch, drop = FALSE],
                                beta[, in_shared_batch, drop = FALSE], sensitivity$matrix)
rm(beta)
results$gene <- NA_character_
results <- results[, c("CpG", "gene", "logFC", "delta_beta", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
utils::write.table(results, file.path(results_dir, "differential_methylation.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(results, file.path(results_dir, "differential_methylation.rds"))
top_results <- head(results[order(results$P.Value, -abs(results$delta_beta),
                                  na.last = NA), ], 100L)
utils::write.table(top_results, file.path(results_dir, "differential_top.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

if (!identical(results$CpG, shared_results$CpG)) {
  stop("The primary and shared-batch fits returned different probe orders.")
}
effect_comparison <- data.frame(
  CpG = results$CpG, primary_logFC = results$logFC, shared_logFC = shared_results$logFC,
  primary_delta_beta = results$delta_beta, shared_delta_beta = shared_results$delta_beta,
  primary_fdr = results$adj.P.Val, shared_fdr = shared_results$adj.P.Val
)
top_comparison <- effect_comparison[match(top_results$CpG, effect_comparison$CpG), ]
significant <- function(result) {
  result$adj.P.Val < 0.05 & abs(result$delta_beta) >= 0.05
}
direction_agreement <- function(first, second) mean(sign(first) == sign(second))
comparison_summary <- data.frame(
  n_probes = nrow(effect_comparison),
  m_effect_pearson_r = cor(results$logFC, shared_results$logFC),
  beta_effect_pearson_r = cor(results$delta_beta, shared_results$delta_beta),
  m_direction_agreement = direction_agreement(results$logFC, shared_results$logFC),
  n_primary_top_probes = nrow(top_comparison),
  top_m_direction_agreement = direction_agreement(top_comparison$primary_logFC, top_comparison$shared_logFC),
  median_abs_delta_beta_change = median(abs(results$delta_beta - shared_results$delta_beta)),
  primary_passing_both_cutoffs = sum(significant(results)),
  shared_passing_both_cutoffs = sum(significant(shared_results)),
  passing_both_cutoffs_in_both = sum(significant(results) & significant(shared_results))
)
utils::write.table(comparison_summary, file.path(results_dir, "differential_sensitivity_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(top_comparison, file.path(results_dir, "differential_sensitivity_top.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(effect_comparison, file.path(results_dir, "differential_sensitivity.rds"))
rm(shared_results, effect_comparison)

pdf(file.path(results_dir, "differential_volcano.pdf"), width = 7, height = 6)
plot_data <- is.finite(results$delta_beta) & is.finite(results$P.Value)
smoothScatter(results$delta_beta[plot_data], -log10(pmax(results$P.Value[plot_data], .Machine$double.xmin)),
              xlab = "Adjusted delta beta (AD - Control)", ylab = "-log10(p-value)",
              main = "Exploratory differential methylation")
abline(v = 0, lty = 2, col = "grey50")
dev.off()

pdf(file.path(results_dir, "differential_qq.pdf"), width = 6, height = 6)
observed <- sort(results$P.Value[is.finite(results$P.Value)])
# Thin the plotted ranks on a log scale while retaining the smallest p-values.
plot_ranks <- unique(round(10^seq(0, log10(length(observed)), length.out = 10000L)))
plot(-log10(ppoints(length(observed))[plot_ranks]),
     -log10(pmax(observed[plot_ranks], .Machine$double.xmin)), pch = 19, cex = 0.35,
     xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)", main = "Primary model QQ plot")
abline(0, 1, lty = 2, col = "grey50")
dev.off()

lambda <- median(qchisq(results$P.Value, df = 1, lower.tail = FALSE), na.rm = TRUE) /
  qchisq(0.5, df = 1)
significant_count <- comparison_summary$primary_passing_both_cutoffs
differential_summary <- data.frame(
  model = primary$diagnostics$model,
  coefficient = coef_name,
  n_samples = nrow(metadata),
  n_probes = nrow(results),
  significant_fdr_005_delta_beta_005 = significant_count,
  genomic_inflation_factor = lambda,
  model_note = model_note,
  stringsAsFactors = FALSE
)
utils::write.table(differential_summary, file.path(results_dir, "differential_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

dmr_summary <- data.frame(
  method = "DMRcate",
  status = "not run",
  n_regions = NA_integer_,
  n_fdr_significant_cpgs = NA_integer_,
  note = "",
  stringsAsFactors = FALSE
)
dmr_ranges <- GenomicRanges::GRanges()
if (requireNamespace("DMRcate", quietly = TRUE)) {
  cpg_annotated <- DMRcate::cpg.annotate(
    datatype = "array",
    object = m_values,
    what = "M",
    arraytype = "EPICv1",
    analysis.type = "differential",
    design = primary$matrix,
    coef = match(coef_name, colnames(primary$matrix)),
    fdr = 0.05
  )
  dmr_summary$n_fdr_significant_cpgs <- sum(cpg_annotated@ranges$is.sig, na.rm = TRUE)
  if (dmr_summary$n_fdr_significant_cpgs >= 2L) {
    dmr_output <- DMRcate::dmrcate(cpg_annotated, lambda = 1000, C = 2,
                                  pcutoff = "fdr", min.cpgs = 2L)
    dmr_ranges <- DMRcate::extractRanges(dmr_output, genome = "hg19")
    dmr_summary$status <- "complete"
    dmr_summary$note <- "Primary covariate design, EPICv1, lambda 1000 bp, minimum two CpGs, FDR cutoff."
  } else {
    dmr_summary$status <- "insufficient significant CpGs"
    dmr_summary$note <- "Fewer than two individually significant CpGs in DMRcate at FDR 0.05."
  }
  dmr_summary$n_regions <- length(dmr_ranges)
} else {
  dmr_summary$note <- "DMRcate was unavailable at runtime."
}
saveRDS(dmr_ranges, file.path(results_dir, "dmr_ranges.rds"))
dmr_table <- as.data.frame(dmr_ranges)
if (!"no.cpgs" %in% names(dmr_table)) {
  dmr_table <- cbind(dmr_table, data.frame(
    no.cpgs = integer(), min_smoothed_fdr = numeric(), Stouffer = numeric(), HMFDR = numeric(),
    Fisher = numeric(), maxdiff = numeric(), meandiff = numeric(), overlapping.genes = character()
  ))
}
utils::write.table(dmr_table, file.path(results_dir, "dmr_ranges.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(dmr_summary, file.path(results_dir, "dmr_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(list(pca_scores = pca_scores, pca_variance = pca_variance,
             differential = results, summary = differential_summary,
             model_design = rbind(primary$diagnostics, sensitivity$diagnostics),
             sensitivity_summary = comparison_summary,
             dmr_summary = dmr_summary),
        file.path(results_dir, "m3_results.rds"))
message("M3 PCA and limma analysis complete. Significant probes: ", significant_count)
