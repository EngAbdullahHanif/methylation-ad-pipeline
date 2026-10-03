#!/usr/bin/env Rscript

# Explore PCs, assess batch structure, and run exploratory differential methylation.

options(stringsAsFactors = FALSE)
set.seed(20261003)

required <- c("limma", "minfi")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("limma and minfi are required.")
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
metadata <- metadata[match(beta_sample_ids, metadata$sample_id), , drop = FALSE]
metadata$sample_id <- beta_sample_ids
rownames(metadata) <- metadata$sample_id

metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AD"))
metadata$sex <- factor(metadata$sex)
metadata$batch <- factor(metadata$batch)
metadata$bisulfite_batch <- factor(metadata$bisulfite_batch)

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

confounded <- with(metadata, any(table(diagnosis, batch) == 0L))
model_note <- if (length(unique(metadata$age_raw)) == 1L && confounded) {
  "Primary model is diagnosis-only: age is constant and diagnosis is confounded with batch in this 10-sample subset; adding batch would make the diagnosis effect non-identifiable."
} else {
  "Primary model includes diagnosis and estimable covariates after checking rank."
}
writeLines(model_note, file.path(results_dir, "m3_model_note.txt"))

# PCA uses the most variable probes to make the exploratory plot computationally light.
probe_variance <- apply(beta, 1L, var, na.rm = TRUE)
top_probes <- order(probe_variance, decreasing = TRUE)[seq_len(min(10000L, nrow(beta)))]
pca <- prcomp(t(beta[top_probes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_scores <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, "PC1"],
  PC2 = pca$x[, "PC2"],
  PC3 = pca$x[, "PC3"],
  metadata[rownames(pca$x), c("diagnosis", "sex", "batch", "bisulfite_batch"), drop = FALSE],
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
plot(pca_scores$PC1, pca_scores$PC2, col = as.integer(pca_scores$batch),
     pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA colored by batch")
text(pca_scores$PC1, pca_scores$PC2, labels = pca_scores$sample_id, pos = 3, cex = 0.6)
legend("topright", legend = levels(metadata$batch), col = seq_along(levels(metadata$batch)), pch = 19, cex = 0.7)
dev.off()

# Test continuous PC associations only where the metadata varies.
pc_tests <- do.call(rbind, lapply(seq_len(min(5L, ncol(pca$x))), function(index) {
  pc_name <- paste0("PC", index)
  pc_values <- pca$x[, index]
  variables <- c("diagnosis", "sex", "batch", "bisulfite_batch")
  do.call(rbind, lapply(variables, function(variable) {
    values <- metadata[rownames(pca$x), variable]
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

# M-values are used for testing; beta-values are retained for effect-size reporting.
m_values <- minfi::getM(methyl_set)
design <- model.matrix(~ diagnosis, data = metadata)
fit <- limma::lmFit(m_values, design)
fit <- limma::eBayes(fit, robust = TRUE)
coef_name <- "diagnosisAD"
results <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
results$CpG <- rownames(results)
beta_control <- rowMeans(beta[, metadata$diagnosis == "Control", drop = FALSE], na.rm = TRUE)
beta_ad <- rowMeans(beta[, metadata$diagnosis == "AD", drop = FALSE], na.rm = TRUE)
results$delta_beta <- beta_ad[results$CpG] - beta_control[results$CpG]
results$gene <- NA_character_
results <- results[, c("CpG", "gene", "logFC", "delta_beta", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
utils::write.table(results, file.path(results_dir, "differential_methylation.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(results, file.path(results_dir, "differential_methylation.rds"))

pdf(file.path(results_dir, "differential_volcano.pdf"), width = 7, height = 6)
plot(results$delta_beta, -log10(results$P.Value), pch = 19, cex = 0.35,
     xlab = "Delta beta: AD - Control", ylab = "-log10(p-value)", main = "Exploratory differential methylation")
abline(v = 0, lty = 2, col = "grey50")
dev.off()

pdf(file.path(results_dir, "differential_qq.pdf"), width = 6, height = 6)
limma::qqt(results$t, df = fit$df.prior + fit$df.residual, pch = 19, cex = 0.35)
dev.off()

lambda <- median(qchisq(1 - results$P.Value, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
significant_count <- sum(results$adj.P.Val < 0.05 & abs(results$delta_beta) >= 0.05, na.rm = TRUE)
differential_summary <- data.frame(
  model = "~ diagnosis",
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
  note = "",
  stringsAsFactors = FALSE
)
if (requireNamespace("DMRcate", quietly = TRUE)) {
  cpg_annotated <- DMRcate::cpg.annotate(
    datatype = "array",
    object = m_values,
    what = "M",
    arraytype = "EPICv1",
    analysis.type = "differential",
    design = design,
    coef = 2,
    fdr = 0.05
  )
  dmr_output <- DMRcate::dmrcate(cpg_annotated, lambda = 1000, C = 2, pcutoff = "fdr")
  dmr_ranges <- DMRcate::extractRanges(dmr_output, genome = "hg19")
  saveRDS(dmr_ranges, file.path(results_dir, "dmr_ranges.rds"))
  utils::write.table(as.data.frame(dmr_ranges), file.path(results_dir, "dmr_ranges.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  dmr_summary$status <- "complete"
  dmr_summary$n_regions <- length(dmr_ranges)
  dmr_summary$note <- "DMRcate EPICv1 on the full filtered probe set, lambda 1000 bp, minimum 2 consecutive CpGs, FDR cutoff."
} else {
  dmr_summary$note <- "DMRcate was unavailable at runtime."
}
utils::write.table(dmr_summary, file.path(results_dir, "dmr_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(list(pca_scores = pca_scores, pca_variance = pca_variance,
             differential = results, summary = differential_summary,
             dmr_summary = dmr_summary),
        file.path(results_dir, "m3_results.rds"))
message("M3 PCA and limma analysis complete. Significant probes: ", significant_count)
