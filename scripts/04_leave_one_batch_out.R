#!/usr/bin/env Rscript

# Measure how the primary methylation contrast changes when a shared batch is omitted.
metadata <- read.delim("data/metadata/sample_sheet_clean.tsv")
model <- subset(read.delim("results/model_design.tsv"), analysis == "primary")
overlap <- read.delim("results/batch_overlap.tsv")
leading <- read.delim("results/differential_top.tsv")
reference <- readRDS("results/differential_methylation.rds")
beta <- readRDS("results/beta_filtered.rds")
m_values <- minfi::getM(readRDS("results/methyl_set_filtered.rds"))
sample_ids <- sub("_.*$", "", colnames(beta))
if (nrow(model) != 1L || anyDuplicated(sample_ids) || anyDuplicated(metadata$sample_id) ||
    !setequal(sample_ids, metadata$sample_id) || !identical(dimnames(beta), dimnames(m_values)) ||
    !identical(rownames(beta), reference$CpG) || !all(leading$CpG %in% reference$CpG)) {
  stop("Expected the same QC samples and masked probes as the primary fit.")
}
metadata <- metadata[match(sample_ids, metadata$sample_id), ]
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AD"))
metadata$batch <- factor(metadata$batch)
metadata$sex <- factor(metadata$sex)
shared <- with(metadata, table(batch, diagnosis))
shared_batches <- rownames(shared)[rowSums(shared > 0L) == 2L]
if (!length(shared_batches) || !setequal(shared_batches, overlap$batch[overlap$shared])) {
  stop("Expected at least one shared batch and matching batch-overlap diagnostics.")
}

formula <- as.formula(model$model)
coefficient <- model$coefficient
leading_index <- match(leading$CpG, reference$CpG)
primary_hits <- leading$adj.P.Val < 0.05 & abs(leading$delta_beta) >= 0.05

fit_without <- function(batch) {
  keep <- if (batch == "none") rep(TRUE, nrow(metadata)) else metadata$batch != batch
  cohort <- droplevels(metadata[keep, , drop = FALSE])
  diagnostics <- data.frame(
    omitted_batch = batch, n_samples = nrow(cohort), n_removed = sum(!keep),
    n_control = sum(cohort$diagnosis == "Control"), n_ad = sum(cohort$diagnosis == "AD"),
    n_shared_batches = sum(shared_batches != batch), design_columns = NA_integer_, design_rank = NA_integer_,
    residual_df = NA_integer_, status = "estimable", m_effect_pearson_r = NA_real_,
    n_leading_same_direction = NA_integer_, leading_median_abs_delta_beta_change = NA_real_,
    n_primary_hits_same_direction = NA_integer_, n_primary_hits_passing_both_cutoffs = NA_integer_,
    genomewide_fdr_005 = NA_integer_
  )
  effects <- data.frame(omitted_batch = batch, primary_rank = seq_len(nrow(leading)), CpG = leading$CpG,
                        primary_passed_both_cutoffs = primary_hits, logFC = NA_real_, delta_beta = NA_real_,
                        ci_95_low = NA_real_, ci_95_high = NA_real_, fdr = NA_real_, same_direction = NA)
  if (any(vapply(cohort[c("diagnosis", "batch", "sex")], nlevels, integer(1)) < 2L)) {
    diagnostics$status <- "insufficient factor levels"
    return(list(diagnostics = diagnostics, effects = effects))
  }
  design <- model.matrix(formula, cohort, na.action = na.fail)
  rank <- qr(design)$rank
  residual_df <- nrow(design) - rank
  diagnostics$design_columns <- ncol(design)
  diagnostics$design_rank <- rank
  diagnostics$residual_df <- residual_df
  if (rank < ncol(design) || residual_df <= 0L || !coefficient %in% colnames(design)) {
    diagnostics$status <- "non-estimable design"
    return(list(diagnostics = diagnostics, effects = effects))
  }

  beta_fit <- limma::lmFit(beta[leading$CpG, keep, drop = FALSE], design)
  effects$delta_beta <- unname(beta_fit$coefficients[, coefficient])
  margin <- qt(0.975, beta_fit$df.residual) * beta_fit$stdev.unscaled[, coefficient] * beta_fit$sigma
  effects$ci_95_low <- effects$delta_beta - margin
  effects$ci_95_high <- effects$delta_beta + margin
  if (batch == "none") {
    if (!isTRUE(all.equal(effects$delta_beta, leading$delta_beta, tolerance = 1e-10))) {
      stop("Recomputed full-cohort beta estimates differ from the primary fit.")
    }
    logfc <- reference$logFC
    fdr <- reference$adj.P.Val
  } else {
    m_fit <- limma::eBayes(limma::lmFit(m_values[, keep, drop = FALSE], design), robust = TRUE)
    logfc <- m_fit$coefficients[, coefficient]
    fdr <- p.adjust(m_fit$p.value[, coefficient], "BH")
  }

  effects$logFC <- logfc[leading_index]
  effects$fdr <- fdr[leading_index]
  effects$same_direction <- sign(effects$delta_beta) == sign(leading$delta_beta)
  diagnostics$m_effect_pearson_r <- cor(reference$logFC, logfc)
  diagnostics$n_leading_same_direction <- sum(effects$same_direction)
  diagnostics$leading_median_abs_delta_beta_change <- median(abs(effects$delta_beta - leading$delta_beta))
  diagnostics$n_primary_hits_same_direction <- sum(primary_hits & effects$same_direction)
  diagnostics$n_primary_hits_passing_both_cutoffs <- sum(primary_hits & effects$fdr < 0.05 & abs(effects$delta_beta) >= 0.05)
  diagnostics$genomewide_fdr_005 <- sum(fdr < 0.05)
  list(diagnostics = diagnostics, effects = effects)
}

fits <- lapply(c("none", shared_batches), fit_without)
diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))
effects <- do.call(rbind, lapply(fits, `[[`, "effects"))
if (diagnostics$status[1L] != "estimable") stop("Could not reconstruct the primary design.")
write.table(diagnostics, "results/leave_one_batch_out_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(effects, "results/leave_one_batch_out_effects.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- effects[effects$primary_rank <= 4L, ]
limits <- range(c(0, plot_data$ci_95_low, plot_data$ci_95_high), finite = TRUE)
labels <- c("All batches", paste("Omit", sub("_", " ", shared_batches)))
pdf("results/leave_one_batch_out.pdf", width = 9, height = 6.5)
par(mfrow = c(2, 2), mar = c(3.2, 7, 2, 1), oma = c(1.5, 0, 0, 0))
for (probe in unique(plot_data$CpG)) {
  rows <- plot_data[plot_data$CpG == probe, ]
  y <- rev(seq_len(nrow(rows)))
  plot(NA_real_, NA_real_, xlim = limits, ylim = range(y) + c(-0.45, 0.45), yaxt = "n",
       xlab = "", ylab = "", main = probe)
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
  abline(v = 0, col = "grey65", lty = 2)
  colors <- ifelse(rows$omitted_batch == "none", "#153F4C", "#287C82")
  segments(rows$ci_95_low, y, rows$ci_95_high, y, col = colors, lwd = 2)
  points(rows$delta_beta, y, pch = 19, col = colors)
  unavailable <- is.na(rows$delta_beta)
  if (any(unavailable)) text(0, y[unavailable], labels = "Not estimable", cex = 0.75, pos = 4)
}
mtext("Adjusted delta beta (AD minus control), 95% confidence intervals", side = 1, outer = TRUE)
dev.off()
