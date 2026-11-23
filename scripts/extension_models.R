read_extension_assays <- function(cohort) {
  directory <- file.path("results", "extensions", cohort)
  metadata <- read.delim(file.path(directory, "sample_sheet_clean.tsv"))
  beta <- readRDS(file.path(directory, "beta_filtered.rds"))
  m <- minfi::getM(readRDS(file.path(directory, "methyl_set_filtered.rds")))
  ids <- sub("_.*$", "", colnames(beta))
  if (anyDuplicated(ids) || anyDuplicated(metadata$sample_id) || anyDuplicated(metadata$donor_id) ||
      !setequal(ids, metadata$sample_id) || !identical(dimnames(m), dimnames(beta))) {
    stop("Expected matching assays, GEO sample IDs, and unique donor IDs in ", cohort, ".")
  }
  if (any(!is.finite(m)) || any(!is.finite(beta))) stop("Non-finite methylation values in ", cohort, ".")
  list(metadata = metadata[match(ids, metadata$sample_id), ], m = m, beta = beta)
}

fit_extension <- function(m, beta, design, coefficient) {
  if (!identical(dimnames(m), dimnames(beta)) || ncol(beta) != nrow(design)) {
    stop("The assay columns must align with the model design.")
  }
  rank <- qr(design)$rank
  if (rank != ncol(design) || nrow(design) <= rank || !coefficient %in% colnames(design)) {
    stop("The specified contrast requires a full-rank design and positive residual degrees of freedom.")
  }
  fit <- limma::eBayes(limma::lmFit(m, design), robust = TRUE)
  p <- fit$p.value[, coefficient]
  beta_fit <- limma::lmFit(beta, design)
  delta <- beta_fit$coefficients[, coefficient]
  margin <- qt(0.975, beta_fit$df.residual) * beta_fit$stdev.unscaled[, coefficient] * beta_fit$sigma
  result <- data.frame(CpG = rownames(m), logFC = fit$coefficients[, coefficient],
                       delta_beta = delta, ci_95_low = delta - margin, ci_95_high = delta + margin,
                       p_value = p, fdr = p.adjust(p, "BH"), row.names = NULL)
  result[order(result$p_value, -abs(result$delta_beta)), ]
}

extension_statistics <- function(result, design) {
  data.frame(design_columns = ncol(design), design_rank = qr(design)$rank,
             residual_df = nrow(design) - qr(design)$rank, n_probes = nrow(result),
             n_fdr_005 = sum(result$fdr < 0.05),
             n_fdr_005_abs_delta_beta_005 = sum(result$fdr < 0.05 & abs(result$delta_beta) >= 0.05),
             genomic_inflation_factor = median(qchisq(result$p_value, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1))
}

write_extension_table <- function(table, name) {
  write.table(table, file.path("results", paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

plot_extension_effects <- function(result, filename, label) {
  top <- head(result, 10L)
  y <- rev(seq_len(nrow(top)))
  pdf(file.path("results", filename), width = 8, height = 5)
  par(mar = c(5, 8, 2, 2))
  plot(top$delta_beta, y, type = "n", xlim = range(c(0, top$ci_95_low, top$ci_95_high)),
       yaxt = "n", xlab = label, ylab = "")
  axis(2, at = y, labels = top$CpG, las = 1)
  abline(v = 0, col = "grey65", lty = 2)
  segments(top$ci_95_low, y, top$ci_95_high, y, col = "#287C82", lwd = 2)
  points(top$delta_beta, y, pch = 19, col = "#153F4C")
  dev.off()
}
