#!/usr/bin/env Rscript

# Compare the masked discovery candidates with an independent frontal cortex cohort.
discovery_path <- "results/differential_top.tsv"
discovery <- read.delim(discovery_path, check.names = FALSE)
excluded <- read.delim("data/reference/EPIC.hg19.MASK_general.20180909.tsv.gz", colClasses = "character")$CpG
if (any(discovery$CpG %in% excluded)) stop("Regenerate the masked discovery results before comparing cohorts.")
prepared <- tempfile("GSE59685-")
status <- system2("python3", c("scripts/prepare_replication.py", shQuote(discovery_path), shQuote(prepared)))
if (status != 0L) stop("GSE59685 input preparation failed.")

samples <- read.delim(file.path(prepared, "samples.tsv"), colClasses = c(chip = "character"))
measurements <- read.delim(file.path(prepared, "beta.tsv"), check.names = FALSE,
                           na.strings = c("", "NA", "NaN"), colClasses = c(CpG = "character"))
unlink(prepared, recursive = TRUE)
if (!identical(names(measurements)[-1L], samples$sample_id) || anyDuplicated(samples$donor_id)) {
  stop("Replication beta columns must match one sample per donor in the cohort table.")
}
samples$diagnosis <- factor(samples$diagnosis, levels = c("Control", "AD"))
samples$sex <- factor(samples$sex)
samples$chip <- factor(samples$chip)
if (anyNA(samples[, c("diagnosis", "age", "sex", "chip")])) stop("Incomplete GSE59685 model metadata.")
design <- model.matrix(~ diagnosis + age + sex + chip, samples)
rank <- qr(design)$rank
residual_df <- nrow(design) - rank
if (rank < ncol(design) || residual_df <= 0L || !"diagnosisAD" %in% colnames(design)) {
  stop("The GSE59685 diagnosis contrast is not estimable with age, sex and chip adjustment.")
}

beta <- as.matrix(measurements[, -1L, drop = FALSE])
rownames(beta) <- measurements$CpG
complete <- rowSums(is.finite(beta) & beta >= 0 & beta <= 1) == ncol(beta)
# A constant beta row has an undefined t statistic.
variances <- matrixStats::rowVars(beta, na.rm = TRUE)
tested <- complete & is.finite(variances) & variances > 0
if (sum(tested) < 2L) stop("Fewer than two candidate probes had complete, variable beta measurements.")
fit <- limma::lmFit(beta[tested, , drop = FALSE], design)
effect <- fit$coefficients[, "diagnosisAD"]
se <- fit$stdev.unscaled[, "diagnosisAD"] * fit$sigma
if (any(!is.finite(se) | se <= 0)) stop("The replication fit returned invalid standard errors.")
critical <- qt(0.975, df = residual_df)
w <- match(discovery$CpG, names(effect))
p <- 2 * pt(-abs(effect / se), df = residual_df)
comparison <- data.frame(
  discovery_rank = seq_len(nrow(discovery)), CpG = discovery$CpG,
  discovery_delta_beta = discovery$delta_beta, discovery_fdr = discovery$adj.P.Val,
  discovery_passed_cutoffs = discovery$adj.P.Val < 0.05 & abs(discovery$delta_beta) >= 0.05,
  status = ifelse(is.na(w), ifelse(discovery$CpG %in% rownames(beta),
                                  "incomplete_or_constant", "not_in_published_matrix"), "tested"),
  replication_delta_beta = unname(effect[w]),
  replication_ci_95_low = unname((effect - critical * se)[w]),
  replication_ci_95_high = unname((effect + critical * se)[w]),
  replication_p = unname(p[w]), replication_fdr = unname(p.adjust(p, "BH")[w])
)
comparison$same_direction <- sign(comparison$discovery_delta_beta) == sign(comparison$replication_delta_beta)
compared <- comparison[comparison$status == "tested", ]
summary <- data.frame(
  cohort = "GSE59685", region = "frontal cortex", platform = "Illumina 450K",
  model = "beta ~ diagnosis + age + sex + chip", n_control = sum(samples$diagnosis == "Control"),
  n_ad = sum(samples$diagnosis == "AD"), design_rank = rank, residual_df = residual_df,
  discovery_candidates = nrow(discovery), probes_in_matrix = nrow(beta), probes_compared = nrow(compared),
  beta_effect_pearson_r = cor(compared$discovery_delta_beta, compared$replication_delta_beta),
  direction_agreement = mean(compared$same_direction),
  median_abs_discovery_delta_beta = median(abs(compared$discovery_delta_beta)),
  median_abs_replication_delta_beta = median(abs(compared$replication_delta_beta)),
  median_abs_effect_difference = median(abs(compared$discovery_delta_beta - compared$replication_delta_beta)),
  same_direction_fdr_below_005 = sum(compared$same_direction & compared$replication_fdr < 0.05),
  discovery_passing_cutoffs = sum(comparison$discovery_passed_cutoffs),
  discovery_passes_tested = sum(compared$discovery_passed_cutoffs),
  discovery_passes_supported = sum(compared$discovery_passed_cutoffs & compared$same_direction &
                                     compared$replication_fdr < 0.05)
)
write.table(summary, "results/replication_GSE59685_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comparison, "results/replication_GSE59685_effects.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("results/replication_GSE59685.pdf", width = 7, height = 6)
plot_range <- range(c(compared$discovery_delta_beta, compared$replication_delta_beta))
plot(compared$discovery_delta_beta, compared$replication_delta_beta, pch = 19, col = "#286A79",
     xlim = plot_range, ylim = plot_range, asp = 1,
     xlab = "Discovery adjusted delta beta (AD - Control)",
     ylab = "GSE59685 adjusted delta beta (AD - Control)", main = "Independent frontal cortex comparison")
abline(h = 0, v = 0, lty = 3, col = "grey60")
abline(0, 1, lty = 2, col = "grey60")
dev.off()
