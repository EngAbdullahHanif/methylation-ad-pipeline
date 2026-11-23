#!/usr/bin/env Rscript

source("scripts/extension_models.R")
mfg <- read_extension_assays("pathology")
dg <- read_extension_assays("dentate")
n_mfg_qc_pass <- nrow(mfg$metadata)
n_dg_qc_pass <- nrow(dg$metadata)
mfg_eligible <- subset(mfg$metadata, diagnosis %in% c("Control", "AD"))
dg_eligible <- subset(dg$metadata, diagnosis %in% c("Control", "AD"))
fields <- c("donor_id", "sample_id", "diagnosis", "sex", "batch", "proportion_neurons")
pairs <- merge(mfg_eligible[, fields], dg_eligible[, fields], by = "donor_id", suffixes = c("_mfg", "_dg"))
if (!nrow(pairs) || anyNA(pairs) || any(pairs$diagnosis_mfg != pairs$diagnosis_dg | pairs$sex_mfg != pairs$sex_dg)) {
  stop("Expected complete paired donor metadata with concordant diagnosis and sex.")
}
pair_data <- data.frame(donor_id = pairs$donor_id, diagnosis = factor(pairs$diagnosis_mfg, c("Control", "AD")),
                        sex = pairs$sex_mfg, neuronal_difference = pairs$proportion_neurons_dg - pairs$proportion_neurons_mfg)

# Within-donor DG minus MFG removes the donor intercept and its time-invariant main effects.
# A diagnosis coefficient on that difference equals the diagnosis-by-region interaction.
batches <- sort(unique(c(pairs$batch_mfg, pairs$batch_dg)))
batch_matrix <- function(batch) model.matrix(~ batch - 1, data.frame(batch = factor(batch, levels = batches)))
batch_difference <- batch_matrix(pairs$batch_dg) - batch_matrix(pairs$batch_mfg)
nuisance <- cbind(model.matrix(~ sex + neuronal_difference, pair_data), batch_difference)
decomposition <- qr(nuisance)
basis <- sort(decomposition$pivot[seq_len(decomposition$rank)])
design <- cbind(nuisance[, basis, drop = FALSE], diagnosisAD = as.integer(pair_data$diagnosis == "AD"))
if (nlevels(droplevels(pair_data$diagnosis)) != 2L) stop("Expected paired donors in both diagnosis groups.")

common <- intersect(rownames(mfg$m), rownames(dg$m))
mfg_columns <- match(pairs$sample_id_mfg, mfg$metadata$sample_id)
dg_columns <- match(pairs$sample_id_dg, dg$metadata$sample_id)
if (!length(common) || anyNA(c(mfg_columns, dg_columns))) stop("Could not match paired arrays and masked probe IDs.")
m_difference <- dg$m[common, dg_columns, drop = FALSE] - mfg$m[common, mfg_columns, drop = FALSE]
beta_difference <- dg$beta[common, dg_columns, drop = FALSE] - mfg$beta[common, mfg_columns, drop = FALSE]
colnames(m_difference) <- colnames(beta_difference) <- as.character(pairs$donor_id)
rm(mfg, dg)
gc()
result <- fit_extension(m_difference, beta_difference, design, "diagnosisAD")
residual_diagnosis <- qr.resid(decomposition, design[, "diagnosisAD"])
pairs$diagnosis_information_share <- residual_diagnosis^2 / sum(residual_diagnosis^2)

summary <- cbind(data.frame(contrast = "AD minus control difference in DG minus MFG",
                            model = "DG minus MFG ~ diagnosis + sex + neuronal_difference + batch_difference",
                            n_mfg_qc_pass = n_mfg_qc_pass, n_dg_qc_pass = n_dg_qc_pass,
                            n_mfg_diagnosis_eligible = nrow(mfg_eligible), n_dg_diagnosis_eligible = nrow(dg_eligible),
                            n_pairs = nrow(pairs), n_control_pairs = sum(pair_data$diagnosis == "Control"),
                            n_ad_pairs = sum(pair_data$diagnosis == "AD"), n_batches = length(batches),
                            nuisance_rank = decomposition$rank,
                            max_donor_diagnosis_information_share = max(pairs$diagnosis_information_share)),
                 extension_statistics(result, design))
write_extension_table(summary, "regional_specificity_summary")
write_extension_table(head(result, 100), "regional_specificity_top")
write_extension_table(pairs, "regional_specificity_pairs")
saveRDS(result, "results/extensions/regional_site_statistics.rds")
plot_extension_effects(result, "regional_specificity.pdf", "Adjusted delta beta interaction (DG minus MFG), 95% CI")
