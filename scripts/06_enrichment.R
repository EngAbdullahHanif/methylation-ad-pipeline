#!/usr/bin/env Rscript

# Run exploratory methylation-gene-set enrichment with missMethyl.

options(stringsAsFactors = FALSE)
set.seed(20261003)

if (!requireNamespace("missMethyl", quietly = TRUE) ||
    !requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  stop("missMethyl and the EPIC annotation package are required.")
}
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

project_dir <- normalizePath(".", mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

differential <- read.delim(
  file.path(results_dir, "differential_methylation.tsv"),
  check.names = FALSE
)
all_cpg <- differential$CpG
ranked <- differential[order(differential$P.Value, -abs(differential$delta_beta), na.last = NA), ]
ranked <- ranked[!duplicated(ranked$CpG), ]
ranked_n <- min(100L, nrow(ranked))
sig_cpg <- ranked$CpG[seq_len(ranked_n)]

run_gometh <- function(collection_name) {
  tryCatch(
    gometh(
      sig.cpg = sig_cpg,
      all.cpg = all_cpg,
      collection = collection_name,
      array.type = "EPIC",
      plot.bias = FALSE,
      prior.prob = TRUE
    ),
    error = function(error) {
      data.frame(
        Pathway = NA_character_,
        P.DE = NA_real_,
        FDR = NA_real_,
        error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    }
  )
}

go_results <- run_gometh("GO")
kegg_results <- run_gometh("KEGG")

if (nrow(go_results) > 0L) {
  go_results$collection <- "GO"
}
if (nrow(kegg_results) > 0L) {
  kegg_results$collection <- "KEGG"
}
combine_results <- function(first, second) {
  columns <- union(names(first), names(second))
  for (column in setdiff(columns, names(first))) {
    first[[column]] <- NA
  }
  for (column in setdiff(columns, names(second))) {
    second[[column]] <- NA
  }
  rbind(first[, columns, drop = FALSE], second[, columns, drop = FALSE])
}
enrichment <- combine_results(go_results, kegg_results)
if (nrow(enrichment) > 0L && "FDR" %in% names(enrichment)) {
  enrichment <- enrichment[order(enrichment$FDR, enrichment$P.DE), ]
}

utils::write.table(
  enrichment,
  file.path(results_dir, "enrichment_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(enrichment, file.path(results_dir, "enrichment_results.rds"))

summary <- data.frame(
  input_cpgs = length(all_cpg),
  ranked_cpgs_used = length(sig_cpg),
  ranking = "ascending limma p-value, then descending absolute delta-beta",
  method = "missMethyl gometh with probe-number bias correction",
  collections = "GO and KEGG",
  interpretation = "Exploratory only; the ranked set is used because the strict significant set is small.",
  stringsAsFactors = FALSE
)
utils::write.table(summary, file.path(results_dir, "enrichment_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- enrichment[!is.na(enrichment$Pathway), , drop = FALSE]
plot_data <- plot_data[order(plot_data$FDR, plot_data$P.DE), ]
plot_data <- head(plot_data, 10L)
pdf(file.path(results_dir, "enrichment_top_pathways.pdf"), width = 9, height = 6)
if (nrow(plot_data) > 0L && "Pathway" %in% names(plot_data)) {
  labels <- plot_data$Pathway
  labels <- abbreviate(labels, minlength = 24)
  barplot(-log10(pmax(plot_data$FDR, .Machine$double.xmin)),
          names.arg = labels, las = 2, cex.names = 0.7,
          ylab = "-log10(FDR)", main = "Top exploratory methylation-gene-set results")
} else {
  plot.new()
  title("No enrichment rows were returned")
}
dev.off()

message("Enrichment complete. Ranked CpGs used: ", length(sig_cpg))
