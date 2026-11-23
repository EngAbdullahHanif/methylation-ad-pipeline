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
  result <- tryCatch(
    gometh(
      sig.cpg = sig_cpg,
      all.cpg = all_cpg,
      collection = collection_name,
      array.type = "EPIC",
      plot.bias = FALSE,
      prior.prob = TRUE
    ),
    error = function(error) {
      stop(collection_name, " enrichment failed: ", conditionMessage(error), call. = FALSE)
    }
  )
  description_column <- if (collection_name == "GO") "TERM" else "Description"
  if (!all(c(description_column, "P.DE", "FDR") %in% names(result))) {
    stop(collection_name, " enrichment returned an unexpected table schema.", call. = FALSE)
  }
  result$Description <- result[[description_column]]
  result$pathway_id <- rownames(result)
  result$collection <- rep(collection_name, nrow(result))
  result
}

go_results <- run_gometh("GO")
kegg_results <- run_gometh("KEGG")

combine_results <- function(first, second) {
  columns <- union(names(first), names(second))
  for (column in setdiff(columns, names(first))) {
    first[[column]] <- rep(NA, nrow(first))
  }
  for (column in setdiff(columns, names(second))) {
    second[[column]] <- rep(NA, nrow(second))
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

plot_data <- enrichment[!is.na(enrichment$Description) & !is.na(enrichment$FDR), , drop = FALSE]
plot_data <- head(plot_data, 10L)
pdf(file.path(results_dir, "enrichment_top_pathways.pdf"), width = 9, height = 7)
if (nrow(plot_data) > 0L) {
  plot_data <- plot_data[rev(seq_len(nrow(plot_data))), , drop = FALSE]
  labels <- vapply(plot_data$Description, function(description) {
    paste(strwrap(description, width = 36), collapse = "\n")
  }, character(1))
  par(mar = c(5, 16, 4, 2))
  barplot(-log10(pmax(plot_data$FDR, .Machine$double.xmin)),
          names.arg = labels, horiz = TRUE, las = 1, cex.names = 0.75,
          xlab = "-log10(FDR)", main = "Top exploratory methylation-gene-set results")
} else {
  plot.new()
  title("No enrichment rows were returned")
}
dev.off()

message("Enrichment complete. Ranked CpGs used: ", length(sig_cpg))
