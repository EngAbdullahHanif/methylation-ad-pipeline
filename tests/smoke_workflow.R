# Run from the repository root with Rscript tests/smoke_workflow.R.
metadata <- read.delim("data/metadata/sample_sheet_clean.tsv")
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
if (file.exists("results/differential_methylation.rds")) {
  full <- readRDS("results/differential_methylation.rds")
  expected <- head(order(full$P.Value, -abs(full$delta_beta), na.last = NA), 100L)
  expected_lambda <- median(qchisq(1 - full$P.Value, df = 1), na.rm = TRUE) /
    qchisq(0.5, df = 1)
  stopifnot(identical(top$CpG, full$CpG[expected]),
            isTRUE(all.equal(summary$genomic_inflation_factor, expected_lambda)))
}

enrichment <- read.delim("results/enrichment_results.tsv")
stopifnot(setequal(unique(enrichment$collection), c("GO", "KEGG")),
          all(nzchar(enrichment$Description, keepNA = TRUE)),
          all(nzchar(enrichment$pathway_id, keepNA = TRUE)))

app <- shiny::shinyAppDir("app")
has_beta_matrix <- file.exists("results/beta_filtered.rds")
shiny::testServer(app, {
  session$setInputs(pca_colour = "batch", fdr_cutoff = 1, delta_cutoff = 0)
  stopifnot(grepl("diagnosis", output$overview_counts, fixed = TRUE),
            !is.null(output$overview_qc), !is.null(output$pca_plot),
            nrow(filtered_diff()) == nrow(top), !is.null(output$volcano_plot),
            !is.null(output$diff_table), !is.null(output$enrichment_plot),
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
message("Workflow smoke test passed.")
