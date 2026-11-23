library(shiny)
library(DT)
library(plotly)

# Load all app data once at startup. Raw IDATs are never processed here.
results_dir <- file.path("..", "results")
metadata_dir <- file.path("..", "data", "metadata")

read_required <- function(path, columns) {
  if (!file.exists(path)) {
    stop("Missing app input: ", path, ". Run the pipeline first.")
  }
  data <- read.delim(path, check.names = FALSE)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop("Missing columns in ", path, ": ", paste(missing_columns, collapse = ", "),
         ". Run the pipeline to regenerate this input.")
  }
  data
}

sample_sheet <- read_required(file.path(metadata_dir, "sample_sheet_clean.tsv"),
                              c("sample_id", "diagnosis"))
pca_scores <- read_required(file.path(results_dir, "pca_scores.tsv"),
                           c("sample_id", "PC1", "PC2", "diagnosis", "sex", "batch", "bisulfite_batch"))
differential <- read_required(file.path(results_dir, "differential_top.tsv"),
                             c("CpG", "delta_beta", "P.Value", "adj.P.Val"))
model_design <- read_required(file.path(results_dir, "model_design.tsv"),
                              c("analysis", "n_samples", "n_control", "n_ad", "n_batches",
                                "design_rank", "residual_df"))
sensitivity <- read_required(file.path(results_dir, "differential_sensitivity_summary.tsv"),
                              c("m_effect_pearson_r", "top_m_direction_agreement",
                                "shared_passing_both_cutoffs"))
enrichment <- read_required(file.path(results_dir, "enrichment_results.tsv"),
                           c("Description", "FDR", "collection", "pathway_id"))
qc_summary <- read_required(file.path(results_dir, "qc_summary.tsv"), "qc_status")
sex_check <- read_required(file.path(results_dir, "sex_check.tsv"), "sex_mismatch")

beta_file <- file.path(results_dir, "beta_filtered.rds")
beta_values <- if (file.exists(beta_file)) readRDS(beta_file) else NULL

pca_scores$diagnosis <- as.character(pca_scores$diagnosis)
enrichment <- enrichment[!is.na(enrichment$Description) & !is.na(enrichment$FDR), , drop = FALSE]

make_pca_plot <- function(data, colour_by) {
  groups <- factor(data[[colour_by]])
  plot_ly(data, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers",
          text = ~sample_id, hoverinfo = "text+x+y", color = groups,
          colors = grDevices::hcl.colors(nlevels(groups), "Dark 3"), marker = list(size = 10)) |>
    layout(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"),
           legend = list(title = list(text = colour_by)))
}

make_volcano_plot <- function(data) {
  data$neg_log10_p <- -log10(pmax(data$P.Value, .Machine$double.xmin))
  plot_ly(data, x = ~delta_beta, y = ~neg_log10_p, type = "scatter", mode = "markers",
          text = ~paste("CpG:", CpG, "<br>FDR:", signif(adj.P.Val, 3)),
          hoverinfo = "text", marker = list(size = 7)) |>
    layout(xaxis = list(title = "Adjusted delta beta (AD - Control)"),
           yaxis = list(title = "-log10(p-value)"))
}

make_probe_plot <- function(query) {
  if (is.null(beta_values)) {
    return(NULL)
  }
  probe_id <- rownames(beta_values)[tolower(rownames(beta_values)) == tolower(query)][1]
  if (is.na(probe_id)) {
    return(NULL)
  }
  values <- data.frame(
    sample_id = sub("_.*$", "", colnames(beta_values)),
    beta = as.numeric(beta_values[probe_id, ]),
    stringsAsFactors = FALSE
  )
  values <- merge(values, sample_sheet[, c("sample_id", "diagnosis")], by = "sample_id")
  plot_ly(values, x = ~diagnosis, y = ~beta, color = ~diagnosis,
          type = "box", jitter = 0.25, pointpos = 0, boxpoints = "all",
          text = ~sample_id, hoverinfo = "text+y") |>
    layout(xaxis = list(title = "Diagnosis"), yaxis = list(title = "Observed beta value"))
}

ui <- navbarPage(
  "AD Methylation Explorer",
  tabPanel("Overview",
           h3("GSE212682 exploratory analysis"),
           p("This app displays precomputed results. It does not process raw IDAT files."),
           tableOutput("overview_counts"),
           tableOutput("overview_qc"),
           h4("Model design"), tableOutput("overview_design"),
           p("The models adjust for array batch, recorded sex, and estimated neuronal proportion. The sensitivity model uses batches with both diagnoses.")),
  tabPanel("Sample explorer",
           selectInput("pca_colour", "Color PCA by", choices = c("diagnosis", "sex", "batch", "bisulfite_batch")),
           plotlyOutput("pca_plot", height = "600px")),
  tabPanel("Differential results",
           sidebarLayout(
             sidebarPanel(
               p("Filters apply to the top 100 probes ranked by limma p-value."),
               sliderInput("fdr_cutoff", "Maximum FDR", min = 0, max = 1, value = 0.2, step = 0.01),
               sliderInput("delta_cutoff", "Minimum absolute adjusted delta beta", min = 0, max = 1, value = 0, step = 0.01),
               downloadButton("download_diff", "Download filtered results")
             ),
             mainPanel(plotlyOutput("volcano_plot", height = "500px"), DTOutput("diff_table"),
                       h4("Shared-batch sensitivity"), tableOutput("sensitivity_summary"))
           )),
  tabPanel("Probe viewer",
           textInput("probe_query", "CpG ID", value = "", placeholder = "Example: cg00000029"),
           uiOutput("probe_note"),
           plotlyOutput("probe_plot", height = "500px")),
  tabPanel("Enrichment",
           plotlyOutput("enrichment_plot", height = "500px"),
           DTOutput("enrichment_table"))
)

server <- function(input, output, session) {
  filtered_diff <- reactive({
    differential[differential$adj.P.Val <= input$fdr_cutoff &
                   abs(differential$delta_beta) >= input$delta_cutoff, , drop = FALSE]
  })

  output$overview_counts <- renderTable({
    as.data.frame(table(sample_sheet$diagnosis), responseName = "samples") |>
      setNames(c("diagnosis", "samples"))
  })

  output$overview_qc <- renderTable({
    data.frame(
      samples = nrow(qc_summary),
      failed_detection_qc = sum(qc_summary$qc_status == "fail"),
      sex_mismatches = sum(sex_check$sex_mismatch, na.rm = TRUE)
    )
  })

  output$overview_design <- renderTable({
    data.frame(Analysis = c("Primary", "Shared batches")[match(model_design$analysis,
                                                               c("primary", "shared_batches"))],
               Samples = model_design$n_samples, Controls = model_design$n_control,
               AD = model_design$n_ad, Batches = model_design$n_batches,
               Rank = model_design$design_rank, `Residual df` = model_design$residual_df,
               check.names = FALSE)
  })

  output$sensitivity_summary <- renderTable({
    data.frame(`M-value effect correlation` = sensitivity$m_effect_pearson_r,
               `Top-100 direction agreement (%)` = 100 * sensitivity$top_m_direction_agreement,
               `Shared-batch probes passing both cutoffs` = sensitivity$shared_passing_both_cutoffs,
               check.names = FALSE)
  })

  output$pca_plot <- renderPlotly({
    make_pca_plot(pca_scores, input$pca_colour)
  })

  output$volcano_plot <- renderPlotly({
    make_volcano_plot(filtered_diff())
  })

  output$diff_table <- renderDT({
    datatable(filtered_diff(), filter = "top", options = list(pageLength = 10, scrollX = TRUE))
  })

  output$download_diff <- downloadHandler(
    filename = function() paste0("differential_results_fdr_", input$fdr_cutoff, ".tsv"),
    content = function(file) write.table(filtered_diff(), file, sep = "\t", quote = FALSE, row.names = FALSE)
  )

  output$probe_note <- renderUI({
    if (is.null(beta_values)) {
      return(p("Run the workflow locally to enable the probe viewer."))
    }
    p("Enter an exact CpG ID to plot the observed beta values by diagnosis.")
  })

  output$probe_plot <- renderPlotly({
    req(input$probe_query)
    plot <- make_probe_plot(input$probe_query)
    validate(need(!is.null(plot), "No exact CpG ID was found in the filtered beta matrix."))
    plot
  })

  output$enrichment_plot <- renderPlotly({
    validate(need(nrow(enrichment) > 0L, "No enrichment results are available."))
    top <- head(enrichment[order(enrichment$FDR), ], 15)
    plot_ly(top, x = ~FDR, y = ~reorder(Description, FDR), color = ~collection,
            type = "bar", orientation = "h", text = ~pathway_id, hoverinfo = "text+x") |>
      layout(xaxis = list(title = "FDR"), yaxis = list(title = "Pathway"), barmode = "group")
  })

  output$enrichment_table <- renderDT({
    datatable(head(enrichment[order(enrichment$FDR), ], 100), filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })
}

shinyApp(ui, server)
