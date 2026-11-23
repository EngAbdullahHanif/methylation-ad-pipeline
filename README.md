# DNA Methylation Array Analysis in Alzheimer's Disease

This repository contains a reproducible R/Bioconductor workflow for exploratory analysis of human brain DNA methylation array data. The project was developed to gain practical experience with methylation-array quality control, preprocessing, statistical analysis, R Markdown, and Shiny.

The analysis is a methods demonstration using public data. It is not a clinical study and does not establish a biological discovery.

## Data

The analysis uses GEO accession [GSE212682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212682), generated with the Illumina Infinium MethylationEPIC array. The series contains samples from eight brain regions and provides sample-level metadata and raw IDAT files.

The analysis uses all 29 eligible middle frontal gyrus samples, including 18 controls (diagnosis code 0) and 11 AD samples (code 2). The intermediate diagnosis group was excluded. GEO reports age as `90+` throughout this cohort and provides estimated proportions for five cell types.

Raw IDAT files are downloaded by `scripts/00_download_data.R` and are excluded from Git. The original candidate records GSE80970 and GSE59685 were inspected first; their series raw archives contained platform resources rather than sample IDAT files, so they were not used for the raw-data workflow.

## Workflow

1. Download GEO metadata and raw IDAT files.
2. Import IDATs with `minfi` and perform sample quality control.
3. Apply Noob normalization and probe filtering.
4. Explore principal components and technical variables.
5. Perform exploratory differential methylation and region-level analysis.
6. Run probe-bias-aware gene-set enrichment.
7. Render the R Markdown report and inspect precomputed results through the Shiny app.

## Methods

Quality control includes mean detection p-values, intensity plots, and methylation-based sex prediction. Samples with mean detection p-values above 0.01 or non-finite means are excluded before normalization. Sex mismatches are flagged for review.

Noob is used for background and dye-bias correction. Filtering uses detection performance, chromosome, CpG or single-base-extension SNP annotations, and the archived Zhou EPIC v1 hg19 `MASK_general` probe IDs. The mask version, source, and checksum are recorded in the [reference file documentation](data/reference/README.md). The final beta values are checked for missingness.

The primary limma model tests the AD versus control contrast on M-values with fixed effects for array batch and recorded sex, plus the estimated neuronal proportion. Array batch and bisulfite batch define the same groups in this cohort, so the models use array batch. The reported adjusted delta beta comes from fitting the same design to beta values. The sensitivity analysis fits the same model and probes on samples from batches with both diagnoses. Both designs must have full column rank and positive residual degrees of freedom.

DMRcate uses the primary design and requires at least two CpGs for a region. Exploratory enrichment applies `missMethyl::gometh` to the 100 CpGs with the smallest primary p-values, accounting for differences in the number of array probes assigned to each gene.

## Results

All 29 samples passed detection QC and no sex mismatches were flagged. The primary fit covers 12 batches and has rank 15 with 14 residual degrees of freedom. The shared-batch fit covers batches 3, 4, and 5 with 4 controls and 6 AD samples, rank 6, and 4 residual degrees of freedom.

The mask removed **70,340** of 815,581 QC-passing probes and left **745,241**. Seven of the original top 100 were excluded, 92 probes appeared in both top-100 lists, and the leading probe remained `cg12586176`. Four primary probes passed both FDR < 0.05 and absolute adjusted delta beta ≥ 0.05 after masking, compared with zero before masking. The primary genomic inflation factor after masking was **2.812**. DMRcate returned zero regions at the specified threshold.

Across the masked probe set, the primary and shared-batch adjusted M-value effects correlated at **r = 0.977**. All 100 primary top-ranked probes had the same effect direction in the shared-batch fit. The shared-batch fit had one probe meeting both cutoffs, with zero probes meeting both cutoffs in both fits.

The report includes the paired estimates for the leading probes and summarizes agreement across all probes. The workflow saves the full comparison locally as `results/differential_sensitivity.rds`.

The direct within-batch diagnosis contrast comes from the shared batches. The two analyses share those samples, so their agreement measures sensitivity to cohort restriction. Elevated genomic inflation warrants caution when interpreting the p-values.

The leave-one-batch-out analysis omits batches 3, 4, and 5 one at a time from the full cohort. It keeps the primary covariates, ranks candidates using the full fit, and recalculates genome-wide M-value FDR for each omission. The [summary](results/leave_one_batch_out_summary.tsv) reports design rank and effect stability. The [probe estimates](results/leave_one_batch_out_effects.tsv) include ordinary 95% beta-scale confidence intervals, and the [plot](results/leave_one_batch_out.pdf) shows the four leading probes.

All three refits were full rank, with 10 to 13 residual degrees of freedom. The top 100 probes kept their effect directions in every refit. Zero of the four primary hits met both discovery cutoffs after any omission. Omitting batch 5 had the largest overall impact, with a genome-wide M-value effect correlation of **r = 0.686** and a median absolute change of **1.28 percentage points** among the top 100 adjusted beta estimates.

The independent GSE59685 frontal-cortex comparison evaluated 47 of the 100 masked discovery candidates. Effects agreed in direction for 43 of 47 and correlated at **r = 0.780**. Median absolute adjusted effects were **13.2 percentage points** in discovery and **2.2** in the independent cohort. Zero candidates had a same-direction replication FDR below 0.05. [Replication analysis](docs/replication.md)

## Repository structure

```text
scripts/00_download_data.R                 Download metadata and raw IDATs
scripts/01_import_and_qc.R                 Import and quality control
scripts/02_normalisation.R                 Normalization and probe filtering
scripts/03_exploration_and_differential.R PCA, batch assessment, limma, DMRcate
scripts/06_enrichment.R                    GO and KEGG enrichment
report/analysis_report.Rmd                 Reproducible HTML report source
app/app.R                                  Shiny results viewer
data/metadata/                             Clean sample metadata
results/                                   Generated tables, figures, and objects
renv.lock                                  Package versions
run_all.sh                                 Command-line workflow
```

## Requirements

- macOS or Linux
- R 4.5.2 or a compatible R version
- Python 3.9 or newer, Bash, and curl
- Internet access for GEO data, KEGG annotations, and package installation
- Approximately 16 GB RAM recommended for the 29-sample raw-array run
- Pandoc 2.8 or newer for HTML report rendering

## Running the workflow

From the repository root:

```bash
Rscript -e 'renv::restore()'
bash run_all.sh
```

The workflow downloads the raw data, regenerates the metadata and results, runs enrichment, and renders `report/analysis_report.html`. Raw data and large local result objects are not committed.

The runner caches completed stages in `.workflow-cache`. It checks stage source and lockfile hashes, input and output file fingerprints, and upstream state before reusing the results. Run `bash run_all.sh --target normalization` to stop after processing the arrays. Run `bash run_all.sh --force-from report` to rebuild the report from valid upstream caches. Change only `app/app.R` and restart Shiny to reuse the same data directly. A cache miss or `--force-from` reruns the affected stage and its dependent stages up to the specified target.

The runner records elapsed seconds and peak resident memory for every executed stage in `.workflow-cache/benchmarks.tsv`. The 29-sample run with the cached raw archive is recorded in [`results/memory_benchmark.tsv`](results/memory_benchmark.tsv). Peak resident memory measures the stage process and is distinct from R's heap size.

For the independent analysis, run `bash run_all.sh --target replication`. The first run downloads the 1.9 GB GSE59685 beta matrix and streams the candidate rows. Cohort selection, models, and results are documented in the [replication analysis](docs/replication.md).

Run `bash run_all.sh --target leave-one-batch-out` to regenerate the batch-omission results from the cached differential run. The default HTML report also includes this analysis.

To start the interactive viewer after the workflow has run:

```r
shiny::runApp("app")
```

The app reads precomputed results. Its differential table and filters use the 100 probes with the smallest limma p-values. The probe viewer uses the locally generated beta matrix when available.

To check the saved outputs, render the report in a temporary directory, and exercise the Shiny server, run `Rscript tests/smoke_workflow.R` from the repository root. Add `--replication` to also check the independent comparison.

## Limitations and future work

This is an exploratory analysis of one brain region with limited diagnosis overlap across batches. Age is reported as `90+` throughout and is excluded from the models. Adjustment for estimated neuronal proportion leaves possible confounding by other cell types. The independent comparison tests only masked discovery candidates present on the replication platform and uses a different covariate set.

An extension to additional regions of GSE212682 would need to account for repeated samples from the same donors.

## Citation

The data source and original study are documented in GEO accession GSE212682. Users should cite the original publication listed on the GEO record in addition to citing the software packages used for the analysis.
