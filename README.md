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

Noob is used for background and dye-bias correction. Probes are filtered using detection performance, chromosome, and CpG or single-base-extension SNP annotations. A validated EPIC cross-reactive probe list was not available in the selected package set, so that limitation is recorded and no unverified list is substituted. Missing beta values are checked after filtering; no imputation is performed when none are present.

The primary limma model tests the AD versus control contrast on M-values with fixed effects for array batch and recorded sex, plus the estimated neuronal proportion. Array batch and bisulfite batch define the same groups in this cohort, so the models use array batch. The reported adjusted delta beta comes from fitting the same design to beta values. The sensitivity analysis fits the same model and probes on samples from batches with both diagnoses. Both designs must have full column rank and positive residual degrees of freedom.

DMRcate uses the primary design and requires at least two CpGs for a region. Exploratory enrichment applies `missMethyl::gometh` to the 100 CpGs with the smallest primary p-values, accounting for differences in the number of array probes assigned to each gene.

## Results

All 29 samples passed detection QC and no sex mismatches were flagged. The primary fit covers 12 batches and has rank 15 with 14 residual degrees of freedom. The shared-batch fit covers batches 3, 4, and 5 with 4 controls and 6 AD samples, rank 6, and 4 residual degrees of freedom.

Filtering retained 815,581 probes with complete beta values. The adjusted M-value effects correlated at **r = 0.977**, and all 100 probes ranked highest by the primary fit had the same effect direction in the shared-batch fit. Zero primary probes and one shared-batch probe met both FDR < 0.05 and absolute adjusted delta beta ≥ 0.05. The primary genomic inflation factor was **2.761**, and DMRcate found fewer than two individually significant CpGs for region calling.

The report includes the paired estimates for the leading probes and summarizes agreement across all probes. The workflow saves the full comparison locally as `results/differential_sensitivity.rds`.

The direct within-batch diagnosis contrast comes from the shared batches. The two analyses share those samples, so their agreement measures sensitivity to cohort restriction. Elevated genomic inflation warrants caution when interpreting the p-values.

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

To start the interactive viewer after the workflow has run:

```r
shiny::runApp("app")
```

The app reads precomputed results. Its differential table and filters use the 100 probes with the smallest limma p-values. The probe viewer uses the locally generated beta matrix when available.

To check the saved outputs, render the report in a temporary directory, and exercise the Shiny server, run `Rscript tests/smoke_workflow.R` from the repository root.

## Limitations and future work

This is an exploratory analysis of one brain region with limited diagnosis overlap across batches. Age is reported as `90+` throughout and is excluded from the models. Adjustment for estimated neuronal proportion leaves possible confounding by other cell types. The probe filtering lacks a cross-reactive probe mask, and the study has no independent replication cohort.

Further work could validate a platform-specific cross-reactive probe mask, test the same contrast in an independent brain cohort, and account for repeated donors when extending to other brain regions.

## Citation

The data source and original study are documented in GEO accession GSE212682. Users should cite the original publication listed on the GEO record in addition to citing the software packages used for the analysis.
