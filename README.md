# DNA Methylation Array Analysis in Alzheimer's Disease

This repository contains a reproducible R/Bioconductor workflow for exploratory analysis of human brain DNA methylation array data. The project was developed to gain practical experience with methylation-array quality control, preprocessing, statistical analysis, R Markdown, and Shiny.

The analysis is a methods demonstration using public data. It is not a clinical study and does not establish a biological discovery.

## Data

The analysis uses GEO accession [GSE212682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212682), generated with the Illumina Infinium MethylationEPIC array. The series contains samples from eight brain regions and provides sample-level metadata and raw IDAT files.

The working analysis is restricted to the middle frontal gyrus. Samples with an intermediate diagnosis code were excluded. A fixed seed (`20261003`) was used to select five Alzheimer's disease samples and five controls. The small subset reflects the memory available for raw EPIC import on the development computer; it is not a biological inclusion criterion.

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

Quality control includes mean detection p-values, intensity plots, and methylation-based sex prediction. The sample-level detection p-value threshold is 0.01.

Noob is used for background and dye-bias correction. Probes are filtered using detection performance, chromosome, and CpG or single-base-extension SNP annotations. A validated EPIC cross-reactive probe list was not available in the selected package set, so that limitation is recorded and no unverified list is substituted. Missing beta values are checked after filtering; no imputation is performed when none are present.

Differential methylation is tested with limma on M-values. Beta-scale delta values are reported as effect sizes. DMRcate is used for exploratory region-level analysis. Gene-set enrichment uses `missMethyl::gometh`, which accounts for the unequal number of array probes assigned to genes.

## Results

The executed analysis contains 10 samples: 5 Alzheimer's disease samples and 5 controls. All samples passed the mean detection p-value threshold, and no sex mismatches were flagged.

Filtering retained 815,651 probes, and no missing beta values were observed after filtering. The exploratory limma analysis identified 4 probes meeting FDR < 0.05 and absolute delta-beta >= 0.05. The genomic inflation factor was 0.625, and DMRcate returned one exploratory region.

These values are descriptive outputs of a small, confounded subset. They should not be interpreted as replicated disease associations.

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
- Internet access for GEO data and package installation
- Approximately 16 GB RAM for the selected raw-data subset
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

The app reads precomputed results and does not process raw IDAT files.

## Limitations and future work

The analysis uses a small subset selected for local raw-array processing. Diagnosis is confounded with the observed batch groups, and age is recorded as `90+` for all selected samples. The data are bulk brain tissue, with no cell-type correction or replication cohort. The analysis is therefore exploratory.

A larger study should use a prespecified covariate design with overlapping batches, estimate cell composition, include an independent replication cohort, and validate important findings experimentally. Further extensions could include genotype integration and comparison with blood or additional brain regions.

## Citation

The data source and original study are documented in GEO accession GSE212682. Users should cite the original publication listed on the GEO record in addition to citing the software packages used for the analysis.
