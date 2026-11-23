# Independent frontal cortex comparison

The discovery analysis uses middle frontal gyrus samples from the 90+ Study and the EPIC v1 array. The independent comparison uses frontal cortex samples from the MRC London Neurodegenerative Disease Brain Bank and the Illumina 450K array. The GEO AD and control labels (`AD` and `C`) define the London contrast. The comparison includes 84 eligible arrays (60 AD, 24 controls), one per donor. [Discovery study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9710143/) · [London study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4410018/) · [GSE59685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685)

## Method

The analysis selects the 100 CpGs with the smallest primary discovery p-values after EPIC v1 masking, without reference to the London outcomes. It aligns GEO accessions and array barcodes to the published beta matrix, then keeps measured CpGs with complete, variable beta values. The London model fits `beta ~ diagnosis + age + sex + chip` and reports the adjusted AD-minus-control beta coefficient, its ordinary 95% confidence interval, and Benjamini–Hochberg FDR within the tested candidate set. The discovery coefficient is on the same beta scale and comes from the discovery model adjusted for batch, sex, and estimated neuronal proportion.

## Outputs

The London model had full rank 15 and 69 residual degrees of freedom. The published beta matrix covered **47 of the 100** candidates, and all 47 passed the measurement checks. Adjusted effects agreed in direction for **43 of 47 (91.5%)** and had a Pearson correlation of **0.780**. The median absolute adjusted delta beta was **0.1323** in discovery and **0.0215** in London, a difference in magnitude from 13.2 to 2.2 percentage points. Zero candidates passed replication FDR < 0.05 in the discovery direction.

The matrix covered one of the four discovery probes that met both discovery cutoffs. At `cg17462978`, the discovery effect was −0.1871 and the London estimate was −0.0116 (95% CI −0.0467 to 0.0236, FDR 0.659).

- [Summary](../results/replication_GSE59685_summary.tsv)
- [Probe-level estimates and availability](../results/replication_GSE59685_effects.tsv)
- [Adjusted effect comparison](../results/replication_GSE59685.pdf)

The analysis compares candidate effects across differently defined frontal regions, age distributions, case definitions, platforms, and covariates. Platform coverage limits which EPIC candidates can be compared with the published 450K matrix. Selecting candidates by their discovery p-values can inflate the discovery effect magnitudes.

## Reproduce

```bash
bash run_all.sh --target replication
Rscript tests/smoke_workflow.R --replication
```

On the first run, the replication stage downloads the GEO metadata and the 1.9 GB processed beta matrix to `data/raw/GSE59685/`. It streams the requested probe rows and caches the completed stage.

| Input | SHA-256 used in this run |
| --- | --- |
| [GSE59685 series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59685/matrix/GSE59685_series_matrix.txt.gz) | `7bb681ebc7398f8e6e371de9877eed96e9ad58728cf6c0ed56ff6fc14e173529` |
| [GSE59685 processed betas](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59685/suppl/GSE59685_betas.csv.gz) | `0704d771d67687754660a4814e1d24e5f19dc8e261faed77839ff5280d4cce20` |
