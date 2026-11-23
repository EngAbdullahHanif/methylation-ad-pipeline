# EPIC v1 probe mask

`EPIC.hg19.MASK_general.20180909.tsv.gz` contains the 99,360 unique probe IDs whose `MASK_general` value is `TRUE` in the Zhou lab's archived EPIC hg19 manifest dated 2018-09-09. The manifest covers the original EPIC array. The workflow matches probe IDs after detection, autosomal, and SNP filters. `MASK_general` includes mapping and other documented probe-quality flags. [Archived manifest directory](https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/) · [Mask field documentation](https://zwdzwd.github.io/InfiniumAnnotation/mask.html) · [Zhou, Laird, and Shen](https://doi.org/10.1093/nar/gkw967)

To regenerate the project file:

```bash
curl --fail --location --output EPIC.hg19.manifest.tsv.gz \
  https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz
python3 scripts/build_epic_mask.py EPIC.hg19.manifest.tsv.gz \
  data/reference/EPIC.hg19.MASK_general.20180909.tsv.gz
```

| File | SHA-256 |
| --- | --- |
| Archived manifest | `81fbb7def0759acb60385eff870f838cbbbb636ead4d4dcc2391d3ad6fb8b72e` |
| Extracted mask | `d410361b2a11ccb73aebda030575cc71d01c0f2d7f33611c78cbe85fa20cf1f6` |
