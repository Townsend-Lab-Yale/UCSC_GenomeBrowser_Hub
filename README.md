# UCSC GenomeBrowser Hub

This repository hosts the UCSC Genome Browser hub files for the Townsend Lab TCGA cancer effect size tracks.

It also includes a minimal reproducible **BRCA** workflow in `scripts/` for:

1. installing dependencies,
2. recreating the BRCA `CESAnalysis` object,
3. exporting the mapped BED track, and
4. converting the BED file to bigBed.

## Files added for BRCA reproducibility

| Path | Purpose |
| --- | --- |
| `scripts/01_install_dependencies.R` | Installs the R package dependencies |
| `scripts/02_run_brca_ces_analysis.R` | Recreates the BRCA `cesa` object |
| `scripts/03_build_bed_tracks.R` | Builds the BRCA BED outputs from recurrent CES results |
| `scripts/04_bed_to_bigbed.sh` | Converts BED to bigBed with UCSC `bedToBigBed` |

## Required local inputs

Place these local files under `data/raw/` before running the full workflow:

- `hg38.chrom.sizes`
- optionally `TCGA-BRCA.maf.gz`

If `TCGA-BRCA.maf.gz` is absent, the workflow downloads it automatically.

## Run order

```bash
Rscript scripts/01_install_dependencies.R
Rscript scripts/02_run_brca_ces_analysis.R
Rscript scripts/03_build_bed_tracks.R
bash scripts/04_bed_to_bigbed.sh \
  outputs/tracks/BRCA_CES_All_mapped_ucsc.bed \
  data/raw/hg38.chrom.sizes \
  hg38/Breast_cancer_stage(0-IV)_TCGA_CES.bb
```

## BRCA workflow notes

1. The analysis uses `cancereffectsizeR` with `ces.refset.hg38`.
2. The hub contains only recurrent-variant cancer effect size results from `ces_variant`; it does not include epistasis analyses.
3. Coordinate-style recurrent variants are parsed directly from `variant_name`.
4. Gene-level variants are mapped to the first matching gene interval in `ces.refset.hg38$gr_genes`.
5. The UCSC score is computed as `round(10 * log10(selection_intensity) - 10)` and clipped to the supported BED score range.
