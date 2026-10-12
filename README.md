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
| `scripts/tcga_ces.as` | autoSql schema of the hub bigBed files |

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
```

## Track contents

Each cohort track contains recurrent variants only (carried by at least 2 TCGA patients of the cohort).
The autoSql schema of the track files is in `scripts/tcga_ces.as`.

Track columns: chrom, chromStart, chromEnd, name, score, gamma, n_patients, consequence.
`score` is the transformed score s = 10*log10(gamma) - 10, rounded; bigBed requires 0-1000, so s < 0 is stored as 0
(shown as gamma "<10"). The two SKCM tracks keep their original five columns.

## Command-line access (copy and paste)

Requires `curl` and the UCSC tool `bigBedToBed`
(download from <https://hgdownload.soe.ucsc.edu/admin/exe/>, then `chmod +x bigBedToBed`).
Each command is a single line. The first line stores the hub address in a variable used by the others.

```bash
HUB=https://townsend-lab-yale.github.io/UCSC_GenomeBrowser_Hub/hg38

# 1. list the data file of every cohort track
curl -s $HUB/trackDb.txt | awk '$1=="bigDataUrl" {print $2}'

# 2. download one cohort (here BRCA) and convert it to text (BED)
curl -O $HUB/TCGA_BRCA_CES.bb
bigBedToBed TCGA_BRCA_CES.bb TCGA_BRCA_CES.bed

# 3. only chromosome 17, or only the TP53 region
bigBedToBed -chrom=chr17 TCGA_BRCA_CES.bb TCGA_BRCA_chr17.bed
bigBedToBed -chrom=chr17 -start=7668401 -end=7687550 TCGA_BRCA_CES.bb TP53.bed

# 4. one variant by name
grep -w "TP53_R175H" TCGA_BRCA_CES.bed
```

File names contain no spaces or parentheses, so they can be pasted into any shell without quoting.

## BRCA workflow notes

1. The analysis uses `cancereffectsizeR` with `ces.refset.hg38`.
2. The hub contains only recurrent-variant cancer effect size results from `ces_variant`; it does not include epistasis analyses.
3. Coordinate-style recurrent variants are parsed directly from `variant_name`.
