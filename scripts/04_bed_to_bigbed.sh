#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: bash scripts/04_bed_to_bigbed.sh <input.bed> <hg38.chrom.sizes> <output.bb>" >&2
  exit 1
fi

input_bed="$1"
chrom_sizes="$2"
output_bb="$3"

if [[ ! -f "$input_bed" ]]; then
  echo "Missing BED file: $input_bed" >&2
  exit 1
fi

if [[ ! -f "$chrom_sizes" ]]; then
  echo "Missing chromosome sizes file: $chrom_sizes" >&2
  exit 1
fi

if ! command -v bedToBigBed >/dev/null 2>&1; then
  echo "bedToBigBed is not on PATH. Install UCSC kentUtils or download the binary from UCSC." >&2
  exit 1
fi

tmp_bed="$(mktemp)"
tmp_sorted="$(mktemp)"
trap 'rm -f "$tmp_bed" "$tmp_sorted"' EXIT

grep -vE '^(track|browser)' "$input_bed" >"$tmp_bed"
sort -k1,1 -k2,2n "$tmp_bed" >"$tmp_sorted"

bedToBigBed "$tmp_sorted" "$chrom_sizes" "$output_bb"
echo "Wrote $output_bb"
