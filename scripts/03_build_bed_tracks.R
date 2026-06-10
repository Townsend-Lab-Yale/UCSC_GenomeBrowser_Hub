suppressPackageStartupMessages({
  library(cancereffectsizeR)
  library(ces.refset.hg38)
  library(data.table)
})

dir.create("outputs/tracks", recursive = TRUE, showWarnings = FALSE)

cesa <- load_cesa("outputs/analysis/brca_ces_analysis.rds")
recurrents <- as.data.table(cesa$selection$recurrents)

parse_coordinate_variants <- function(dt) {
  coord_rows <- grepl("^[0-9XY]+:", dt$variant_name)
  coord_variants <- dt[coord_rows]

  parsed <- coord_variants[, tstrsplit(variant_name, "[: ]", perl = TRUE)]
  if (ncol(parsed) < 2) {
    stop("Failed to parse coordinate-style variant names.")
  }

  setnames(parsed, c("chr", "start", paste0("extra_", seq_len(max(0, ncol(parsed) - 2)))))

  data.table(
    chr = paste0("chr", parsed$chr),
    start = as.integer(parsed$start) - 1L,
    end = as.integer(parsed$start),
    variant_name = coord_variants$variant_name,
    selection_intensity = as.numeric(coord_variants$selection_intensity)
  )
}

map_gene_level_variants <- function(dt) {
  coord_rows <- grepl("^[0-9XY]+:", dt$variant_name)
  gene_variants <- copy(dt[!coord_rows])
  gene_variants[, gene := sub(" .*", "", variant_name)]

  gene_gr <- ces.refset.hg38$gr_genes
  gene_dt <- as.data.table(gene_gr)
  gene_dt <- gene_dt[, .SD[1], by = gene]

  setkey(gene_dt, gene)
  setkey(gene_variants, gene)
  mapped <- gene_dt[gene_variants]

  data.table(
    chr = paste0("chr", mapped$seqnames),
    start = as.integer(mapped$start) - 1L,
    end = as.integer(mapped$start),
    variant_name = mapped$variant_name,
    selection_intensity = as.numeric(mapped$selection_intensity)
  )
}

make_browser_score <- function(selection_intensity) {
  score <- round((10 * log10(selection_intensity)) - 10)
  score[!is.finite(score)] <- 0
  pmin(1000L, pmax(0L, as.integer(score)))
}

coord_bed <- parse_coordinate_variants(recurrents)
gene_bed <- map_gene_level_variants(recurrents)

track_values <- rbindlist(list(coord_bed, gene_bed), use.names = TRUE, fill = TRUE)
track_values[, variant_name := gsub(" ", "_", variant_name)]
track_values[, log_selection_intensity := log(selection_intensity)]
track_values[, browser_score := make_browser_score(selection_intensity)]
setorder(track_values, chr, start, end, variant_name)

fwrite(
  track_values,
  file = "outputs/tracks/BRCA_CES_All_mapped_values.tsv",
  sep = "\t"
)

ucsc_track_line <- paste0(
  'track name="BRCA_CES_All" ',
  'description="All BRCA mapped CES" ',
  'visibility=2 useScore=1 color=0,0,255'
)

writeLines(ucsc_track_line, con = "outputs/tracks/BRCA_CES_All_mapped_ucsc.bed")
write.table(
  track_values[, .(chr, start, end, variant_name, browser_score)],
  file = "outputs/tracks/BRCA_CES_All_mapped_ucsc.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)

write.table(
  track_values[, .(chr, start, end, variant_name, browser_score)],
  file = "outputs/tracks/BRCA_CES_All_mapped_headerless.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Saved mapped BED outputs to outputs/tracks/")
