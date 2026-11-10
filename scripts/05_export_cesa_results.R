# Export the recurrent-variant cancer effect sizes from saved CESA objects to one table per cohort.
# Usage: Rscript 05_export_cesa_results.R <folder with CESA_*_all_recurrent_SNV.rds> <output folder>
suppressMessages(library(data.table))
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]
files <- sort(list.files(IN, pattern = "^CESA_.*\\.rds$", full.names = TRUE))
summ <- list(); hist <- list()
for (f in files) {
  co <- sub("^CESA_(.*)_all_recurrent_SNV\\.rds$", "\\1", basename(f))
  a <- attributes(readRDS(f))
  runs <- names(a$selection_results)
  sel <- as.data.table(a$selection_results[["recurrent_variants"]])
  snv <- as.data.table(a$mutations$snv)[, .(snv_id, chr, pos)]
  key <- as.data.table(a$mutations$aac_snv_key)[, .(snv_id = snv_id[1]), by = aac_id]
  aac <- merge(key, snv, by = "snv_id")[, .(variant_id = aac_id, chr, pos)]
  both <- rbind(snv[, .(variant_id = snv_id, chr, pos)], aac)
  res <- merge(sel, both, by = "variant_id", all.x = TRUE, sort = FALSE)
  out <- res[, .(variant_name, variant_type, gene, chr, pos, selection_intensity, ci_low_95, ci_high_95,
                 included_with_variant, included_total)]
  setorder(out, chr, pos, variant_name)
  fwrite(out, file.path(OUT, paste0("TCGA_", co, "_cancer_effects.tsv")), sep = "\t")
  summ[[co]] <- data.table(cohort = co, runs = paste(runs, collapse = ","), samples = nrow(a$samples),
    recurrent_rows = nrow(out), n_ge2 = sum(out$included_with_variant >= 2), n_lt2 = sum(out$included_with_variant < 2),
    no_pos = sum(is.na(out$pos)), gamma_le1 = sum(out$selection_intensity <= 1),
    min_gamma = signif(min(out$selection_intensity), 3), max_gamma = signif(max(out$selection_intensity), 3),
    all_variants_rows = if ("all_variants" %in% runs) nrow(a$selection_results[["all_variants"]]) else NA)
  rh <- a$run_history; hist[[co]] <- if (is.null(rh)) "none" else paste(unlist(rh), collapse = " || ")
}
s <- rbindlist(summ); print(s, nrows = 50)
fwrite(s, file.path(OUT, "summary.tsv"), sep = "\t")
writeLines(paste(names(hist), unlist(hist), sep = " :: "), file.path(OUT, "run_history.txt"))
