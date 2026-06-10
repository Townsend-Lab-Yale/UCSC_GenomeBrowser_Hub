suppressPackageStartupMessages({
  library(cancereffectsizeR)
  library(ces.refset.hg38)
  library(data.table)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/analysis", recursive = TRUE, showWarnings = FALSE)

tcga_maf_file <- "data/raw/TCGA-BRCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "BRCA", filename = tcga_maf_file)
}

tcga_clinical <- fread(system.file("tutorial/TCGA_BRCA_clinical.txt", package = "cancereffectsizeR"))
setnames(tcga_clinical, "patient_id", "Unique_Patient_Identifier")

tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
tgs_maf_file <- system.file("tutorial/metastatic_breast_2021_hg38.maf", package = "cancereffectsizeR")
tgs_maf <- preload_maf(maf = tgs_maf_file, refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa = cesa, maf = tcga_maf, maf_name = "BRCA")
cesa <- load_sample_data(cesa = cesa, sample_data = tcga_clinical)

top_tgs_genes <- c(
 "TP53", "PIK3CA", "ESR1", "CDH1", "GATA3", "KMT2C",
 "MAP3K1", "AKT1", "ARID1A", "FOXA1", "TBX3", "PTEN"
)
tgs_coverage <- ces.refset.hg38$gr_genes[ces.refset.hg38$gr_genes$gene %in% top_tgs_genes]
tgs_maf$pM <- "M1"

cesa <- load_maf(
  cesa = cesa,
  maf = tgs_maf,
  sample_data_cols = "pM",
  maf_name = "MBC",
  coverage = "targeted",
  covered_regions = tgs_coverage,
  covered_regions_name = "top_genes",
  covered_regions_padding = 10
)

signature_exclusions <- suggest_cosmic_signature_exclusions(
  cancer_type = "BRCA",
  treatment_naive = TRUE
)

cesa <- trinuc_mutation_rates(
  cesa,
  signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)

cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)
cesa <- ces_variant(cesa = cesa, run_name = "recurrents")

group1 <- cesa$variants[c("PIK3CA E545K", "AKT1 E17K"), variant_id, on = "variant_name"]
group2 <- cesa$variants[c("PIK3CA E545K", "PIK3CA E542K"), variant_id, on = "variant_name"]

cesa <- ces_epistasis(
  cesa = cesa,
  variants = list(group1, group2),
  conf = 0.95,
  run_name = "variant_epistasis_example"
)

top_pik3ca <- cesa$variants[gene == "PIK3CA" & maf_prevalence > 1]
top_akt1 <- cesa$variants[variant_name == "AKT1 E17K"]
compound_variants <- define_compound_variants(
  cesa = cesa,
  variant_table = rbind(top_pik3ca, top_akt1),
  by = "gene",
  merge_distance = Inf
)

cesa <- ces_epistasis(
  cesa = cesa,
  variants = compound_variants,
  run_name = "AKT1_E17K_vs_PIK3CA"
)

genes <- c("AKT1", "PIK3CA", "TP53")
combined_coverage <- intersect(
  cesa$coverage_ranges$exome$`exome+`,
  cesa$coverage_ranges$targeted$top_genes
)
selected_variants <- select_variants(cesa, genes = genes, gr = combined_coverage)
cesa <- ces_gene_epistasis(
  cesa = cesa,
  genes = genes,
  variants = selected_variants,
  run_name = "gene_epistasis_example"
)

save_cesa(cesa = cesa, filename = "outputs/analysis/brca_ces_analysis.rds")

fwrite(
  cesa$variants[order(-maf_prevalence)],
  file = "outputs/analysis/brca_variants_by_prevalence.tsv",
  sep = "\t"
)

if (!is.null(cesa$dNdScv_results) && length(cesa$dNdScv_results) > 0) {
  dndscv_results <- cesa$dNdScv_results[[1]]
  fwrite(
    dndscv_results[order(qallsubs_cv)],
    file = "outputs/analysis/brca_dndscv_results.tsv",
    sep = "\t"
  )
}

writeLines(capture.output(sessionInfo()), con = "outputs/analysis/session_info.txt")

message("Saved analysis outputs to outputs/analysis/")
