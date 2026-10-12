table tcgaCES
"TCGA cancer effect size for recurrent somatic SNVs (seen in at least 2 patients)"
    (
    string chrom;       "Reference sequence chromosome"
    uint   chromStart;  "Start position of the variant (0-based)"
    uint   chromEnd;    "End position of the variant (exclusive)"
    string name;        "Variant identifier: GENE_AACHANGE (coding) or chr:pos_ref>alt (noncoding)"
    uint   score;       "Transformed score s = 10*log10(gamma) - 10, rounded; bigBed requires 0-1000, so s < 0 is stored as 0"
    string gamma;       "Cancer effect size (scaled selection coefficient) gamma = 10^((s+10)/10); '<10' when s is stored as 0"
    uint   n_patients;  "Number of TCGA patients in the cohort carrying the variant (GDC open-access mutation calls)"
    string consequence; "missense, silent, nonsense, stop_lost or noncoding"
    )
