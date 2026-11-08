# PUBLIC RELEASE CANDIDATE COPY
# Source: tka_inflammation_redesign/STAGE7M_R_JOURNAL_GRADE_FIGURE_REDESIGN/write_hash_manifest.R
# Source SHA256: 4034b7c13563ebbd157773b1aa242a9bb8c90c92ed7f304c827036fa273ae470
# Purpose: Figure/source hash utility
# This copy does not overwrite the source and contains no study input data.

root <- getwd()
outdir <- file.path(root, "tka_inflammation_redesign", "STAGE7M_R_JOURNAL_GRADE_FIGURE_REDESIGN")
files <- c(
  "RA_TKA_80例.xlsx",
  "RA_TKA论文分析/临床分析/结果/04_CRP各时间点描述.csv",
  "RA_TKA论文分析/临床分析/结果/05_CRP_LMM_contrasts.csv",
  "tka_inflammation_redesign/stage4_clinical_crp_vas_20261026/STAGE4_CRP_VAS_MODEL_RESULTS.csv",
  "tka_inflammation_redesign/stage4_clinical_crp_vas_20261026/data/STAGE4_SENSITIVITY_SUMMARY.csv",
  "tka_inflammation_redesign/stage4_clinical_crp_vas_20261026/data/STAGE4_LOO_PRIMARY.csv",
  "tka_inflammation_redesign/STAGE7J_FIGURE3_FINAL_LOCK/FIGURE3_SAMPLE_SCORES.csv",
  "tka_inflammation_redesign/STAGE7J_FIGURE3_FINAL_LOCK/SUPPLEMENTARY_MODULE_RESULTS.csv",
  "tka_inflammation_redesign/STAGE7J_FIGURE3_FINAL_LOCK/ECM64_FUNCTIONAL_PANEL_MAP.csv",
  "tka_inflammation_redesign/STAGE7K_VERSION_C_FINAL_INTEGRATION/FIGURE4_SOURCE_DATA.csv",
  "tka_inflammation_redesign/STAGE7K_VERSION_C_FINAL_INTEGRATION/ECM_GENE_LEVEL_CELL_SOURCE.csv",
  "tka_inflammation_redesign/STAGE7H_R_CRP_ECM_SPATIAL_BRIDGE/STRUCTURAL_ECM_FINAL_DONOR_EFFECTS.csv",
  "tka_inflammation_redesign/STAGE7H_R_CRP_ECM_SPATIAL_BRIDGE/PROCESSING_ECM_DONOR_EFFECTS_SENSITIVITY.csv"
)
sha256 <- function(path) unname(tools::md5sum(path))
# Use system sha256 for a portable, explicit SHA256 rather than R's MD5 helper.
hash_one <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}
res <- data.frame(file=files, sha256=vapply(files, function(x) hash_one(file.path(root, x)), character(1)))
write.table(res, file.path(outdir, "FIGURE_INPUT_SHA256.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("Wrote input SHA256 manifest\n")
