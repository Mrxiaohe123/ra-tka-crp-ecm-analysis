# PUBLIC RELEASE CANDIDATE COPY
# Source: VERSION_C_PRE_REVISION_TARGETED_CHECKS/CODE/targeted_checks.R
# Source SHA256: 1bb7e4cf9c9882b7a148f76d0ecb2cfdfd30bcc55dd3c2a3a1b3c5ab77b1ad17
# Purpose: Fixed GSE198520 component audit source; public inputs and frozen mapping required
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(digest)
})

root <- normalizePath(Sys.getenv('VERSION_C_PROJECT_ROOT', '.'), mustWork = TRUE)
out <- file.path(root, "VERSION_C_PRE_REVISION_TARGETED_CHECKS")
res <- file.path(out, "RESULTS")
dir.create(res, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) digest::digest(file = path, algo = "sha256")
fisher_ci <- function(r, n) {
  if (!is.finite(r) || n < 6L || abs(r) >= 1) return(c(NA_real_, NA_real_))
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)
  tanh(z + c(-1, 1) * qnorm(0.975) * se)
}
spearman_result <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 4L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.table(N = n, rho = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P_value = NA_real_))
  }
  r <- suppressWarnings(cor(x, y, method = "spearman"))
  p <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  ci <- fisher_ci(r, n)
  data.table(N = n, rho = r, CI_low = ci[1], CI_high = ci[2], P_value = p)
}

# 1. Fixed ECM components in GSE198520
gmt_file <- file.path(root, "tka_inflammation_redesign/rescue_discovery_20261026/data/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
g2_file <- file.path(root, "RA_TKA论文分析/gate3/raw/GSE198520_Raw_gene_count_matrix.txt.gz")
g2_manifest_file <- file.path(root, "RA_TKA论文分析/SSIP分析/系统桥接/临床桥接_GSE198520/结果/01_GSE198520_sample_manifest.csv")

gmt <- lapply(readLines(gmt_file, warn = FALSE), function(line) {
  p <- strsplit(line, "\t", fixed = TRUE)[[1]]
  list(name = p[1], genes = unique(toupper(p[-c(1, 2)])))
})
gset <- function(name) gmt[[which(vapply(gmt, function(x) x$name, character(1)) == name)[1]]]$genes
ecm64 <- gset("REACTOME_COLLAGEN_DEGRADATION")
stopifnot(length(ecm64) == 64L)
structural <- ecm64[grepl("^COL", ecm64)]
processing <- setdiff(ecm64, structural)
stopifnot(length(structural) == 38L, length(processing) == 26L)

g2 <- fread(g2_file, check.names = FALSE, showProgress = FALSE)
g2_gene <- toupper(trimws(as.character(g2[[1]])))
g2_id <- names(g2)[-1]
g2_counts <- as.matrix(g2[, -1, with = FALSE])
storage.mode(g2_counts) <- "numeric"
gm <- fread(g2_manifest_file, na.strings = c("", "NA"), showProgress = FALSE)
gm <- gm[match(sub("_bx$", "", g2_id), Sample_ID)]
stopifnot(all(gm$Sample_ID == sub("_bx$", "", g2_id)), ncol(g2_counts) == 92L)

g2_dge <- calcNormFactors(DGEList(counts = g2_counts), method = "TMM")
g2_lcpm <- log2(cpm(g2_dge, normalized.lib.sizes = TRUE) + 1)
g2_rs <- rowsum(g2_lcpm, group = g2_gene, reorder = FALSE)
g2_pre <- gm$Timepoint == "Pre"
g2_mu <- rowMeans(g2_rs[, g2_pre, drop = FALSE])
g2_sd <- apply(g2_rs[, g2_pre, drop = FALSE], 1, sd)
g2_valid <- is.finite(g2_sd) & g2_sd > 0
g2_z <- (g2_rs[g2_valid, , drop = FALSE] - g2_mu[g2_valid]) / g2_sd[g2_valid]

component_name <- ifelse(ecm64 %in% structural, "Structural_COL", "NonCOL_Processing")
component_map <- data.table(
  gene_symbol = ecm64,
  component = component_name,
  raw_row_present = ecm64 %in% rownames(g2_rs),
  baseline_sd = as.numeric(g2_sd[match(ecm64, rownames(g2_rs))]),
  included_in_fixed_62_score = ecm64 %in% rownames(g2_z),
  exclusion_reason = fifelse(ecm64 %in% rownames(g2_z), "None", "Zero baseline SD after fixed transform")
)
fwrite(component_map, file.path(res, "GSE198520_ECM_COMPONENT_MAPPING.csv"))

score_component <- function(genes) {
  p <- intersect(genes, rownames(g2_z))
  if (!length(p)) return(rep(NA_real_, ncol(g2_z)))
  s <- colMeans(g2_z[p, , drop = FALSE])
  names(s) <- colnames(g2_z)
  s
}
scores <- list(Structural = score_component(structural), NonCOL_Processing = score_component(processing))
score_dt <- rbindlist(lapply(names(scores), function(comp) {
  s <- scores[[comp]]
  data.table(Sample_ID = sub("_bx$", "", names(s)), Patient_ID = gm$Patient_ID,
             Timepoint = gm$Timepoint, Score = as.numeric(s), CRP = as.numeric(gm$CRP), component = comp)
}))

assoc_rows <- list()
for (comp in names(scores)) {
  z <- score_dt[component == comp]
  baseline <- z[Timepoint == "Pre"]
  b <- spearman_result(baseline$Score, baseline$CRP)
  ptab <- dcast(z, Patient_ID ~ Timepoint, value.var = "Score")
  ctab <- dcast(z, Patient_ID ~ Timepoint, value.var = "CRP")
  setnames(ptab, c("Pre", "Post"), c("Score_Pre", "Score_Post"))
  setnames(ctab, c("Pre", "Post"), c("CRP_Pre", "CRP_Post"))
  w <- merge(ptab, ctab, by = "Patient_ID")
  w[, delta_score := Score_Post - Score_Pre]
  w[, delta_crp := CRP_Post - CRP_Pre]
  ok <- complete.cases(w[, .(delta_score, delta_crp)])
  dlt <- spearman_result(w$delta_score[ok], w$delta_crp[ok])
  use <- component_map$component == ifelse(comp == "Structural", "Structural_COL", "NonCOL_Processing") & component_map$included_in_fixed_62_score
  genes_used <- paste(component_map$gene_symbol[use], collapse = ";")
  for (item in list(list(name = "Baseline_score_vs_baseline_CRP", stat = b),
                    list(name = "Within_patient_delta_score_vs_delta_CRP_12wk", stat = dlt))) {
    assoc_rows[[length(assoc_rows) + 1L]] <- data.table(
      dataset = "GSE198520", analysis = item$name, component = comp,
      score_gene_count = sum(use), score_genes = genes_used,
      N = item$stat$N, rho = item$stat$rho, CI_low = item$stat$CI_low,
      CI_high = item$stat$CI_high, P_value = item$stat$P_value,
      CI_method = "Fisher z approximate 95% CI",
      score_method = "TMM -> log2(CPM+1) -> baseline gene-wise z-score -> unweighted within-component mean",
      standardization_reference = "46 baseline samples; same baseline means/SDs applied to post samples",
      zero_handling = "CRP values retained as recorded; no pseudocount",
      analysis_status = "Post-hoc fixed-component analysis; not prespecified; no gene replacement"
    )
  }
}
component_results <- rbindlist(assoc_rows, fill = TRUE)
fwrite(component_results, file.path(res, "GSE198520_ECM_COMPONENT_ASSOCIATIONS.csv"))
historical_module_file <- file.path(root, "tka_inflammation_redesign/STAGE7H_R_CRP_ECM_SPATIAL_BRIDGE/GSE198520_FIXED_MODULE_REPLICATION.csv")
if (file.exists(historical_module_file)) {
  old <- fread(historical_module_file, showProgress = FALSE)
  current_cmp <- component_results[analysis == "Baseline_score_vs_baseline_CRP",
                                   .(component, current_N = N, current_rho = rho, current_P = P_value)]
  old_cmp <- old[signature %in% c("Structural_all", "Processing_all"),
                 .(signature, historical_N = baseline_N, historical_rho = baseline_rho, historical_P = baseline_P)]
  current_cmp[, signature := fifelse(component == "Structural", "Structural_all", "Processing_all")]
  cmp <- merge(current_cmp, old_cmp, by = "signature", all = TRUE)
  cmp[, c("delta_N", "delta_rho", "delta_P", "status") :=
        list(current_N - historical_N, current_rho - historical_rho, current_P - historical_P,
             ifelse(abs(current_rho - historical_rho) < 1e-12 &
                      abs(current_P - historical_P) < 1e-12, "NO_CHANGE", "DIFFERENCE_REQUIRES_REVIEW"))]
  fwrite(cmp, file.path(res, "GSE198520_COMPONENT_RECONCILIATION.csv"))
}

# 2. Clinical laterality and index-time aggregate audit
x <- read_excel(file.path(root, "RA_TKA_80例.xlsx"), sheet = "门诊", col_names = FALSE, .name_repair = "minimal")
d <- as.data.frame(x[-c(1, 2), , drop = FALSE])
names(d) <- as.character(x[2, ])
for (v in setdiff(names(d), c("Patient_ID", "Exclude_Reason"))) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$Surg_Date <- as.Date(d$Surg_Date, origin = "1899-12-30")
d$Pre_Date <- as.Date(d$Pre_Date, origin = "1899-12-30")
d$days_pre_to_surgery <- as.numeric(d$Surg_Date - d$Pre_Date)

requested <- data.table(raw_code = c(0, 1, 2),
                        requested_label = c("right unilateral", "left unilateral", "staged bilateral"),
                        requested_N = c(67L, 8L, 5L))
actual <- data.table(raw_code = sort(unique(d$Bilateral)))
actual[, c("source_workbook_N", "date_min_days", "date_median_days", "date_max_days",
           "unique_surgery_dates", "unique_pre_dates") :=
         list(as.integer(tabulate(match(d$Bilateral, raw_code))),
              NA_real_, NA_real_, NA_real_, NA_integer_, NA_integer_)]
for (i in seq_len(nrow(actual))) {
  g <- d$Bilateral == actual$raw_code[i]
  actual[i, c("date_min_days", "date_median_days", "date_max_days",
              "unique_surgery_dates", "unique_pre_dates") :=
            list(min(d$days_pre_to_surgery[g], na.rm = TRUE),
                 median(d$days_pre_to_surgery[g], na.rm = TRUE),
                 max(d$days_pre_to_surgery[g], na.rm = TRUE),
                 length(unique(d$Surg_Date[g])),
                 length(unique(d$Pre_Date[g])))]
}
lat <- merge(requested, actual, by = "raw_code", all = TRUE)
lat[, count_match_to_requested := source_workbook_N == requested_N]
lat[, coding_interpretation_status := fifelse(
  count_match_to_requested,
  "COUNT_MATCH_ONLY; side semantics still require source codebook",
  "COUNT_MISMATCH; do not map 0/1 to right/left without author record"
)]
fwrite(lat, file.path(res, "CLINICAL_LATERALITY_AND_INDEX_TIME_AGGREGATE.csv"))

writeLines(c(
  "# Clinical laterality and index-time source-fact audit",
  "",
  paste0("Workbook: RA_TKA_80例.xlsx; sheet=门诊; unique records=", nrow(d), "."),
  "",
  "VERIFIED from the workbook:",
  "- Bilateral has raw codes 0, 1 and 2 with counts 49, 26 and 5, respectively.",
  "- The workbook contains one Surg_Date field and one Pre_Date field per record.",
  "- For raw code 2, Pre_Date precedes Surg_Date by 3–4 days (median 3 days).",
  "- There is no second-operation date field and no machine-readable side dictionary in this workbook.",
  "",
  "CORRECTION_REQUIRED / AUTHOR CHECK:",
  "- The requested 67 right unilateral / 8 left unilateral / 5 staged bilateral distribution is not reproduced by the workbook's raw Bilateral counts.",
  "- The raw 0/1 values must not be assigned to right/left by assumption.",
  "- Whether Surg_Date represents the index or last TKA for the five raw-code-2 records cannot be verified from Excel alone; operative records are required.",
  "- No patient-level identifiers are exported in this audit."
), file.path(res, "CLINICAL_LATERALITY_SOURCE_FACTS.md"))

# 3. M12 intermediate-follow-up auxiliary missingness check
d$m12_available <- as.integer(!is.na(d$M12_CRP))
d$m3_available <- as.integer(!is.na(d$M3_CRP))
d$m6_available <- as.integer(!is.na(d$M6_CRP))
d$log2_Pre_CRP <- log2(d$Pre_CRP)

sum_rows <- list()
for (grp in c(1L, 0L)) {
  for (v in c("M3_CRP", "M6_CRP")) {
    g <- d$m12_available == grp
    xval <- d[[v]][g]
    sum_rows[[length(sum_rows) + 1L]] <- data.table(
      M12_group = ifelse(grp == 1L, "M12_CRP_observed", "M12_CRP_missing"),
      variable = v, group_N = sum(g), observed_N = sum(is.finite(xval)),
      observed_percent = 100 * mean(is.finite(xval)), median = median(xval, na.rm = TRUE),
      Q1 = as.numeric(quantile(xval, 0.25, na.rm = TRUE)),
      Q3 = as.numeric(quantile(xval, 0.75, na.rm = TRUE)), mean = mean(xval, na.rm = TRUE),
      min = min(xval, na.rm = TRUE), max = max(xval, na.rm = TRUE)
    )
  }
}
fwrite(rbindlist(sum_rows), file.path(res, "M12_INTERMEDIATE_FOLLOWUP_SUMMARY.csv"))

avail_tab <- rbindlist(lapply(c("m3_available", "m6_available"), function(v) {
  tab <- as.data.table(table(
    M12_group = ifelse(d$m12_available == 1, "M12_CRP_observed", "M12_CRP_missing"),
    available = d[[v]]
  ))
  setnames(tab, "N", "count")
  tab[, variable := v]
  tab
}))
fwrite(avail_tab, file.path(res, "M12_INTERMEDIATE_AVAILABILITY_TABLE.csv"))

fit_base <- glm(m12_available ~ log2_Pre_CRP, family = binomial(), data = d)
fit_aux <- glm(m12_available ~ log2_Pre_CRP + m3_available + m6_available, family = binomial(), data = d)
extract_glm <- function(fit, model_name) {
  sm <- summary(fit)$coefficients
  ci <- suppressWarnings(confint.default(fit))
  interp <- c("intercept", "OR per CRP doubling", "OR per M3 availability indicator",
              "OR per M6 availability indicator")[seq_len(nrow(sm))]
  data.table(model = model_name, term = rownames(sm), beta_log_odds = sm[, 1], SE = sm[, 2],
             OR = exp(sm[, 1]), CI_low_OR = exp(ci[, 1]), CI_high_OR = exp(ci[, 2]),
             P_value = sm[, 4], N = nobs(fit), AIC = AIC(fit),
             residual_deviance = deviance(fit), interpretation = interp)
}
aux_models <- rbind(extract_glm(fit_base, "baseline_availability_model"),
                    extract_glm(fit_aux, "auxiliary_intermediate_availability_model"))
lrt <- anova(fit_base, fit_aux, test = "LRT")
aux_models[, added_terms_LRT_P := as.numeric(lrt[2, "Pr(>Chi)"])]
fwrite(aux_models, file.path(res, "M12_AUXILIARY_MAR_MODEL.csv"))

writeLines(c(
  "# M12 auxiliary MAR / intermediate-follow-up check",
  "",
  "## Descriptive facts",
  paste0("M12 CRP observed: ", sum(d$m12_available == 1), "/", nrow(d),
         "; missing: ", sum(d$m12_available == 0), "."),
  "Among M12-observed records: M3 CRP observed in 45/46 and M6 CRP in 41/46.",
  "Among M12-missing records: M3 CRP observed in 30/34 and M6 CRP in 16/34.",
  "The M12-missing group had higher observed M3 and M6 CRP distributions, but these are descriptive and post-baseline observations.",
  "",
  "## Auxiliary model",
  "The single predeclared auxiliary availability model was: M12 availability ~ log2(Pre CRP) + M3 availability indicator + M6 availability indicator.",
  "This model uses intermediate follow-up availability only to diagnose selection; it does not impute M12 outcomes and does not replace the frozen MAR/MNAR estimands.",
  "Because M3/M6 are post-baseline and are themselves incompletely observed, the auxiliary model is not treated as a fully identified outcome-MAR model.",
  "No new delta grid or outcome-favourable missingness scenario was introduced.",
  "",
  "## Limitation",
  "The workbook does not document a complete survival registry or definitive reasons for 29 M12-missing records. Intermediate follow-up information therefore supports an auxiliary missingness diagnosis, not a causal or fully identified missing-data correction."
), file.path(res, "M12_AUXILIARY_MAR_METHODS.md"))

# 4. Statistics, software and reference provenance audit
pathway_file <- file.path(root, "tka_inflammation_redesign/stage6_version_b_integration_20261026/STAGE6_SUPPLEMENTARY_TABLES/GSE225731_PATHWAY_DISCOVERY_ALL.csv")
pathways <- fread(pathway_file, showProgress = FALSE)
keep_pathways <- pathways[toupper(pathway) %in% toupper(c(
  "REACTOME_COLLAGEN_DEGRADATION",
  "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN",
  "REACTOME_GPCR_LIGAND_BINDING"
))]
keep_pathways[, audit_source := "GSE225731_PATHWAY_DISCOVERY_ALL.csv; original 1,389 tested Reactome pathways; BH across tested family"]
keep_pathways[, status := "VERIFIED"]
fwrite(keep_pathways, file.path(res, "REACTOME_THREE_BH_HITS_VERIFIED.csv"))

prov <- rbind(
  data.table(
    domain = "M6_CRP_VAS_expanded_model",
    item = c("formula", "N", "coefficient", "HC3_CI", "P_value", "software"),
    value = c(
      "M6_VAS ~ Pre_VAS + log2(Pre_CRP) + log2(M6_CRP) + Age + factor(Sex) + RA_Duration_y + BMI + factor(Base_IL6i)",
      "50", "0.7684063385", "0.0625328257 to 1.4742798514", "0.0389034077",
      "R 4.5.2; readxl 1.4.5; sandwich 3.1-1; lmtest 0.9-40; data.table 1.18.0; MASS 7.3-65"
    ),
    status = "VERIFIED",
    source = "tka_inflammation_redesign/code/15_stage4_crp_vas_20261026.R; STAGE4_CRP_VAS_MODEL_RESULTS.csv"
  ),
  data.table(
    domain = "Reactome_pathway_screen",
    item = c("tested_family_N", "pathways", "BH_scope", "score_method"),
    value = c(
      "1,389",
      paste(keep_pathways$pathway, collapse = " | "),
      "BH q values were adjusted over the tested 1,389 pathways",
      "TMM log2(CPM+1) -> gene-wise z -> within-set mean"
    ),
    status = "VERIFIED",
    source = "GSE225731_PATHWAY_DISCOVERY_ALL.csv; rescue_discovery_20261026/01_GSE225731_pathway_discovery_sessionInfo.txt"
  ),
  data.table(
    domain = "Reactome_snapshot",
    item = c("file", "version", "SHA256", "release_manifest"),
    value = c(
      basename(gmt_file), "v2023.2.Hs local symbols GMT", sha256(gmt_file),
      "external release manifest not recovered"
    ),
    status = c("VERIFIED", "VERIFIED", "VERIFIED", "UNRESOLVED"),
    source = "tka_inflammation_redesign/rescue_discovery_20261026/data/c2.cp.reactome.v2023.2.Hs.symbols.gmt"
  ),
  data.table(
    domain = "GSE198520_baseline_definition",
    item = c("patients", "baseline_time", "sampling", "treatment_context", "source"),
    value = c(
      "46 RA patients", "before anti-TNF initiation",
      "minimally invasive ultrasound-guided synovial biopsy of the most inflamed joint",
      "eligible to start anti-TNF; 12-week paired post-treatment samples",
      "GEO GSE198520"
    ),
    status = "VERIFIED_FROM_GEO_RECORD",
    source = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198520"
  ),
  data.table(
    domain = "Xenium_reference",
    item = c("reference_list_position", "first_in_text_mention", "paper"),
    value = c(
      "4 in current manuscript reference list",
      "Introduction line 38 after Richard 2025",
      "Bhamidipati et al., Nat Immunol 2026; DOI 10.1038/s41590-025-02386-2"
    ),
    status = c("VERIFIED_LOCAL", "VERIFIED_LOCAL", "VERIFIED_EXTERNAL"),
    source = "VERSION_C_FULL_MANUSCRIPT_EN.md; PubMed/Publisher record"
  )
)
fwrite(prov, file.path(res, "STATISTICS_SOFTWARE_REFERENCE_AUDIT.csv"))

writeLines(c(
  "# Reference order and provenance audit",
  "",
  "## GSE198520 baseline wording",
  "The official GEO series record describes 46 RA patients with synovial RNA-seq before and 12 weeks after anti-TNF therapy. Baseline therefore means pre-treatment synovial sampling in an anti-TNF-starting cohort; it must not be described as post-treatment baseline.",
  "",
  "## Current manuscript reference order",
  "The current local manuscript reference list is: 1 Triaille/GSE225731; 2 Wang/GSE198520; 3 Richard/ECM; 4 Bhamidipati/Xenium; 5 Sisó/spatial prior work.",
  "The first Introduction citation order is Richard (reference 3), Bhamidipati (reference 4), then Sisó (reference 5). The GSE225731 and GSE198520 sources enter the Methods without explicit numbered citations. Thus the current numeric list is not a strict first-in-text-appearance order if numeric citation style is required.",
  "This is a manuscript QA finding only. No manuscript file was modified.",
  "",
  "## Software/database provenance",
  "The expanded M6 CRP–VAS model was executed under R 4.5.2 with readxl 1.4.5, sandwich 3.1-1, lmtest 0.9-40, data.table 1.18.0 and MASS 7.3-65; the coefficient and HC3 interval are retained from the original result file.",
  "The Reactome screen used the local v2023.2.Hs symbols GMT (hash recorded in the CSV). The Stage 2 session records edgeR 4.8.2 and limma 3.66.0, but the historical external Reactome release manifest and the exact edgeR/TMM version used for every frozen earlier run are not fully recoverable.",
  "Plotting versions are not presented as upstream bulk/Xenium analysis versions."
), file.path(res, "REFERENCE_ORDER_AND_PROVENANCE_AUDIT.md"))

writeLines(capture.output(sessionInfo()), file.path(out, "LOGS", "targeted_checks_sessionInfo.txt"))
writeLines(c(
  paste0("Started: ", Sys.time()),
  "Completed all four targeted checks.",
  paste0("GSE198520 component output: ", nrow(component_results), " rows."),
  paste0("Clinical raw Bilateral codes: ", paste(names(table(d$Bilateral)), table(d$Bilateral), collapse = "; ")),
  paste0("M12 auxiliary model rows: ", nrow(aux_models), "."),
  paste0("Completed: ", Sys.time())
), file.path(out, "LOGS", "targeted_checks_execution.log"))

hash_files <- c(
  "RA_TKA_80例.xlsx",
  "RA_TKA论文分析/gate3/raw/GSE198520_Raw_gene_count_matrix.txt.gz",
  "RA_TKA论文分析/SSIP分析/系统桥接/临床桥接_GSE198520/结果/01_GSE198520_sample_manifest.csv",
  "tka_inflammation_redesign/rescue_discovery_20261026/data/c2.cp.reactome.v2023.2.Hs.symbols.gmt",
  "tka_inflammation_redesign/stage6_version_b_integration_20261026/STAGE6_SUPPLEMENTARY_TABLES/GSE225731_PATHWAY_DISCOVERY_ALL.csv",
  "tka_inflammation_redesign/code/15_stage4_crp_vas_20261026.R",
  "tka_inflammation_redesign/stage4_clinical_crp_vas_20261026/STAGE4_CRP_VAS_MODEL_RESULTS.csv",
  "STAGE7T_PRE_SUBMISSION_EVIDENCE_LOCK/CODE/M12_MISSINGNESS_ANALYSIS.R",
  "STAGE7T_PRE_SUBMISSION_EVIDENCE_LOCK/CODE/stage7r_clinical_validity_reanalysis.R",
  "VERSION_C_PRE_REVISION_TARGETED_CHECKS/CODE/targeted_checks.R"
)
hash_dt <- data.table(relative_path = hash_files,
                      exists = file.exists(file.path(root, hash_files)),
                      size_bytes = ifelse(file.exists(file.path(root, hash_files)),
                                          file.info(file.path(root, hash_files))$size, NA),
                      SHA256 = vapply(file.path(root, hash_files),
                                      function(p) if (file.exists(p)) sha256(p) else NA_character_,
                                      character(1)))
fwrite(hash_dt, file.path(out, "INPUT_HASHES", "INPUT_SHA256.tsv"), sep = "\t", na = "NA")

prior_idx <- file.path(root, "STAGE7T_PRE_SUBMISSION_EVIDENCE_LOCK/INPUT_SOURCE_INDEX.csv")
if (file.exists(prior_idx)) {
  prior <- fread(prior_idx, showProgress = FALSE)
  current_clin <- hash_dt[relative_path == "RA_TKA_80例.xlsx"]
  prior_clin <- prior[relative_path == "RA_TKA_80例.xlsx"]
  drift <- data.table(
    file = "RA_TKA_80例.xlsx",
    prior_stage7t_size_bytes = prior_clin$size_bytes,
    current_size_bytes = current_clin$size_bytes,
    prior_stage7t_SHA256 = prior_clin$SHA256,
    current_SHA256 = current_clin$SHA256,
    status = ifelse(identical(as.character(prior_clin$SHA256), as.character(current_clin$SHA256)) &&
                      identical(as.numeric(prior_clin$size_bytes), as.numeric(current_clin$size_bytes)),
                    "NO_DRIFT", "SOURCE_HASH_DRIFT_REQUIRES_AUTHOR_RECONCILIATION")
  )
  fwrite(drift, file.path(out, "INPUT_HASHES", "INPUT_HASH_COMPARISON.csv"))
}



