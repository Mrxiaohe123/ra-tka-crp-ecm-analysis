# PUBLIC RELEASE CANDIDATE COPY
# Source: tka_inflammation_redesign/rescue_discovery_20261026/code/01_GSE225731_pathway_discovery.R
# Source SHA256: c7443f8b444a14193d176780813e9752780005a46494ed06b4b87f8c86431760
# Purpose: Public-data pathway discovery source; not rerun during packaging
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, encoding = "UTF-8")
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

root <- Sys.getenv("VERSION_C_PROJECT_ROOT", ".")
out <- Sys.getenv("VERSION_C_DISCOVERY_OUT", file.path(root, "rescue_discovery_20261026"))
matrix_file <- file.path(root, "gse225731_external_20261025", "data", "GSE225731_raw_counts_74RA_samples.txt.gz")
meta_file <- file.path(root, "GSE225731_BACKUP_READINESS_sample_metadata.csv")
hallmark_file <- file.path(out, "data", "h.all.v2023.2.Hs.symbols.gmt")
reactome_file <- file.path(out, "data", "c2.cp.reactome.v2023.2.Hs.symbols.gmt")

stopifnot(file.exists(matrix_file), file.exists(meta_file), file.exists(hallmark_file), file.exists(reactome_file))

read_gmt <- function(path, collection) {
  ln <- readLines(path, warn = FALSE)
  x <- lapply(ln, function(z) {
    p <- strsplit(z, "\t", fixed = FALSE)[[1L]]
    list(name = p[[1L]], genes = unique(p[-c(1L, 2L)]))
  })
  data.table(collection = collection,
             pathway = vapply(x, `[[`, character(1), "name"),
             genes = I(lapply(x, `[[`, "genes")))
}

gmt <- rbind(read_gmt(hallmark_file, "Hallmark"), read_gmt(reactome_file, "Reactome"))
gmt <- gmt[!duplicated(pathway)]

mat <- fread(cmd = paste("gzip -cd", shQuote(matrix_file)), check.names = FALSE, showProgress = FALSE)
gene_id <- sub("\\..*$", "", as.character(mat[[1L]]))
sample_ids <- names(mat)[-1L]
counts <- as.matrix(mat[, -1L, with = FALSE])
storage.mode(counts) <- "numeric"
stopifnot(ncol(counts) == 74L, !anyDuplicated(gene_id), !anyDuplicated(sample_ids),
          all(is.finite(counts)), all(counts >= 0), max(abs(counts - round(counts))) < 1e-8)

meta <- fread(meta_file, na.strings = c("", "NA"), showProgress = FALSE)
stopifnot("sample_title" %in% names(meta), "crp_bl (mg/dl)" %in% names(meta))
meta <- meta[match(sample_ids, sample_title)]
stopifnot(all(meta$sample_title == sample_ids))
crp <- suppressWarnings(as.numeric(meta[["crp_bl (mg/dl)"]]))
stopifnot(all(is.finite(crp)))

# Fixed transformation: TMM -> log2(CPM+1) -> Ensembl-to-symbol aggregation by mean
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
log_cpm1 <- log2(cpm(dge, normalized.lib.sizes = TRUE) + 1)

ann <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(gene_id),
                             columns = "SYMBOL", keytype = "ENSEMBL")
ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), ]
ann <- ann[!duplicated(paste(ann$ENSEMBL, ann$SYMBOL)), ]
map_sym <- ann$SYMBOL[match(gene_id, ann$ENSEMBL)]
keep <- !is.na(map_sym) & nzchar(map_sym)
expr <- rowsum(log_cpm1[keep, , drop = FALSE], group = map_sym[keep], reorder = FALSE) /
  as.numeric(table(map_sym[keep])[rownames(rowsum(log_cpm1[keep, , drop = FALSE], group = map_sym[keep], reorder = FALSE))])
expr <- as.matrix(expr)
rownames(expr) <- make.unique(rownames(expr))
gene_sd <- apply(expr, 1L, sd)
z <- expr[gene_sd > 0, , drop = FALSE]
z <- t(scale(t(z)))

min_genes <- 10L
min_fraction <- 0.20
rank_one <- function(x) {
  if (length(unique(x)) < 2L) return(NA_real_)
  suppressWarnings(cor(x, crp, method = "spearman"))
}
fisher_ci <- function(r, n) {
  if (!is.finite(r) || n < 6L || abs(r) >= 1) return(c(NA_real_, NA_real_))
  zr <- atanh(r); se <- 1 / sqrt(n - 3L)
  tanh(zr + c(-1, 1) * qnorm(0.975) * se)
}

rows <- vector("list", nrow(gmt))
for (i in seq_len(nrow(gmt))) {
  genes <- gmt$genes[[i]]
  present <- intersect(genes, rownames(z))
  coverage <- length(present) / length(unique(genes))
  eligible <- length(present) >= min_genes && coverage >= min_fraction
  if (!eligible) {
    rows[[i]] <- data.table(collection = gmt$collection[i], pathway = gmt$pathway[i],
      pathway_gene_N = length(unique(genes)), mapped_gene_N = length(present),
      mapped_fraction = coverage, analysis_status = "EXCLUDED_LOW_COVERAGE",
      N = NA_integer_, rho = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      P_value = NA_real_, loo_rho_min = NA_real_, loo_rho_max = NA_real_,
      loo_same_direction_N = NA_integer_, score_method = "TMM_log2CPM1_gene_z_mean")
    next
  }
  score <- colMeans(z[present, , drop = FALSE])
  rho <- rank_one(score)
  p <- suppressWarnings(cor.test(score, crp, method = "spearman", exact = FALSE)$p.value)
  ci <- fisher_ci(rho, length(crp))
  loo <- vapply(seq_along(score), function(j) {
    if (length(unique(score[-j])) < 2L || length(unique(crp[-j])) < 2L) return(NA_real_)
    suppressWarnings(cor(score[-j], crp[-j], method = "spearman"))
  }, numeric(1L))
  same_dir <- if (is.finite(rho) && rho != 0) sum(sign(loo) == sign(rho), na.rm = TRUE) else NA_integer_
  rows[[i]] <- data.table(collection = gmt$collection[i], pathway = gmt$pathway[i],
    pathway_gene_N = length(unique(genes)), mapped_gene_N = length(present),
    mapped_fraction = coverage, analysis_status = "TESTED", N = length(crp), rho = rho,
    CI_low = ci[[1L]], CI_high = ci[[2L]], P_value = p,
    loo_rho_min = min(loo, na.rm = TRUE), loo_rho_max = max(loo, na.rm = TRUE),
    loo_same_direction_N = same_dir, score_method = "TMM_log2CPM1_gene_z_mean")
}
res <- rbindlist(rows, fill = TRUE)
res[, q_value_BH := NA_real_]
res[analysis_status == "TESTED", q_value_BH := p.adjust(P_value, method = "BH")]
res[, CRP_zero_N := sum(crp == 0)]
res[, CRP_missing_N := sum(!is.finite(crp))]
setorder(res, q_value_BH, P_value)
fwrite(res, file.path(out, "GSE225731_PATHWAY_DISCOVERY_ALL.csv"))

tested <- res[analysis_status == "TESTED"]
stable <- tested[loo_same_direction_N >= ceiling(0.80 * N)]
top <- stable[order(q_value_BH, -abs(rho), -mapped_fraction)][1:min(20L, .N)]

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
fmtp <- function(x) ifelse(is.na(x), "NA", ifelse(x < .001, "<0.001", formatC(x, digits = 4, format = "f")))

report <- c(
  "# Exploratory multiple-testing report",
  "",
  "## Material Passport",
  "- Type: exploratory pathway discovery execution",
  "- Status: complete; no historical output overwritten",
  paste0("- Input: GSE225731 raw counts, ", nrow(counts), " Ensembl rows x ", ncol(counts), " samples"),
  paste0("- Gene-set inputs: MSigDB v2023.2.Hs Hallmark + Reactome; total sets = ", nrow(gmt)),
  paste0("- Score: TMM-normalized log2(CPM+1), gene-wise z-score across 74 samples, mean over mapped genes"),
  paste0("- Eligibility: mapped genes >= ", min_genes, " and mapped fraction >= ", min_fraction),
  "- Association: two-sided Spearman correlation with continuous baseline CRP; raw CRP zeros retained; no pseudocount",
  "- Multiple testing: BH adjustment across tested gene sets jointly",
  "- CI: Fisher-transformed approximate 95% CI for the Spearman coefficient",
  "- Stability: leave-one-sample-out rho range and same-direction count; diagnostic only",
  "",
  "## Summary",
  paste0("- Tested = ", nrow(tested), "; excluded for low coverage = ", nrow(res[analysis_status == "EXCLUDED_LOW_COVERAGE"]), "."),
  paste0("- CRP zero values retained = ", sum(crp == 0), "; CRP missing = ", sum(!is.finite(crp)), "."),
  "- All tested pathways, including null and opposite-direction results, are stored in the CSV.",
  "- This is discovery, not an independent validation or preregistered primary analysis.",
  "",
  "## Stable triage list",
  "The following is a transparent post-hoc priority triage requiring independent validation; it is not a new confirmatory endpoint.",
  "",
  "| Collection | Pathway | Mapped/N | rho | 95% CI | P | q | LOO same direction |",
  "|---|---|---:|---:|---|---:|---:|---:|",
  vapply(seq_len(nrow(top)), function(i) paste0("| ", top$collection[i], " | ", top$pathway[i], " | ", top$mapped_gene_N[i], "/", top$pathway_gene_N[i], " | ", fmt(top$rho[i]), " | ", fmt(top$CI_low[i]), "–", fmt(top$CI_high[i]), " | ", fmtp(top$P_value[i]), " | ", fmtp(top$q_value_BH[i]), " | ", top$loo_same_direction_N[i], "/", top$N[i], " |"), character(1L)),
  "",
  "## Interpretation boundary",
  "Top-ranked pathways are candidates for biological review only. They are not automatically drug targets, spatially validated programs, or evidence that TKA changes the pathway. The candidate list cannot be used to rewrite the fixed IL6/JAK/STAT3 result.")
writeLines(report, file.path(out, "EXPLORATORY_MULTIPLE_TESTING_REPORT.md"))

priority <- c(
  "# Candidate pathway priority",
  "",
  "This priority list is an exploratory triage created after the complete prespecified scan. It combines BH q-value, effect magnitude, coverage and leave-one-out direction stability; no pathway was selected solely by a single nominal P value.",
  "",
  paste0("Ranked stable candidates retained = ", nrow(top), "; all other pathways remain in `GSE225731_PATHWAY_DISCOVERY_ALL.csv`.")
)
if (nrow(top)) priority <- c(priority, "", paste0(seq_len(nrow(top)), ". `", top$pathway, "` — ", top$collection, "; rho=", fmt(top$rho), "; q=", fmtp(top$q_value_BH), "; mapped ", top$mapped_gene_N, "/", top$pathway_gene_N, "; LOO same direction ", top$loo_same_direction_N, "/", top$N, "."))
priority <- c(priority, "", "Next-stage use requires: external biological evidence, cell-type/spatial assessment, target tractability review, and independent validation. No docking or virtual intervention is authorized by this file alone.")
writeLines(priority, file.path(out, "CANDIDATE_PATHWAY_PRIORITY.md"))

writeLines(capture.output(sessionInfo()), file.path(out, "01_GSE225731_pathway_discovery_sessionInfo.txt"))
