#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: apa_landscape_fresh.R <sce_rds> <output_dir>")
}

sce_path <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(sce_path)
rd <- as.data.frame(rowData(sce))
counts_mat <- assay(sce, "counts")

required_cols <- c(
  "gene_name", "utr_name", "is_blacklisted", "is_ipa",
  "utr_position", "utr_length"
)
missing_cols <- setdiff(required_cols, names(rd))
if (length(missing_cols) > 0) {
  stop("Missing required rowData columns: ", paste(missing_cols, collapse = ", "))
}

optional_cols <- c(
  "atlas.utr_type", "mws.cs_usage_class", "is_consistent",
  "is_improper_utr_length", "is_lu", "is_distal",
  "atlas.ncelltypes_gene", "atlas.ncelltypes_utr",
  "atlas.n_utrs_no_ipa", "atlas.pct_utr_no_ipa", "atlas.pct_utr_total",
  "atlas.rank_utr_total", "mws.n_celltypes_gene", "mws.n_celltypes_tx",
  "mws.cs_usage_score", "wt_raw_all", "wt_raw_no_ipa",
  "wt_atlas_all", "wt_atlas_no_ipa", "gene_id", "transcript_id", "merged_txs"
)
for (col in optional_cols) {
  if (!col %in% names(rd)) rd[[col]] <- NA
}

rd$gene_name <- as.character(rd$gene_name)
rd$utr_name <- as.character(rd$utr_name)

# -----------------------------
# Per-sample summary from the SCE
# -----------------------------
cell_counts <- data.frame(
  cell_id = colData(sce)$cell_id,
  sample_id = colData(sce)$sample_id.sq,
  cell_total_umis = as.numeric(Matrix::colSums(counts_mat)),
  stringsAsFactors = FALSE
)

sample_summary <- cell_counts %>%
  group_by(sample_id) %>%
  summarise(
    n_cells = n(),
    total_umis = sum(cell_total_umis, na.rm = TRUE),
    median_umis_per_cell = median(cell_total_umis, na.rm = TRUE),
    mean_umis_per_cell = mean(cell_total_umis, na.rm = TRUE),
    min_umis_per_cell = min(cell_total_umis, na.rm = TRUE),
    max_umis_per_cell = max(cell_total_umis, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_umis))

# -----------------------------
# UTR-level summary
# -----------------------------
utr_tbl <- rd %>%
  mutate(
    total_counts = as.numeric(Matrix::rowSums(counts_mat)),
    detected = total_counts > 0,
    valid_feature = !is_blacklisted & !is_ipa
  ) %>%
  group_by(gene_name) %>%
  mutate(
    gene_total_counts = sum(total_counts, na.rm = TRUE),
    usage_frac = ifelse(gene_total_counts > 0, total_counts / gene_total_counts, NA_real_),
    usage_bin = case_when(
      is.na(usage_frac) ~ NA_character_,
      usage_frac < 0.01 ~ "<1%",
      usage_frac < 0.05 ~ "1-5%",
      usage_frac < 0.10 ~ "5-10%",
      usage_frac < 0.25 ~ "10-25%",
      usage_frac < 0.50 ~ "25-50%",
      TRUE ~ ">=50%"
    ),
    passes_10pct = !is.na(usage_frac) & usage_frac >= 0.10,
    detected_ge10 = total_counts >= 10
  ) %>%
  ungroup()

# -----------------------------
# Gene-level summary
# -----------------------------
gene_tbl <- utr_tbl %>%
  group_by(gene_name) %>%
  summarise(
    n_utrs = n(),
    gene_total_counts = sum(total_counts, na.rm = TRUE),
    n_expressed_utrs = sum(detected, na.rm = TRUE),
    n_blacklisted_utrs = sum(is_blacklisted, na.rm = TRUE),
    n_ipa_utrs = sum(is_ipa, na.rm = TRUE),
    n_valid_utrs = sum(valid_feature, na.rm = TRUE),
    n_valid_expressed_utrs = sum(valid_feature & detected, na.rm = TRUE),
    n_utrs_ge10pct_all = sum(passes_10pct, na.rm = TRUE),
    n_utrs_ge10pct_valid = sum(passes_10pct & valid_feature, na.rm = TRUE),
    n_utrs_ge10_counts = sum(detected_ge10, na.rm = TRUE),
    max_usage_frac = if (all(is.na(usage_frac))) NA_real_ else max(usage_frac, na.rm = TRUE),
    min_usage_frac = if (all(is.na(usage_frac))) NA_real_ else min(usage_frac, na.rm = TRUE),
    any_blacklisted = any(is_blacklisted, na.rm = TRUE),
    any_ipa = any(is_ipa, na.rm = TRUE),
    any_valid_feature = any(valid_feature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(gene_total_counts), desc(n_utrs_ge10pct_valid), desc(n_valid_utrs))

# -----------------------------
# Overall summary
# -----------------------------
overall_summary <- data.frame(
  metric = c(
    "total_utrs",
    "total_detected_utrs",
    "total_blacklisted_utrs",
    "total_ipa_utrs",
    "total_valid_utrs",
    "total_genes",
    "genes_with_ge2_utrs",
    "genes_with_ge2_valid_utrs",
    "genes_with_ge2_valid_and_ge2_ge10pct",
    "genes_with_any_blacklisted",
    "genes_with_any_ipa",
    "median_utrs_per_gene",
    "max_utrs_per_gene",
    "median_gene_total_counts",
    "median_usage_frac",
    "median_detected_utrs_per_gene_ge10counts"
  ),
  value = c(
    nrow(utr_tbl),
    sum(utr_tbl$detected, na.rm = TRUE),
    sum(utr_tbl$is_blacklisted, na.rm = TRUE),
    sum(utr_tbl$is_ipa, na.rm = TRUE),
    sum(utr_tbl$valid_feature, na.rm = TRUE),
    nrow(gene_tbl),
    sum(gene_tbl$n_utrs >= 2, na.rm = TRUE),
    sum(gene_tbl$n_valid_utrs >= 2, na.rm = TRUE),
    sum(gene_tbl$n_valid_utrs >= 2 & gene_tbl$n_utrs_ge10pct_valid >= 2, na.rm = TRUE),
    sum(gene_tbl$any_blacklisted, na.rm = TRUE),
    sum(gene_tbl$any_ipa, na.rm = TRUE),
    median(gene_tbl$n_utrs, na.rm = TRUE),
    max(gene_tbl$n_utrs, na.rm = TRUE),
    median(gene_tbl$gene_total_counts, na.rm = TRUE),
    median(utr_tbl$usage_frac, na.rm = TRUE),
    median(gene_tbl$n_utrs_ge10_counts, na.rm = TRUE)
  )
)

# -----------------------------
# Distribution summaries
# -----------------------------
isoform_distribution_all <- gene_tbl %>%
  count(n_utrs, name = "n_genes") %>%
  arrange(n_utrs)

isoform_distribution_valid <- gene_tbl %>%
  count(n_valid_utrs, name = "n_genes") %>%
  arrange(n_valid_utrs)

usage_bin_summary_all <- utr_tbl %>%
  filter(!is.na(usage_bin)) %>%
  count(usage_bin, name = "n_utrs") %>%
  arrange(factor(usage_bin, levels = c("<1%", "1-5%", "5-10%", "10-25%", "25-50%", ">=50%")))

usage_bin_summary_valid <- utr_tbl %>%
  filter(valid_feature, !is.na(usage_bin)) %>%
  count(usage_bin, name = "n_utrs") %>%
  arrange(factor(usage_bin, levels = c("<1%", "1-5%", "5-10%", "10-25%", "25-50%", ">=50%")))

feature_flag_summary <- data.frame(
  flag = c(
    "is_blacklisted",
    "is_ipa",
    "is_consistent",
    "is_improper_utr_length",
    "is_lu",
    "is_distal"
  ),
  n_features = c(
    sum(utr_tbl$is_blacklisted, na.rm = TRUE),
    sum(utr_tbl$is_ipa, na.rm = TRUE),
    sum(utr_tbl$is_consistent, na.rm = TRUE),
    sum(utr_tbl$is_improper_utr_length, na.rm = TRUE),
    sum(utr_tbl$is_lu, na.rm = TRUE),
    sum(utr_tbl$is_distal, na.rm = TRUE)
  )
)

# -----------------------------
# Gene candidate tables
# -----------------------------
gene_candidates_all <- gene_tbl %>%
  filter(n_utrs >= 2)

gene_candidates_valid <- gene_tbl %>%
  filter(n_valid_utrs >= 2, n_utrs_ge10pct_valid >= 2)

gene_blacklisted_summary <- gene_tbl %>%
  filter(any_blacklisted) %>%
  arrange(desc(n_blacklisted_utrs), desc(gene_total_counts))

gene_ipa_summary <- gene_tbl %>%
  filter(any_ipa) %>%
  arrange(desc(n_ipa_utrs), desc(gene_total_counts))

top_all_apa_genes <- gene_candidates_all %>%
  arrange(desc(gene_total_counts), desc(n_utrs_ge10pct_all), desc(max_usage_frac)) %>%
  head(50)

top_valid_apa_genes <- gene_candidates_valid %>%
  arrange(desc(gene_total_counts), desc(n_utrs_ge10pct_valid), desc(max_usage_frac)) %>%
  head(50)

top_blacklisted_genes <- gene_blacklisted_summary %>%
  head(50)

top_ipa_genes <- gene_ipa_summary %>%
  head(50)

# -----------------------------
# UTR-level tables
# -----------------------------
utr_all <- utr_tbl %>%
  arrange(gene_name, desc(total_counts), utr_position)

utr_valid <- utr_tbl %>%
  filter(valid_feature) %>%
  arrange(gene_name, desc(total_counts), utr_position)

utr_blacklisted <- utr_tbl %>%
  filter(is_blacklisted) %>%
  arrange(gene_name, desc(total_counts), utr_position)

utr_ipa <- utr_tbl %>%
  filter(is_ipa) %>%
  arrange(gene_name, desc(total_counts), utr_position)

# -----------------------------
# Write outputs
# -----------------------------
write.csv(overall_summary, file.path(outdir, "overall_summary.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(outdir, "sample_summary.csv"), row.names = FALSE)
write.csv(feature_flag_summary, file.path(outdir, "feature_flag_summary.csv"), row.names = FALSE)

write.csv(gene_tbl, file.path(outdir, "gene_summary_all.csv"), row.names = FALSE)
write.csv(utr_all, file.path(outdir, "utr_summary_all.csv"), row.names = FALSE)

write.csv(isoform_distribution_all, file.path(outdir, "isoform_distribution_all.csv"), row.names = FALSE)
write.csv(isoform_distribution_valid, file.path(outdir, "isoform_distribution_valid.csv"), row.names = FALSE)

write.csv(usage_bin_summary_all, file.path(outdir, "usage_bin_summary_all.csv"), row.names = FALSE)
write.csv(usage_bin_summary_valid, file.path(outdir, "usage_bin_summary_valid.csv"), row.names = FALSE)

write.csv(gene_candidates_all, file.path(outdir, "gene_candidates_all_ge2isoforms.csv"), row.names = FALSE)
write.csv(gene_candidates_valid, file.path(outdir, "gene_candidates_valid_ge2isoforms_ge10pct.csv"), row.names = FALSE)

write.csv(gene_blacklisted_summary, file.path(outdir, "gene_summary_blacklisted.csv"), row.names = FALSE)
write.csv(gene_ipa_summary, file.path(outdir, "gene_summary_ipa.csv"), row.names = FALSE)

write.csv(utr_valid, file.path(outdir, "utr_summary_valid.csv"), row.names = FALSE)
write.csv(utr_blacklisted, file.path(outdir, "utr_summary_blacklisted.csv"), row.names = FALSE)
write.csv(utr_ipa, file.path(outdir, "utr_summary_ipa.csv"), row.names = FALSE)

write.csv(top_all_apa_genes, file.path(outdir, "top_apa_genes_all.csv"), row.names = FALSE)
write.csv(top_valid_apa_genes, file.path(outdir, "top_apa_genes_valid.csv"), row.names = FALSE)
write.csv(top_blacklisted_genes, file.path(outdir, "top_blacklisted_genes.csv"), row.names = FALSE)
write.csv(top_ipa_genes, file.path(outdir, "top_ipa_genes.csv"), row.names = FALSE)

cat("Saved outputs to: ", outdir, "\n", sep = "")
cat("\n--- Overall summary ---\n")
print(overall_summary)
cat("\n--- Sample summary ---\n")
print(sample_summary)
cat("\n--- Feature flag summary ---\n")
print(feature_flag_summary)
cat("\n--- Top valid APA genes ---\n")
print(top_valid_apa_genes)
cat("\n--- Top IPA genes ---\n")
print(top_ipa_genes)
cat("\n--- Top blacklisted genes ---\n")
print(top_blacklisted_genes)
