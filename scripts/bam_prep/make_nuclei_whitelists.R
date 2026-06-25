library(readxl)

# Paths
nuclei_dir <- "/sc/arion/projects/ad-omics/Jennifer/in_house_mg_data/nuclei"
meta_file  <- file.path(nuclei_dir, "metadata", "nuclei_all_metadata.csv")
master_file <- file.path(nuclei_dir, "metadata", "Nuclei_metadata.xlsx")

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta$barcode <- sub("\\.1$", "-1", meta$barcode)

# Keep only usable rows
meta <- meta[!is.na(meta$run) &
               !is.na(meta$donor_id_relaxed) &
               meta$donor_id_relaxed != "unassigned" &
               !is.na(meta$barcode),
             c("barcode", "run", "donor_id_relaxed")]

# Read master run/sample sheet
master <- as.data.frame(read_excel(master_file))
master <- unique(master[, c("donor_key", "donor_key_region", "run")])

# Merge per-nucleus assignments with run/sample mapping
merged <- merge(meta, master,
                by.x = c("run", "donor_id_relaxed"),
                by.y = c("run", "donor_key"),
                all = FALSE)

# Prepare output base
barcodes_root <- file.path(nuclei_dir, "bam_prep")

# Create one whitelist per run + donor_key_region
samples <- unique(merged[, c("run", "donor_key_region")])
samples <- samples[order(samples$run, samples$donor_key_region), ]

for (i in seq_len(nrow(samples))) {
  run_i <- samples$run[i]
  region_i <- samples$donor_key_region[i]

  out_dir <- file.path(barcodes_root, run_i, "barcodes")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_file <- file.path(out_dir, paste0(run_i, "__", region_i, ".txt"))

  # Skip if the file already exists and is non-empty
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    message("Skipping existing file: ", out_file)
    next
  }

  bcs <- sort(unique(merged$barcode[merged$run == run_i &
                                    merged$donor_key_region == region_i]))

  writeLines(bcs, out_file)
  message("Wrote ", length(bcs), " barcodes to ", out_file)
}

# Write a quick summary
summary_file <- file.path(barcodes_root, "manifests", "whitelist_summary.tsv")
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)

summary_df <- aggregate(barcode ~ run + donor_key_region, data = merged,
                        FUN = function(x) length(unique(x)))
names(summary_df)[names(summary_df) == "barcode"] <- "n_barcodes"
write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote summary to: ", summary_file)
