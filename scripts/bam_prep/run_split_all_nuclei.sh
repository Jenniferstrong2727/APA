#!/bin/bash
set -euo pipefail

SAMTOOLS="/hpc/packages/minerva-rocky9/samtools/1.21/bin/samtools"
NUCLEI_DIR="/sc/arion/projects/ad-omics/Jennifer/in_house_mg_data/nuclei"
MANIFEST="$NUCLEI_DIR/bam_prep/manifests/nuclei_bam_manifest.tsv"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run donor_key donor_key_region n_barcodes bam_in barcode_file out_bam log_file; do
  mkdir -p "$(dirname "$out_bam")" "$(dirname "$log_file")"

  echo "[$(date)] Starting ${run} ${donor_key_region}" | tee "$log_file"
  echo "Run:          $run" | tee -a "$log_file"
  echo "Donor key:    $donor_key" | tee -a "$log_file"
  echo "Region:       $donor_key_region" | tee -a "$log_file"
  echo "Barcodes:     $n_barcodes" | tee -a "$log_file"
  echo "Input BAM:    $NUCLEI_DIR/$bam_in" | tee -a "$log_file"
  echo "Barcode file: $NUCLEI_DIR/$barcode_file" | tee -a "$log_file"
  echo "Output BAM:   $NUCLEI_DIR/$out_bam" | tee -a "$log_file"

  if [[ -f "$NUCLEI_DIR/$out_bam" && -f "$NUCLEI_DIR/${out_bam}.bai" ]]; then
    echo "[$(date)] Skipping ${run} ${donor_key_region} (BAM and BAI already exist)" | tee -a "$log_file"
    echo "" | tee -a "$log_file"
    continue
  fi

  "$SAMTOOLS" view -b -D CB:"$NUCLEI_DIR/$barcode_file" "$NUCLEI_DIR/$bam_in" -o "$NUCLEI_DIR/$out_bam" >> "$log_file" 2>&1
  "$SAMTOOLS" index "$NUCLEI_DIR/$out_bam" >> "$log_file" 2>&1
  "$SAMTOOLS" quickcheck -v "$NUCLEI_DIR/$out_bam" >> "$log_file" 2>&1

  echo "[$(date)] Done ${run} ${donor_key_region}" | tee -a "$log_file"
  echo "" | tee -a "$log_file"
done
