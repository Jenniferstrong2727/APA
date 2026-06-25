#!/bin/bash
set -euo pipefail

SAMTOOLS="/hpc/packages/minerva-rocky9/samtools/1.21/bin/samtools"

# Run-specific settings
RUN="GEX-10-1"
NUCLEI_DIR="/sc/arion/projects/ad-omics/Jennifer/in_house_mg_data/nuclei"
BAM_IN="$NUCLEI_DIR/bams/${RUN}.bam"
BARCODE_DIR="$NUCLEI_DIR/bam_prep/${RUN}/barcodes"
OUT_DIR="$NUCLEI_DIR/bam_prep/${RUN}/donor_bams"
LOG_DIR="$NUCLEI_DIR/bam_prep/${RUN}/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Only the two remaining donor-region samples for GEX-10-1
SAMPLES=(
  "NBB19_071_AMY"
  "NBB19_074_AMY"
)

for SAMPLE in "${SAMPLES[@]}"; do
  BC_FILE="$BARCODE_DIR/${RUN}__${SAMPLE}.txt"
  OUT_BAM="$OUT_DIR/${RUN}__${SAMPLE}.bam"
  LOG_FILE="$LOG_DIR/${RUN}__${SAMPLE}.log"

  echo "[$(date)] Starting $SAMPLE" | tee "$LOG_FILE"
  echo "Barcode file: $BC_FILE" | tee -a "$LOG_FILE"
  echo "Output BAM:   $OUT_BAM" | tee -a "$LOG_FILE"

  if [[ -f "$OUT_BAM" && -f "${OUT_BAM}.bai" ]]; then
    echo "[$(date)] Skipping $SAMPLE (BAM and BAI already exist)" | tee -a "$LOG_FILE"
    continue
  fi

  "$SAMTOOLS" view -b -D CB:"$BC_FILE" "$BAM_IN" -o "$OUT_BAM" >> "$LOG_FILE" 2>&1
  "$SAMTOOLS" index "$OUT_BAM" >> "$LOG_FILE" 2>&1
  "$SAMTOOLS" quickcheck -v "$OUT_BAM" >> "$LOG_FILE" 2>&1

  echo "[$(date)] Done $SAMPLE" | tee -a "$LOG_FILE"
  echo "" | tee -a "$LOG_FILE"
done
