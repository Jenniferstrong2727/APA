#!/bin/bash
set -euo pipefail

module purge
module load R/4.2.0

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Using Rscript: $(which Rscript)"
Rscript --version

echo "Running pseudobulk APA landscape with:"
echo "  input:  /sc/arion/projects/ad-omics/Jennifer/APA_MG/scUTRquant/data/sce/utrome_hg38_v1/fresh_SN_pilot.txs.Rds"
echo "  output: /sc/arion/projects/ad-omics/Jennifer/APA_MG/results/Cells_SN"

Rscript --vanilla /sc/arion/projects/ad-omics/Jennifer/APA_MG/scripts/apa_landscape_pseudobulk.R \
  /sc/arion/projects/ad-omics/Jennifer/APA_MG/scUTRquant/data/sce/utrome_hg38_v1/fresh_SN_pilot.txs.Rds \
  /sc/arion/projects/ad-omics/Jennifer/APA_MG/results/Cells_SN
