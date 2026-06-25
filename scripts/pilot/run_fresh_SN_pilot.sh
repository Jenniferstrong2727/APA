#!/bin/bash

# Load conda
source /hpc/users/stronj04/miniconda3/etc/profile.d/conda.sh

# Activate environment
conda activate scutr_smk

# Thread limits
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Go to scUTRquant repo
cd /sc/arion/projects/ad-omics/Jennifer/APA_MG/scUTRquant

# Run pipeline
snakemake --use-conda \
          --conda-frontend conda \
          --conda-prefix /sc/arion/work/stronj04/pineda_apa_test/scUTRquant/.snakemake/conda \
          --configfile /sc/arion/projects/ad-omics/Jennifer/APA_MG/configs/fresh_SN_pilot.yaml \
          --cores 8 \
          --printshellcmds
