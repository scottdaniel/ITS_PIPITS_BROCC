#!/usr/bin/env bash

#SBATCH --mem=2G
#SBATCH -n 1
#SBATCH --export=ALL
#SBATCH --mail-user=danielsg@chop.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --no-requeue
#SBATCH -t 12:00:00
#SBATCH --output=slurm_%x_%j.out

#Uncomment the next two lines if you want to 'qsub' this script
source ~/.bashrc.conda #needed to make "conda" command to work
conda activate ITS_PIPITS_BROCC

set -xeuo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: bash $0 PATH_TO_CONFIG"
    exit 1
fi

CONFIG_FP=$1

snakemake \
    --jobs 100 \
    --configfile ${CONFIG_FP} \
    --cluster-config cluster.json \
    --keep-going \
    --latency-wait 90 \
    --notemp \
    --printshellcmds \
    --cluster \
    "sbatch --no-requeue --export=ALL --mem={cluster.mem_free} -n {threads}"
