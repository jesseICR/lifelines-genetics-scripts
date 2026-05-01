#!/bin/bash
# run_kinship.sh
#
# Kinship pipeline entrypoint. Run on the LOGIN node, AFTER run_all.sh
# has finished:
#     bash run_kinship.sh
#
# Steps:
#   1. Downloads ukb_snp_qc.txt from biobank.ndph.ox.ac.uk (login,
#      idempotent).
#   2. If merged_kinship.kin0 already exists, skips the SLURM submission.
#      Otherwise submits run_kinship.sbatch which:
#        - subsets each panel's dense pfile to UKB relatedness rsids
#          and converts to bfile (plink2)
#        - merges the two bfiles with plink1 --bmerge (auto-exclude
#          allele-mismatched SNPs and retry, mirroring the user's
#          merge_wgs_imputed.sh pattern)
#        - runs plink2 --make-king-table on the merged bfile
#
# Output: merged_kinship.{bed,bim,fam,kin0}
# (Cross-panel relationships are visible in the kinship table — they
# would otherwise be missed if GSA and Affymetrix were run separately.)

set -euo pipefail
trap 'echo "[run_kinship] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

mkdir -p logs

bash 09_download_ukb_qc.sh

if [[ -f merged_kinship.kin0 ]]; then
    echo "[skip] kinship: merged_kinship.kin0 already exists"
    echo "All outputs already present. Nothing submitted."
    exit 0
fi

JOB=$(sbatch --parsable run_kinship.sbatch)
[[ -n "$JOB" ]] || { echo "ERROR: failed to submit run_kinship.sbatch" >&2; exit 1; }
echo "Submitted kinship job: $JOB"

echo
echo "Monitor:  squeue -u \$USER"
echo "Log:      logs/kinship.${JOB}.log"
