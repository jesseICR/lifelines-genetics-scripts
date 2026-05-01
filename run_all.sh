#!/bin/bash
# run_all.sh
#
# Convenience entrypoint. Run on the LOGIN node:
#     bash run_all.sh
#
# Steps:
#   1. Downloads snp.info, the dense rsid list, and the public-statgen
#      PCA reference (login node, idempotent). The EUR-balanced PCA
#      reference and centroids_pc6.tsv are shipped with the repo.
#   2. Submits SLURM jobs chained with --dependency=afterok:
#          run_extract.sbatch -> run_relabel.sbatch -> run_dense.sbatch -> run_project.sbatch
#                                                                       +-> run_eur_project.sbatch  (parallel)
#                                                  +-> run_freq.sbatch  (parallel branch off relabel)
#
# Note: 09 + 10 (kinship) are NOT in this pipeline. Run them separately
# after this finishes:
#       bash 09_download_ukb_qc.sh
#       sbatch run_kinship.sbatch
#
# Recommended smoke test before committing to the full pipeline:
#     bash 01_download_snp_info.sh
#     sbatch run_extract.sbatch 22       # chr 22 only, ~5–10 min on compute
#     # inspect logs/extract.<jobid>.log; if happy, then:
#     bash run_all.sh

set -euo pipefail
trap 'echo "[run_all] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

mkdir -p logs

bash 01_download_snp_info.sh
bash 04_download_dense_rsids.sh
bash 06_download_pca_ref.sh

JOB1=$(sbatch --parsable run_extract.sbatch)
[[ -n "$JOB1" ]] || { echo "ERROR: failed to submit run_extract.sbatch" >&2; exit 1; }
echo "Submitted extract job:      $JOB1"

JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" run_relabel.sbatch)
[[ -n "$JOB2" ]] || { echo "ERROR: failed to submit run_relabel.sbatch" >&2; exit 1; }
echo "Submitted relabel job:      $JOB2 (after $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:"$JOB2" run_dense.sbatch)
[[ -n "$JOB3" ]] || { echo "ERROR: failed to submit run_dense.sbatch" >&2; exit 1; }
echo "Submitted dense job:        $JOB3 (after $JOB2)"

JOB4=$(sbatch --parsable --dependency=afterok:"$JOB3" run_project.sbatch)
[[ -n "$JOB4" ]] || { echo "ERROR: failed to submit run_project.sbatch" >&2; exit 1; }
echo "Submitted project job:      $JOB4 (after $JOB3)"

JOB5=$(sbatch --parsable --dependency=afterok:"$JOB2" run_freq.sbatch)
[[ -n "$JOB5" ]] || { echo "ERROR: failed to submit run_freq.sbatch" >&2; exit 1; }
echo "Submitted freq job:         $JOB5 (after $JOB2, parallel to dense/project)"

JOB6=$(sbatch --parsable --dependency=afterok:"$JOB3" run_eur_project.sbatch)
[[ -n "$JOB6" ]] || { echo "ERROR: failed to submit run_eur_project.sbatch" >&2; exit 1; }
echo "Submitted eur-project job:  $JOB6 (after $JOB3, parallel to project)"

echo
echo "Monitor:  squeue -u \$USER"
echo "Logs:     logs/extract.${JOB1}.log"
echo "          logs/relabel.${JOB2}.log"
echo "          logs/dense.${JOB3}.log"
echo "          logs/project.${JOB4}.log"
echo "          logs/freq.${JOB5}.log"
echo "          logs/eur_project.${JOB6}.log"
