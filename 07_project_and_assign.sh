#!/bin/bash
# 07_project_and_assign.sh
#
# Step 7 of 7. For each dense panel (gsa, affymetrix):
#   (a) Project samples onto the public-statgen reference PCA via
#       plink2 --score (using pca_pcs.eigenvec.allele + pca_counts.acount).
#   (b) Assign each IID to the nearest super-population centroid by
#       Euclidean distance in PC1..PC20 space.
#
# Run via SLURM:
#     sbatch run_project.sbatch
#
# Inputs (from previous steps):
#   pca_ref/pca_pcs.eigenvec.allele       (from 06)
#   pca_ref/pca_counts.acount             (from 06)
#   pca_ref/centroids_top20_supervised.tsv (from 06)
#   gsa/dense.{pgen,pvar,psam}            (from 05)
#   affymetrix/dense.{pgen,pvar,psam}     (from 05)
#   ./assign_superpop.py                  (in repo)
#
# Outputs:
#   gsa/dense_projected.sscore
#   gsa/dense_assigned_superpop.csv
#   affymetrix/dense_projected.sscore
#   affymetrix/dense_assigned_superpop.csv
#
# CSV columns: FID,IID,supervised_group,euclidean_distance

set -euo pipefail
trap 'echo "[07] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

PCA_REF=pca_ref
ALLELE_FILE="$PCA_REF/pca_pcs.eigenvec.allele"
COUNTS_FILE="$PCA_REF/pca_counts.acount"
CENTROIDS_FILE="$PCA_REF/centroids_top20_supervised.tsv"
PCA_SNPS="$PCA_REF/pca_snps.txt"
N_PCS=20

THREADS=${SLURM_CPUS_PER_TASK:-2}
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && (( SLURM_MEM_PER_NODE > 2048 )); then
    MEM_MB=$(( SLURM_MEM_PER_NODE - 1024 ))
else
    MEM_MB=${PLINK_MEM_MB:-8000}
fi

module load PLINK/2.0-alpha6.20-20250707
module load Python/3.12.3-GCCcore-13.3.0

echo "=== Pre-flight ==="
echo "  hostname:    $(hostname)"
echo "  job:         ${SLURM_JOB_ID:-<not under SLURM>}"
echo "  cwd:         $(pwd)"

[[ -s "$ALLELE_FILE"     ]] || { echo "ERROR: $ALLELE_FILE missing — run 06_download_pca_ref.sh"     >&2; exit 1; }
[[ -s "$COUNTS_FILE"     ]] || { echo "ERROR: $COUNTS_FILE missing — run 06_download_pca_ref.sh"     >&2; exit 1; }
[[ -s "$CENTROIDS_FILE"  ]] || { echo "ERROR: $CENTROIDS_FILE missing — run 06_download_pca_ref.sh"  >&2; exit 1; }
[[ -f assign_superpop.py ]] || { echo "ERROR: assign_superpop.py missing in cwd"                     >&2; exit 1; }
command -v plink2  >/dev/null 2>&1 || { echo "ERROR: plink2 not on PATH after module load"  >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not on PATH after module load" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/dense.pgen" ]] \
        || { echo "ERROR: $panel/dense.pgen missing — run 05_extract_dense.sh first" >&2; exit 1; }
done
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo "  python3:     $(python3 --version)"
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo

# Build the PCA SNP list (one rsid per variant) once
if [[ ! -s "$PCA_SNPS" ]]; then
    awk 'NR>1 {print $2}' "$ALLELE_FILE" | sort -u > "$PCA_SNPS"
fi
echo "  PCA reference: $(wc -l < "$PCA_SNPS") SNPs (centroids on PC1..PC${N_PCS})"

# Locate column numbers in eigenvec.allele dynamically
A1_COL=$(head -1 "$ALLELE_FILE" | tr '\t' '\n' | grep -n '^A1$' | cut -d: -f1)
ID_COL=2
FIRST_PC=$(( A1_COL + 1 ))
LAST_PC=$((  A1_COL + N_PCS ))
echo "  --score cols: ID=${ID_COL}, A1=${A1_COL}, scores=${FIRST_PC}-${LAST_PC}"
echo

for panel in gsa affymetrix; do
    OUT_CSV="$panel/dense_assigned_superpop.csv"
    SSCORE="$panel/dense_projected.sscore"

    if [[ -f "$OUT_CSV" ]]; then
        echo "[$panel] $OUT_CSV already exists, skipping"
        continue
    fi

    if [[ ! -f "$SSCORE" ]]; then
        echo "[$panel] projecting onto reference PCA..."
        plink2 \
            --pfile "$panel/dense" \
            --extract "$PCA_SNPS" \
            --read-freq "$COUNTS_FILE" \
            --score "$ALLELE_FILE" "$ID_COL" "$A1_COL" header-read no-mean-imputation variance-standardize \
            --score-col-nums "${FIRST_PC}-${LAST_PC}" \
            --threads "$THREADS" \
            --memory "$MEM_MB" \
            --out "$panel/dense_projected"

        [[ -f "$SSCORE" ]] \
            || { echo "ERROR: projection did not produce $SSCORE" >&2; exit 1; }

        N_PROJ=$(( $(wc -l < "$SSCORE") - 1 ))
        echo "  [$panel] projected $N_PROJ samples"
    fi

    echo "[$panel] assigning super-population by nearest centroid..."
    python3 assign_superpop.py \
        "$CENTROIDS_FILE" \
        "$SSCORE" \
        "$OUT_CSV" \
        "$N_PCS"
done

echo
echo "=== Done ==="
echo "  GSA:        gsa/dense_projected.sscore"
echo "              gsa/dense_assigned_superpop.csv"
echo "  Affymetrix: affymetrix/dense_projected.sscore"
echo "              affymetrix/dense_assigned_superpop.csv"
