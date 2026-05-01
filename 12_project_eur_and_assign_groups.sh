#!/bin/bash
# 12_project_eur_and_assign_groups.sh
#
# For each dense panel (gsa, affymetrix):
#   (a) Project samples onto the EUR-balanced reference PCA via
#       plink2 --score (using eur_pca_balanced.eigenvec.allele +
#       eur_pca_balanced.acount). All 20 PCs preserved in the .sscore.
#   (b) Compute Euclidean distance in PC1..PC6 space from each sample
#       to each gaussian-fit centroid (countries + populations from
#       UKBB / 1000G / HGDP; giab dataset excluded).
#   (c) Write per-sample CSV with the top-10 nearest groups.
#
# Run via SLURM:
#     sbatch run_eur_project.sbatch
#
# Inputs:
#   ./eur_pca_balanced.eigenvec.allele  (shipped with repo)
#   ./eur_pca_balanced.acount           (shipped with repo)
#   ./centroids_pc6.tsv                 (shipped with repo)
#   ./assign_eur_groups.py              (shipped with repo)
#   gsa/dense.{pgen,pvar,psam}          (from 05)
#   affymetrix/dense.{pgen,pvar,psam}   (from 05)
#
# Outputs:
#   gsa/dense_eur_projected.sscore
#   gsa/dense_eur_top10_groups.csv
#   affymetrix/dense_eur_projected.sscore
#   affymetrix/dense_eur_top10_groups.csv
#
# CSV: IID, group_1, distance_1, group_2, distance_2, ..., group_10, distance_10
#       (21 columns total)

set -euo pipefail
trap 'echo "[12] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

ALLELE_FILE=eur_pca_balanced.eigenvec.allele   # shipped with repo
COUNTS_FILE=eur_pca_balanced.acount            # shipped with repo
CENTROIDS_TSV=centroids_pc6.tsv                # shipped with repo
PCA_SNPS=eur_pca_snps.txt                      # generated at runtime (gitignored)
N_PCS_DIST=6      # number of PCs used for distance to centroids
N_PCS_PROJ=20     # number of PCs preserved in .sscore output
TOP_K=10

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
[[ -s "$ALLELE_FILE"   ]] || { echo "ERROR: $ALLELE_FILE missing in cwd (shipped with repo)"   >&2; exit 1; }
[[ -s "$COUNTS_FILE"   ]] || { echo "ERROR: $COUNTS_FILE missing in cwd (shipped with repo)"   >&2; exit 1; }
[[ -s "$CENTROIDS_TSV" ]] || { echo "ERROR: $CENTROIDS_TSV missing in cwd (shipped with repo)" >&2; exit 1; }
[[ -f assign_eur_groups.py ]] || { echo "ERROR: assign_eur_groups.py missing in cwd" >&2; exit 1; }
command -v plink2  >/dev/null 2>&1 || { echo "ERROR: plink2 not on PATH"  >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/dense.pgen" ]] \
        || { echo "ERROR: $panel/dense.pgen missing — run 05 first" >&2; exit 1; }
done
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo "  python3:     $(python3 --version)"
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo

# Build the EUR PCA SNP list once
if [[ ! -s "$PCA_SNPS" ]]; then
    awk 'NR>1 {print $2}' "$ALLELE_FILE" | sort -u > "$PCA_SNPS"
fi
echo "  EUR PCA SNPs: $(wc -l < "$PCA_SNPS") (subset of dense)"

N_GROUPS=$(( $(wc -l < "$CENTROIDS_TSV") - 1 ))
echo "  centroids:    $N_GROUPS groups (from $CENTROIDS_TSV)"

# Locate column numbers in eigenvec.allele dynamically
A1_COL=$(head -1 "$ALLELE_FILE" | tr '\t' '\n' | grep -n '^A1$' | cut -d: -f1)
ID_COL=2
FIRST_PC=$(( A1_COL + 1 ))
LAST_PC=$((  A1_COL + N_PCS_PROJ ))
echo "  --score cols: ID=${ID_COL}, A1=${A1_COL}, scores=${FIRST_PC}-${LAST_PC} (PC1..PC${N_PCS_PROJ})"
echo

for panel in gsa affymetrix; do
    SSCORE="$panel/dense_eur_projected.sscore"
    OUT_CSV="$panel/dense_eur_top${TOP_K}_groups.csv"

    if [[ -f "$OUT_CSV" ]]; then
        echo "[$panel] $OUT_CSV already exists, skipping"
        continue
    fi

    if [[ ! -f "$SSCORE" ]]; then
        echo "[$panel] projecting onto EUR-balanced PCA..."
        plink2 \
            --pfile "$panel/dense" \
            --extract "$PCA_SNPS" \
            --read-freq "$COUNTS_FILE" \
            --score "$ALLELE_FILE" "$ID_COL" "$A1_COL" header-read no-mean-imputation variance-standardize \
            --score-col-nums "${FIRST_PC}-${LAST_PC}" \
            --threads "$THREADS" \
            --memory "$MEM_MB" \
            --out "$panel/dense_eur_projected"

        [[ -f "$SSCORE" ]] \
            || { echo "ERROR: projection did not produce $SSCORE" >&2; exit 1; }
        N_PROJ=$(( $(wc -l < "$SSCORE") - 1 ))
        echo "  [$panel] projected $N_PROJ samples"
    fi

    echo "[$panel] computing top-${TOP_K} group distances (PC1..PC${N_PCS_DIST})..."
    python3 assign_eur_groups.py \
        "$CENTROIDS_TSV" \
        "$SSCORE" \
        "$OUT_CSV" \
        "$TOP_K" \
        "$N_PCS_DIST"
done

echo
echo "=== Done ==="
echo "  GSA:        gsa/dense_eur_projected.sscore"
echo "              gsa/dense_eur_top${TOP_K}_groups.csv"
echo "  Affymetrix: affymetrix/dense_eur_projected.sscore"
echo "              affymetrix/dense_eur_top${TOP_K}_groups.csv"
