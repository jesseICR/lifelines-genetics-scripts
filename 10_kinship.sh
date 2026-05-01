#!/bin/bash
# 10_kinship.sh
#
# For each dense panel (gsa, affymetrix):
#   (a) Extract from ukb_snp_qc.txt the rsids with in_Relatedness == 1.
#   (b) Subset the dense pfile to those rsids.
#   (c) Run plink2 --make-king-table to produce a .kin0 file.
#
# Run via SLURM:
#     sbatch run_kinship.sbatch
#
# Inputs:
#   ./ukb_snp_qc.txt                  (from 09_download_ukb_qc.sh)
#   gsa/dense.{pgen,pvar,psam}        (from 05_extract_dense.sh)
#   affymetrix/dense.{pgen,pvar,psam}
#
# Outputs (per panel):
#   $panel/dense_kinship.{pgen,pvar,psam}   subsetted pfile
#   $panel/dense_kinship.kin0               KING kinship table
#
# Tunable via env:
#   KING_CUTOFF   minimum kinship reported in .kin0 (default 0.0884,
#                 the standard 3rd-degree-relative threshold)

set -euo pipefail
trap 'echo "[10] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

UKB_QC=ukb_snp_qc.txt
RELATEDNESS_RSIDS=ukb_relatedness_rsids.txt
KING_CUTOFF=${KING_CUTOFF:-0.0884}

THREADS=${SLURM_CPUS_PER_TASK:-2}
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && (( SLURM_MEM_PER_NODE > 2048 )); then
    MEM_MB=$(( SLURM_MEM_PER_NODE - 1024 ))
else
    MEM_MB=${PLINK_MEM_MB:-12000}
fi

module load PLINK/2.0-alpha6.20-20250707

echo "=== Pre-flight ==="
echo "  hostname:    $(hostname)"
echo "  job:         ${SLURM_JOB_ID:-<not under SLURM>}"
echo "  cwd:         $(pwd)"
[[ -s "$UKB_QC" ]] \
    || { echo "ERROR: $UKB_QC missing — run 09_download_ukb_qc.sh first" >&2; exit 1; }
command -v plink2 >/dev/null 2>&1 \
    || { echo "ERROR: plink2 not on PATH after module load" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/dense.pgen" ]] \
        || { echo "ERROR: $panel/dense.pgen missing — run 05_extract_dense.sh first" >&2; exit 1; }
done
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo "  king cutoff: $KING_CUTOFF"
echo

# -- step 1: build the relatedness rsid list (idempotent) --
if [[ -s "$RELATEDNESS_RSIDS" ]]; then
    echo "=== Relatedness rsid list already built ==="
    echo "  $RELATEDNESS_RSIDS: $(wc -l < "$RELATEDNESS_RSIDS") rsids"
else
    echo "=== Extracting in_Relatedness rsids from $UKB_QC ==="
    REL_COL=$(head -1 "$UKB_QC" | tr ' ' '\n' | grep -n '^in_Relatedness$' | cut -d: -f1)
    [[ -n "$REL_COL" ]] \
        || { echo "ERROR: in_Relatedness column not found in $UKB_QC" >&2; exit 1; }
    echo "  in_Relatedness is column $REL_COL"

    awk -v col="$REL_COL" 'NR > 1 && $col == 1 {print $1}' "$UKB_QC" \
        | sort -u > "$RELATEDNESS_RSIDS"

    [[ -s "$RELATEDNESS_RSIDS" ]] \
        || { echo "ERROR: no rsids with in_Relatedness == 1" >&2; exit 1; }
    echo "  wrote $(wc -l < "$RELATEDNESS_RSIDS") rsids to $RELATEDNESS_RSIDS"
fi
echo

# -- step 2 + 3: subset each dense pfile and compute KING kinship --
for panel in gsa affymetrix; do
    BASE="$panel/dense_kinship"
    KIN0="$BASE.kin0"

    echo "=== $panel ==="

    if [[ -f "$BASE.pgen" ]]; then
        echo "  [$panel] $BASE.pgen already exists, skipping subset"
    else
        echo "  [$panel] subsetting dense to UKB relatedness rsids..."
        plink2 \
            --pfile "$panel/dense" \
            --extract "$RELATEDNESS_RSIDS" \
            --threads "$THREADS" \
            --memory "$MEM_MB" \
            --make-pgen \
            --out "$BASE"

        [[ -f "$BASE.pgen" ]] \
            || { echo "ERROR: $BASE.pgen not produced" >&2; exit 1; }
        KEPT=$(( $(wc -l < "$BASE.pvar") - $(grep -c '^#' "$BASE.pvar") ))
        echo "  [$panel] kinship pfile has $KEPT variants"
    fi

    if [[ -f "$KIN0" ]]; then
        echo "  [$panel] $KIN0 already exists, skipping kinship"
    else
        echo "  [$panel] computing KING kinship (cutoff=$KING_CUTOFF)..."
        plink2 \
            --pfile "$BASE" \
            --make-king-table \
            --king-table-filter "$KING_CUTOFF" \
            --threads "$THREADS" \
            --memory "$MEM_MB" \
            --out "$BASE"

        [[ -f "$KIN0" ]] \
            || { echo "ERROR: $KIN0 not produced" >&2; exit 1; }
        N_PAIRS=$(( $(wc -l < "$KIN0") - 1 ))
        echo "  [$panel] $N_PAIRS sample pairs with kinship >= $KING_CUTOFF"
    fi
    echo
done

echo "=== Done ==="
echo "  GSA:        gsa/dense_kinship.{pgen,pvar,psam,kin0}"
echo "  Affymetrix: affymetrix/dense_kinship.{pgen,pvar,psam,kin0}"
