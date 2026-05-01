#!/bin/bash
# 05_extract_dense.sh
#
# Step 5 of 5. Extract the dense rsid subset from the rsid-keyed merged
# pfiles produced by step 03 (relabel).
#
# Run via SLURM:
#     sbatch run_dense.sbatch
#
# Inputs:
#   ./rsids_dense_chr1_22.txt              (from 04_download_dense_rsids.sh)
#   gsa/merged_rsid.{pgen,pvar,psam}       (from 03_relabel_rsids.sh)
#   affymetrix/merged_rsid.{pgen,pvar,psam}
#
# Outputs:
#   gsa/dense.{pgen,pvar,psam}
#   affymetrix/dense.{pgen,pvar,psam}

set -euo pipefail
trap 'echo "[05] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

DENSE_LIST=rsids_dense_chr1_22.txt
THREADS=${SLURM_CPUS_PER_TASK:-2}
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && (( SLURM_MEM_PER_NODE > 2048 )); then
    MEM_MB=$(( SLURM_MEM_PER_NODE - 1024 ))
else
    MEM_MB=${PLINK_MEM_MB:-8000}
fi

module load PLINK/2.0-alpha6.20-20250707

echo "=== Pre-flight ==="
echo "  hostname:    $(hostname)"
echo "  job:         ${SLURM_JOB_ID:-<not under SLURM>}"
echo "  cwd:         $(pwd)"
[[ -s "$DENSE_LIST" ]] || {
    echo "ERROR: $DENSE_LIST missing — run 04_download_dense_rsids.sh first" >&2
    exit 1
}
command -v plink2 >/dev/null 2>&1 \
    || { echo "ERROR: plink2 not on PATH after module load" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/merged_rsid.pgen" ]] || {
        echo "ERROR: $panel/merged_rsid.pgen missing — run 03_relabel_rsids.sh first" >&2
        exit 1
    }
done
DENSE_N=$(wc -l < "$DENSE_LIST")
echo "  dense list:  $DENSE_N rsids"
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo

for panel in gsa affymetrix; do
    if [[ -f "$panel/dense.pgen" ]]; then
        echo "[$panel] dense.pgen already exists, skipping"
        continue
    fi

    echo "[$panel] extracting dense subset..."
    plink2 \
        --pfile "$panel/merged_rsid" \
        --extract "$DENSE_LIST" \
        --threads "$THREADS" \
        --memory "$MEM_MB" \
        --make-pgen \
        --out "$panel/dense"

    [[ -f "$panel/dense.pgen" ]] \
        || { echo "ERROR: $panel/dense.pgen not produced" >&2; exit 1; }

    KEPT=$(( $(wc -l < "$panel/dense.pvar") - $(grep -c '^#' "$panel/dense.pvar") ))
    echo "  [$panel] kept $KEPT / $DENSE_N dense rsids"
    if (( KEPT < DENSE_N )); then
        echo "  [$panel] note: $((DENSE_N - KEPT)) dense rsids absent from merged_rsid"
        echo "             (most likely lost during step 02 position match or step 02 merge)"
    fi
done

echo
echo "=== Done ==="
echo "  GSA:        gsa/dense.{pgen,pvar,psam}"
echo "  Affymetrix: affymetrix/dense.{pgen,pvar,psam}"
