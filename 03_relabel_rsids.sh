#!/bin/bash
# 03_relabel_rsids.sh
#
# Step 3 of 5. Renames variant IDs in the merged pfiles from
# chr:pos:ref:alt to the rsid given in snp.info.
#
# Run via SLURM:
#     sbatch run_relabel.sbatch
#
# Inputs (from 02_extract_plink.sh):
#   gsa/merged.{pgen,pvar,psam}
#   affymetrix/merged.{pgen,pvar,psam}
#
# Outputs:
#   gsa/merged_rsid.{pgen,pvar,psam}
#   affymetrix/merged_rsid.{pgen,pvar,psam}

set -euo pipefail
trap 'echo "[03] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

SNPINFO=snp.info
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
[[ -s "$SNPINFO" ]] || { echo "ERROR: $SNPINFO missing" >&2; exit 1; }
command -v plink2 >/dev/null 2>&1 \
    || { echo "ERROR: plink2 not on PATH after module load" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/merged.pgen" ]] || {
        echo "ERROR: $panel/merged.pgen missing — run 02_extract_plink.sh first" >&2
        exit 1
    }
done
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo

echo "=== Building rename map ==="
awk 'NR>1 {print $1":"$5":"$7":"$6"\t"$2}' "$SNPINFO" > rename_ids.txt
echo "  rename_ids.txt: $(wc -l < rename_ids.txt) entries"
echo

for panel in gsa affymetrix; do
    if [[ -f "$panel/merged_rsid.pgen" ]]; then
        echo "[$panel] merged_rsid.pgen already exists, skipping"
        continue
    fi

    echo "[$panel] relabeling variant IDs..."
    plink2 \
        --pfile "$panel/merged" \
        --update-name rename_ids.txt \
        --threads "$THREADS" \
        --memory "$MEM_MB" \
        --make-pgen \
        --out "$panel/merged_rsid"

    [[ -f "$panel/merged_rsid.pgen" ]] \
        || { echo "ERROR: $panel/merged_rsid.pgen not produced" >&2; exit 1; }
    echo "  [$panel] $(wc -l < "$panel/merged_rsid.pvar") variant lines (incl. header)"
done

echo
echo "=== Done ==="
echo "  GSA:        gsa/merged_rsid.{pgen,pvar,psam}"
echo "  Affymetrix: affymetrix/merged_rsid.{pgen,pvar,psam}"
echo "  Next step:  sbatch run_dense.sbatch"
