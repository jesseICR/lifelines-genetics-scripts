#!/bin/bash
# 08_compute_freq.sh
#
# Step 8. For each panel:
#   (a) plink2 --freq counts on the full rsid-keyed merged pfile (~7M
#       SNPs from snp.info), producing $panel/freq.acount.
#   (b) Compare each variant's A1 (=ALT) frequency in Lifelines against
#       A1Freq in snp.info and write two flagged CSVs per panel:
#         freq_diff_ge_0.1.csv  -- |A1_freq_lifelines - A1_freq_snpinfo| >= 0.1
#         freq_diff_ge_0.2.csv  -- |...| >= 0.2
#       Each CSV: rsid, a1_freq_lifelines, a1_freq_snpinfo, abs_diff
#
# Output is summary statistics only (no individual-level data) and may
# be exported off the cluster freely.
#
# Run via SLURM:
#     sbatch run_freq.sbatch
#
# Inputs (from 03_relabel_rsids.sh):
#   gsa/merged_rsid.{pgen,pvar,psam}
#   affymetrix/merged_rsid.{pgen,pvar,psam}
#   ./snp.info
#
# Outputs:
#   gsa/freq.acount               affymetrix/freq.acount
#   gsa/freq_diff_ge_0.1.csv      affymetrix/freq_diff_ge_0.1.csv
#   gsa/freq_diff_ge_0.2.csv      affymetrix/freq_diff_ge_0.2.csv

set -euo pipefail
trap 'echo "[08] ERROR on line $LINENO (last command: $BASH_COMMAND)" >&2' ERR

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
[[ -s "$SNPINFO" ]] || { echo "ERROR: $SNPINFO missing — run 01_download_snp_info.sh" >&2; exit 1; }
command -v plink2 >/dev/null 2>&1 \
    || { echo "ERROR: plink2 not on PATH after module load" >&2; exit 1; }
for panel in gsa affymetrix; do
    [[ -f "$panel/merged_rsid.pgen" ]] \
        || { echo "ERROR: $panel/merged_rsid.pgen missing — run 03_relabel_rsids.sh first" >&2; exit 1; }
done
echo "  resources:   ${THREADS} threads, ${MEM_MB} MB"
echo "  plink2:      $(plink2 --version 2>&1 | head -n 1)"
echo

compute_freq() {
    local panel=$1
    local out_acount="$panel/freq.acount"

    if [[ -f "$out_acount" ]]; then
        echo "  [$panel] $out_acount already exists, skipping plink2"
        return 0
    fi

    echo "  [$panel] computing allele counts..."
    plink2 \
        --pfile "$panel/merged_rsid" \
        --freq counts \
        --threads "$THREADS" \
        --memory "$MEM_MB" \
        --out "$panel/freq"

    [[ -f "$out_acount" ]] || { echo "ERROR: $out_acount not produced" >&2; exit 1; }
    local n=$(( $(wc -l < "$out_acount") - 1 ))
    echo "  [$panel] wrote $n variant rows to $out_acount"
}

# Single-pass awk: snp.info into a hash, then stream freq.acount and
# write rows to the two threshold CSVs in one go. ~30 s per panel.
compare_freq() {
    local panel=$1
    local acount="$panel/freq.acount"
    local out01="$panel/freq_diff_ge_0.1.csv"
    local out02="$panel/freq_diff_ge_0.2.csv"

    if [[ -f "$out01" && -f "$out02" ]]; then
        echo "  [$panel] both freq_diff CSVs already exist, skipping comparison"
        return 0
    fi

    echo "  [$panel] comparing A1 frequency to snp.info..."
    awk -v OUT01="$out01" -v OUT02="$out02" '
        BEGIN { OFS = "," }

        # First pass: snp.info -> hash by rsid (skip header)
        # Columns: 1=Chrom 2=ID 3=Index 4=GenPos 5=PhysPos 6=A1 7=A2 8=A1Freq
        NR == FNR {
            if (FNR == 1) next
            A1[$2] = $6
            A2[$2] = $7
            F[$2]  = $8 + 0
            next
        }

        # Second file: freq.acount header
        # Columns: #CHROM ID REF ALT PROVISIONAL_REF? ALT_CTS OBS_CT
        FNR == 1 {
            HDR = "rsid,a1_freq_lifelines,a1_freq_snpinfo,abs_diff"
            print HDR > OUT01
            print HDR > OUT02
            next
        }

        {
            n_total++
            rsid    = $2
            ref     = $3
            alt     = $4
            alt_cts = $6 + 0
            obs_ct  = $7 + 0

            if (obs_ct == 0)        { n_no_calls++;        next }
            if (!(rsid in A1))      { n_not_in_snpinfo++;  next }

            a1_si    = A1[rsid]
            a2_si    = A2[rsid]
            a1f_si   = F[rsid]
            a1f_ll   = alt_cts / obs_ct        # plink A1 == VCF ALT

            if      (a1_si == alt && a2_si == ref) { a1f_si_oriented = a1f_si       }
            else if (a1_si == ref && a2_si == alt) { a1f_si_oriented = 1 - a1f_si   }
            else                                   { n_allele_mismatch++; next     }

            n_matched++
            diff     = a1f_ll - a1f_si_oriented
            abs_diff = (diff < 0) ? -diff : diff

            if (abs_diff >= 0.1) {
                printf "%s,%.6f,%.6f,%.6f\n", rsid, a1f_ll, a1f_si_oriented, abs_diff > OUT01
                n_ge_01++
            }
            if (abs_diff >= 0.2) {
                printf "%s,%.6f,%.6f,%.6f\n", rsid, a1f_ll, a1f_si_oriented, abs_diff > OUT02
                n_ge_02++
            }
        }

        END {
            printf "    total variants in freq.acount: %d\n",   n_total            > "/dev/stderr"
            printf "    matched (alleles oriented ok): %d\n",   n_matched          > "/dev/stderr"
            printf "    diff >= 0.1:                   %d\n",   n_ge_01 + 0        > "/dev/stderr"
            printf "    diff >= 0.2:                   %d\n",   n_ge_02 + 0        > "/dev/stderr"
            printf "    allele mismatch (skipped):     %d\n",   n_allele_mismatch + 0 > "/dev/stderr"
            printf "    rsid not in snp.info:          %d\n",   n_not_in_snpinfo + 0  > "/dev/stderr"
            printf "    no-call variants (OBS_CT=0):   %d\n",   n_no_calls + 0     > "/dev/stderr"
        }
    ' "$SNPINFO" "$acount"

    echo "  [$panel] $(( $(wc -l < "$out01") - 1 )) rows in $out01"
    echo "  [$panel] $(( $(wc -l < "$out02") - 1 )) rows in $out02"
}

for panel in gsa affymetrix; do
    echo "=== $panel ==="
    compute_freq "$panel"
    compare_freq "$panel"
    echo
done

echo "=== Done ==="
echo "  gsa/freq.acount             affymetrix/freq.acount"
echo "  gsa/freq_diff_ge_0.1.csv    affymetrix/freq_diff_ge_0.1.csv"
echo "  gsa/freq_diff_ge_0.2.csv    affymetrix/freq_diff_ge_0.2.csv"
echo "  All output is summary-level (no individual data) and exportable."
