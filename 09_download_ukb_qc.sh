#!/bin/bash
# 09_download_ukb_qc.sh
#
# Run on the LOGIN node. Downloads the public UKB SNP QC file
# (ukb_snp_qc.txt) which lists per-SNP QC flags including
# `in_Relatedness` (the SNP set used for UKB kinship/relatedness).
# Idempotent.

set -euo pipefail
trap 'echo "[09] ERROR on line $LINENO" >&2' ERR

URL=https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt
OUT=ukb_snp_qc.txt

verify() {
    [[ -s "$OUT" ]] || return 1
    head -n 1 "$OUT" | grep -q 'in_Relatedness' || return 1
    return 0
}

if verify; then
    echo "$OUT already present and looks valid ($(wc -l < "$OUT") lines). Skipping."
    exit 0
fi

curl -L --fail --retry 3 --retry-delay 5 -o "$OUT" "$URL"

verify || { echo "ERROR: downloaded $OUT does not contain 'in_Relatedness' header" >&2; exit 1; }
echo "Downloaded $(wc -l < "$OUT") lines to ./$OUT"
