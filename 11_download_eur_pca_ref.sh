#!/bin/bash
# 11_download_eur_pca_ref.sh
#
# Run on the LOGIN node. Downloads the EUR-balanced PCA reference into
# ./pca_ref_eur/. Idempotent.
#
# Files fetched:
#   eur_pca_balanced.eigenvec.allele   (PC loadings, used for --score)
#   eur_pca_balanced.acount            (allele counts, used by --read-freq)
#
# (centroids_pc6.tsv is shipped with this repo, not downloaded.)
#
# NOTE: as of writing, these two files must be present at:
#   https://raw.githubusercontent.com/jesseICR/public-statgen/main/outputs/pca/
# If they live elsewhere, edit BASE below.

set -euo pipefail
trap 'echo "[11] ERROR on line $LINENO" >&2' ERR

BASE=https://raw.githubusercontent.com/jesseICR/public-statgen/main/outputs/pca
DIR=pca_ref_eur
mkdir -p "$DIR"

FILES=(
    eur_pca_balanced.eigenvec.allele
    eur_pca_balanced.acount
)

for f in "${FILES[@]}"; do
    if [[ -s "$DIR/$f" ]]; then
        echo "$DIR/$f already present ($(wc -l < "$DIR/$f") lines). Skipping."
        continue
    fi
    curl -L --fail --retry 3 --retry-delay 5 -o "$DIR/$f" "$BASE/$f"
    echo "Downloaded $DIR/$f ($(wc -l < "$DIR/$f") lines)"
done
