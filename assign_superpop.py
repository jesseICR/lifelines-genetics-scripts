#!/usr/bin/env python3
"""
assign_superpop.py

Assigns each sample in a plink2 --score output (.sscore) to its nearest
super-population centroid by Euclidean distance in PC space.

Usage:
    python3 assign_superpop.py CENTROIDS_TSV SSCORE OUT_CSV [N_PCS]

CENTROIDS_TSV  TSV with header: supervised_group n PC1_AVG ... PCk_AVG
SSCORE         plink2 --score output (.sscore) with PC1_AVG ... PCk_AVG cols
OUT_CSV        output CSV: FID,IID,supervised_group,euclidean_distance
N_PCS          number of PCs to use (default 20; must be <= cols available)

Pure stdlib. No numpy/pandas needed.
"""

import csv
import math
import sys


def main():
    if len(sys.argv) < 4:
        sys.stderr.write(__doc__)
        sys.exit(2)

    centroids_path = sys.argv[1]
    sscore_path = sys.argv[2]
    out_path = sys.argv[3]
    n_pcs = int(sys.argv[4]) if len(sys.argv) > 4 else 20

    # --- centroids ---------------------------------------------------
    centroids = []  # list of (name, [pc1, ..., pcN])
    with open(centroids_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row["supervised_group"]
            coords = [float(row[f"PC{i}_AVG"]) for i in range(1, n_pcs + 1)]
            centroids.append((name, coords))

    if not centroids:
        sys.exit(f"ERROR: no centroids parsed from {centroids_path}")

    # --- samples + assignment ---------------------------------------
    with open(sscore_path) as fin, open(out_path, "w", newline="") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        try:
            fid_idx = header.index("#FID")
        except ValueError:
            fid_idx = header.index("FID")
        iid_idx = header.index("IID")
        pc_idxs = [header.index(f"PC{i}_AVG") for i in range(1, n_pcs + 1)]

        writer = csv.writer(fout)
        writer.writerow(["FID", "IID", "supervised_group", "euclidean_distance"])

        n_samples = 0
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            sample_pcs = [float(fields[i]) for i in pc_idxs]

            best_name = None
            best_d2 = float("inf")
            for name, c in centroids:
                d2 = sum((s - x) ** 2 for s, x in zip(sample_pcs, c))
                if d2 < best_d2:
                    best_d2 = d2
                    best_name = name

            writer.writerow([
                fields[fid_idx],
                fields[iid_idx],
                best_name,
                f"{math.sqrt(best_d2):.6f}",
            ])
            n_samples += 1

    print(
        f"Assigned {n_samples} samples to nearest of {len(centroids)} centroids "
        f"using PC1..PC{n_pcs}. Output: {out_path}"
    )


if __name__ == "__main__":
    main()
