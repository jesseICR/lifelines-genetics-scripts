#!/usr/bin/env python3
"""
assign_eur_groups.py

For each sample in a plink2 --score output (.sscore), compute Euclidean
distance to each EUR group centroid in PC1..PC{N_PCS} space and write a
CSV listing the top-K closest groups (label + distance).

Usage:
    python3 assign_eur_groups.py CENTROIDS_TSV SSCORE OUT_CSV [TOP_K] [N_PCS]

CENTROIDS_TSV  TSV: group_label  center_PC1 ... center_PC{N_PCS}
SSCORE         plink2 --score output (.sscore) with PC1_AVG ... PCk_AVG
OUT_CSV        IID, group_1, distance_1, ..., group_K, distance_K  (1 + 2K cols)
TOP_K          number of nearest groups to report (default 10)
N_PCS          number of PCs to use (default 6)

Pure stdlib.
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
    top_k = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    n_pcs = int(sys.argv[5]) if len(sys.argv) > 5 else 6

    # --- centroids ---
    centroids = []
    with open(centroids_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            label = row["group_label"]
            coords = [float(row[f"center_PC{i}"]) for i in range(1, n_pcs + 1)]
            centroids.append((label, coords))

    if not centroids:
        sys.exit(f"ERROR: no centroids parsed from {centroids_path}")
    if top_k > len(centroids):
        sys.stderr.write(
            f"WARN: top_k ({top_k}) > number of centroids ({len(centroids)}); "
            f"reducing to {len(centroids)}\n"
        )
        top_k = len(centroids)

    # --- samples ---
    with open(sscore_path) as fin, open(out_path, "w", newline="") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        iid_idx = header.index("IID")
        pc_idxs = [header.index(f"PC{i}_AVG") for i in range(1, n_pcs + 1)]

        writer = csv.writer(fout)
        out_header = ["IID"]
        for k in range(1, top_k + 1):
            out_header.extend([f"group_{k}", f"distance_{k}"])
        writer.writerow(out_header)

        n_samples = 0
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            iid = fields[iid_idx]
            sample_pcs = [float(fields[i]) for i in pc_idxs]

            dists = []
            for label, c in centroids:
                d2 = sum((s - x) ** 2 for s, x in zip(sample_pcs, c))
                dists.append((d2, label))
            dists.sort()  # ascending by squared distance == ascending by distance

            row = [iid]
            for d2, label in dists[:top_k]:
                row.append(label)
                row.append(f"{math.sqrt(d2):.6f}")
            writer.writerow(row)
            n_samples += 1

    print(
        f"Assigned {n_samples} samples; top {top_k} of {len(centroids)} "
        f"centroids per sample (PC1..PC{n_pcs}). Output: {out_path}"
    )


if __name__ == "__main__":
    main()
