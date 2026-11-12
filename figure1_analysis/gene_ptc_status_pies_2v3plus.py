#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Load transcript-level table
# -----------------------------
iso = pd.read_csv(
    "AtRTD2_reannotated_CDS_transcript_summary.tsv",
    sep="\t"
)
iso.columns = iso.columns.str.strip()

# -----------------------------
# Define coding transcripts
# -----------------------------
iso["is_coding"] = (
    iso["cds_start_1based"].notna() &
    iso["cds_end_1based"].notna() &
    (iso["cds_start_1based"] > 0) &
    (iso["cds_end_1based"] > 0)
)

# -----------------------------
# Nuclear, coding transcripts only
# -----------------------------
coding_iso = iso[
    iso["is_coding"] &
    ~iso["gene_id"].str.contains("ChrM|ChrC", na=False)
].copy()

# -----------------------------
# Per-gene isoform / PTC counts
# -----------------------------
gene_status = (
    coding_iso
    .groupby("gene_id")
    .apply(lambda g: pd.Series({
        "Total_Isoforms": len(g),
        "PTC_Isoforms":   g["has_ptc"].sum(),
        "Non_PTC_Isoforms": (~g["has_ptc"]).sum()
    }), include_groups=False)
    .reset_index()
)

def classify(row):
    if row["PTC_Isoforms"] == 0:
        return "All No_PTC"
    if row["Non_PTC_Isoforms"] == 0:
        return "All PTC"
    return "Mixed"

gene_status["Status"] = gene_status.apply(classify, axis=1)

# -----------------------------
# Split: exactly 2 isoforms vs 3-or-more
# -----------------------------
two_iso   = gene_status[gene_status["Total_Isoforms"] == 2]
three_plus = gene_status[gene_status["Total_Isoforms"] >= 3]

print(f"Genes with 2 isoforms: {len(two_iso)}")
print(two_iso["Status"].value_counts())
print()
print(f"Genes with 3+ isoforms: {len(three_plus)}")
print(three_plus["Status"].value_counts())

# -----------------------------
# Plot
# -----------------------------
category_order  = ["All PTC", "All No_PTC", "Mixed"]
category_colors = {
    "All PTC":     "darkred",
    "All No_PTC":  "yellowgreen",
    "Mixed":       "darkorange",
}

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
max_n = max(len(two_iso), len(three_plus))

for ax, subset, title in zip(
    axes,
    [two_iso, three_plus],
    ["Genes with 2 isoforms", "Genes with 3+ isoforms"]
):
    counts = subset["Status"].value_counts().reindex(category_order, fill_value=0)
    counts = counts[counts > 0]
    colors = [category_colors[c] for c in counts.index]
    radius = np.sqrt(len(subset) / max_n) if max_n else 1

    ax.pie(
        counts,
        colors=colors,
        radius=radius,
        wedgeprops=dict(edgecolor="white", linewidth=1)
    )
    ax.set_title(title, fontsize=11)

handles = [
    plt.Rectangle((0, 0), 1, 1, color=category_colors[c])
    for c in category_order
]
axes[0].legend(
    handles, category_order,
    loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=9, frameon=True
)

plt.tight_layout()
plt.savefig("gene_ptc_status_pies_2v3plus.png", dpi=300, bbox_inches="tight")
plt.show()
