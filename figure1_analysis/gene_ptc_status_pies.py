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

# -----------------------------
# Classify each gene
# -----------------------------
def classify(row):
    if row["PTC_Isoforms"] == 0:
        return "All No_PTC"
    if row["Non_PTC_Isoforms"] == 0:
        return "All PTC"
    return "Mixed"

gene_status["Status"] = gene_status.apply(classify, axis=1)

# -----------------------------
# Split single-transcript vs alternatively-spliced genes
# -----------------------------
single_tx = gene_status[gene_status["Total_Isoforms"] == 1]
multi_tx  = gene_status[gene_status["Total_Isoforms"] >= 2]

print(f"Genes with 1 transcript: {len(single_tx)}")
print(single_tx["Status"].value_counts())
print()
print(f"Genes with alternative isoforms: {len(multi_tx)}")
print(multi_tx["Status"].value_counts())

# -----------------------------
# Plot: two pie charts, sized by relative gene count
# -----------------------------
category_order  = ["All PTC", "All No_PTC", "Mixed"]
category_colors = {
    "All PTC":     "darkred",
    "All No_PTC":  "yellowgreen",
    "Mixed":       "darkorange",
}

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

max_n = max(len(single_tx), len(multi_tx))

for ax, subset, title in zip(
    axes,
    [single_tx, multi_tx],
    ["Genes with 1 transcript", "Genes with alternative isoforms"]
):
    counts = subset["Status"].value_counts().reindex(category_order, fill_value=0)
    counts = counts[counts > 0]  # drop empty wedges (e.g. no Mixed in single-tx genes)
    colors = [category_colors[c] for c in counts.index]

    radius = np.sqrt(len(subset) / max_n) if max_n else 1

    ax.pie(
        counts,
        colors=colors,
        radius=radius,
        wedgeprops=dict(edgecolor="white", linewidth=1)
    )
    ax.set_title(title, fontsize=11)

# Single shared legend (all three categories, in fixed order/colors)
handles = [
    plt.Rectangle((0, 0), 1, 1, color=category_colors[c])
    for c in category_order
]
axes[0].legend(
    handles, category_order,
    loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=9, frameon=True
)

plt.tight_layout()
plt.savefig("gene_ptc_status_pies.png", dpi=300, bbox_inches="tight")
plt.show()
