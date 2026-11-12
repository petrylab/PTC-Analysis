#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import spearmanr
import numpy as np

# -----------------------------
# 1. Load transcript-level table
# -----------------------------
iso = pd.read_csv(
    "AtRTD2_reannotated_CDS_transcript_summary.tsv",
    sep="\t"
)
iso.columns = iso.columns.str.strip()

# -----------------------------
# 2. Define coding transcripts
# -----------------------------
iso["is_coding"] = (
    iso["cds_start_1based"].notna() &
    iso["cds_end_1based"].notna() &
    (iso["cds_start_1based"] > 0) &
    (iso["cds_end_1based"] > 0)
)

# -----------------------------
# 3. Keep nuclear, coding transcripts only
# -----------------------------
coding_iso = iso[
    iso["is_coding"] &
    ~iso["gene_id"].str.contains("ChrM|ChrC", na=False)
].copy()

# -----------------------------
# 4. Count isoforms per gene
# -----------------------------
gene_counts = (
    coding_iso
    .groupby("gene_id")
    .size()
    .reset_index(name="Total_Isoforms")
)

# -----------------------------
# 5. Merge counts back
# -----------------------------
coding_iso = coding_iso.merge(gene_counts, on="gene_id")

# -----------------------------
# 6. Keep genes with ≥2 isoforms
# -----------------------------
multi_isoform = coding_iso[coding_iso["Total_Isoforms"] >= 2].copy()

# -----------------------------
# 7. Bin genes (≥20 collapsed)
# -----------------------------
multi_isoform["Isoform_bin"] = multi_isoform["Total_Isoforms"].apply(
    lambda x: x if x < 20 else 20
)

# -----------------------------
# 8. Summary per bin
# -----------------------------
summary = (
    multi_isoform
    .groupby("Isoform_bin")
    .apply(lambda g: pd.Series({
        "Total_isoforms": len(g),
        "Non_PTC_isoforms": (~g["has_ptc"]).sum(),
        "Fraction_non_PTC": (~g["has_ptc"]).sum() / len(g)
    }), include_groups=False)
    .reset_index()
    .sort_values("Isoform_bin")
)

print(summary)

# -----------------------------
# 9. Spearman correlation (2–19)
# -----------------------------
filtered_summary = summary[summary["Isoform_bin"] < 20]
corr, pval = spearmanr(
    filtered_summary["Isoform_bin"],
    filtered_summary["Fraction_non_PTC"]
)

print(
    f"Spearman correlation (2–19 isoforms): "
    f"{corr:.3f}, p-value: {pval:.3e}"
)

# -----------------------------
# 10. Dual plots
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# -------- Left panel (numeric, continuous) --------
axes[0].plot(
    summary["Isoform_bin"],
    summary["Fraction_non_PTC"],
    marker="o",
    color="steelblue"
)
axes[0].set_xlabel("Number of isoforms per gene (20 = ≥20)")
axes[0].set_ylabel("Fraction of non-PTC isoforms")
axes[0].set_title("Non-PTC fraction vs gene isoform count")
axes[0].set_ylim(0, 1)
axes[0].grid(True)
axes[0].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

axes[0].text(
    0.05, 0.05,
    f"ρ = {corr:.2f} p = {pval:.1e}",
    transform=axes[0].transAxes,
    fontsize=10,
    bbox=dict(facecolor="white", alpha=0.8)
)

# -------- Right panel (categorical bins) --------
bin_labels = [
    str(b) if b < 20 else "≥20"
    for b in summary["Isoform_bin"]
]

x = np.arange(len(bin_labels))
bar_width = 0.4

axes[1].bar(
    x - bar_width / 2,
    summary["Total_isoforms"],
    width=bar_width,
    color="lightgray",
    label="Total isoforms"
)

axes[1].bar(
    x + bar_width / 2,
    summary["Non_PTC_isoforms"],
    width=bar_width,
    color="steelblue",
    label="Non-PTC isoforms"
)

axes[1].set_xticks(x)
axes[1].set_xticklabels(bin_labels)

axes[1].set_xlabel("Isoform-count bin per gene")
axes[1].set_ylabel("Isoform count")
axes[1].set_title("Total vs non-PTC isoforms per bin")
axes[1].legend()

plt.tight_layout()
plt.show()
