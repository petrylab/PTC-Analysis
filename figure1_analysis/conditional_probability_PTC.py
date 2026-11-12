#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import matplotlib.ticker as ticker

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
# Nuclear genes only
# -----------------------------
iso = iso[~iso["gene_id"].str.contains("ChrM|ChrC", na=False)]

# -----------------------------
# Keep only coding transcripts
# -----------------------------
coding_iso = iso[iso["is_coding"]].copy()

# -----------------------------
# Count PTC vs non-PTC isoforms per gene
# -----------------------------
counts = (
    coding_iso
    .groupby("gene_id")
    .apply(lambda g: pd.Series({
        "Non_PTC_Isoforms": (~g["has_ptc"]).sum(),
        "PTC_Isoforms": (g["has_ptc"]).sum()
    }))
    .reset_index()
)

# -----------------------------
# Limit dynamic range (≤ 6)
# -----------------------------
counts = counts[
    (counts["Non_PTC_Isoforms"] <= 6) &
    (counts["PTC_Isoforms"] <= 6)
]

# -----------------------------
# Conditional probability:
# Non-PTC → PTC
# -----------------------------
results_nptc2ptc = []

for n in range(1, 7):
    genes_with_n_or_more = (counts["Non_PTC_Isoforms"] >= n).sum()
    genes_with_n_and_ptc = (
        (counts["Non_PTC_Isoforms"] >= n) &
        (counts["PTC_Isoforms"] >= 1)
    ).sum()

    prob = (
        genes_with_n_and_ptc / genes_with_n_or_more
        if genes_with_n_or_more else None
    )

    results_nptc2ptc.append({
        "Isoforms": n,
        "Prob_next_PTC": prob,
        "Gene_count": genes_with_n_or_more
    })

summary_nptc2ptc = pd.DataFrame(results_nptc2ptc).dropna()

# -----------------------------
# Conditional probability:
# PTC → non-PTC
# -----------------------------
results_ptc2nptc = []

for n in range(1, 7):
    genes_with_n_or_more = (counts["PTC_Isoforms"] >= n).sum()
    genes_with_n_and_non_ptc = (
        (counts["PTC_Isoforms"] >= n) &
        (counts["Non_PTC_Isoforms"] >= 1)
    ).sum()

    prob = (
        genes_with_n_and_non_ptc / genes_with_n_or_more
        if genes_with_n_or_more else None
    )

    results_ptc2nptc.append({
        "Isoforms": n,
        "Prob_next_non_PTC": prob,
        "Gene_count": genes_with_n_or_more
    })

summary_ptc2nptc = pd.DataFrame(results_ptc2nptc).dropna()

# -----------------------------
# Spearman correlations
# -----------------------------
corr_1, pval_1 = spearmanr(
    summary_nptc2ptc["Isoforms"],
    summary_nptc2ptc["Prob_next_PTC"]
)

corr_2, pval_2 = spearmanr(
    summary_ptc2nptc["Isoforms"],
    summary_ptc2nptc["Prob_next_non_PTC"]
)

# -----------------------------
# Plot
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Non-PTC → PTC
axes[0].plot(
    summary_nptc2ptc["Isoforms"],
    summary_nptc2ptc["Prob_next_PTC"],
    marker="o",
    color="darkred"
)

for _, row in summary_nptc2ptc.iterrows():
    axes[0].text(
        row["Isoforms"],
        row["Prob_next_PTC"] + 0.02,
        str(row["Gene_count"]),
        ha="center",
        fontsize=8
    )

axes[0].set_xlabel("Non-PTC isoforms already present (max = 6)")
axes[0].set_ylabel("Probability next isoform is PTC-containing")
axes[0].set_title("Conditional probability: non-PTC → PTC")
axes[0].set_ylim(0, 1)
axes[0].grid(True)
axes[0].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
axes[0].text(
    0.05, 0.05,
    f"ρ = {corr_1:.2f}, p = {pval_1:.1e}",
    transform=axes[0].transAxes,
    bbox=dict(facecolor="white", alpha=0.7)
)

# PTC → non-PTC
axes[1].plot(
    summary_ptc2nptc["Isoforms"],
    summary_ptc2nptc["Prob_next_non_PTC"],
    marker="o",
    color="steelblue"
)

for _, row in summary_ptc2nptc.iterrows():
    axes[1].text(
        row["Isoforms"],
        row["Prob_next_non_PTC"] + 0.02,
        str(row["Gene_count"]),
        ha="center",
        fontsize=8
    )

axes[1].set_xlabel("PTC-containing isoforms already present (max = 6)")
axes[1].set_ylabel("Probability next isoform is non-PTC")
axes[1].set_title("Conditional probability: PTC → non-PTC")
axes[1].set_ylim(0, 1)
axes[1].grid(True)
axes[1].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
axes[1].text(
    0.05, 0.05,
    f"ρ = {corr_2:.2f}, p = {pval_2:.1e}",
    transform=axes[1].transAxes,
    bbox=dict(facecolor="white", alpha=0.7)
)

plt.tight_layout()
plt.show()
