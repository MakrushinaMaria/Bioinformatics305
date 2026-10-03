#!/usr/bin/env python3
"""Analyze exported soil QIIME2 tables and build report figures/summaries."""
from __future__ import annotations

import json
import math
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

COLORS = {"k0": "#4C78A8", "k5": "#F2A541", "k25": "#D64F4F"}
TIME_ORDER = ["d3", "d90", "d180", "d360"]
COND_ORDER = ["k0", "k5", "k25"]


def rank_name(taxon: str, prefix: str) -> str:
    for item in str(taxon).split(";"):
        item = item.strip()
        if item.startswith(prefix):
            value = item[len(prefix):].strip()
            return value if value else "Unclassified"
    return "Unclassified"


def load_pcoa(path: Path):
    lines = path.read_text().splitlines()
    eigvals = np.array([float(x) for x in lines[1].split("\t")])
    prop_idx = lines.index(next(x for x in lines if x.startswith("Proportion explained")))
    props = np.array([float(x) for x in lines[prop_idx + 1].split("\t")])
    site_idx = next(i for i, x in enumerate(lines) if x.startswith("Site\t"))
    n_sites, n_axes = (int(x) for x in lines[site_idx].split("\t")[1:3])
    rows = []
    for line in lines[site_idx + 1: site_idx + 1 + n_sites]:
        parts = line.split("\t")
        rows.append([parts[0], *map(float, parts[1:1 + n_axes])])
    cols = ["sample-id"] + [f"Axis{i + 1}" for i in range(n_axes)]
    return pd.DataFrame(rows, columns=cols), eigvals, props


metadata = pd.read_csv(TABLES / "soil_metadata_full.tsv", sep="\t", comment="#")
dada = pd.read_csv(TABLES / "dada2_stats.tsv", sep="\t", comment="#")
dada = dada.rename(columns={dada.columns[0]: "sample-id"})
shannon = pd.read_csv(TABLES / "shannon.tsv", sep="\t")
shannon = shannon.rename(columns={shannon.columns[0]: "sample-id"})

feature = pd.read_csv(TABLES / "feature-table.tsv", sep="\t", skiprows=1)
feature = feature.rename(columns={feature.columns[0]: "Feature ID"}).set_index("Feature ID")
taxonomy = pd.read_csv(TABLES / "taxonomy.tsv", sep="\t").set_index("Feature ID")
taxonomy["Class"] = taxonomy["Taxon"].map(lambda x: rank_name(x, "c__"))
taxonomy["Genus"] = taxonomy["Taxon"].map(lambda x: rank_name(x, "g__"))

input_total = int(dada["input"].sum())
filtered_total = int(dada["filtered"].sum())
denoised_total = int(dada["denoised"].sum())
nonchim_total = int(dada["non-chimeric"].sum())
dada_summary = {
    "input_total": input_total,
    "filtered_total": filtered_total,
    "denoised_total": denoised_total,
    "nonchimeric_total": nonchim_total,
    "filtered_percent": 100 * filtered_total / input_total,
    "nonchimeric_percent": 100 * nonchim_total / input_total,
    "smallest_nonchimeric": int(dada["non-chimeric"].min()),
    "smallest_sample": dada.loc[dada["non-chimeric"].idxmin(), "sample-id"],
    "largest_nonchimeric": int(dada["non-chimeric"].max()),
    "largest_sample": dada.loc[dada["non-chimeric"].idxmax(), "sample-id"],
}
pd.DataFrame([dada_summary]).to_csv(TABLES / "dada2_summary.csv", index=False)

# FastQC figure extracted verbatim from the selected sample report.
fastqc_zip = FIGURES / "SRR17307316_1_fastqc.zip"
with zipfile.ZipFile(fastqc_zip) as archive:
    member = next(x for x in archive.namelist() if x.endswith("/Images/per_base_quality.png"))
    (FIGURES / "task1_fastqc_per_base_quality.png").write_bytes(archive.read(member))

# DADA2 aggregate retention figure.
stages = ["Input", "Filtered", "Denoised", "Non-chimeric"]
totals = [input_total, filtered_total, denoised_total, nonchim_total]
fig, ax = plt.subplots(figsize=(7.2, 4.2))
bars = ax.bar(stages, totals, color=["#808080", "#5B8FF9", "#61DDAA", "#65789B"])
ax.set_ylabel("Reads across 36 samples")
ax.set_title("DADA2 read retention")
ax.spines[["top", "right"]].set_visible(False)
for bar, value in zip(bars, totals):
    ax.text(bar.get_x() + bar.get_width()/2, value + max(totals)*0.012,
            f"{value:,}\n({100*value/input_total:.1f}%)", ha="center", va="bottom", fontsize=9)
ax.set_ylim(0, max(totals)*1.15)
fig.tight_layout()
fig.savefig(FIGURES / "task3_dada2_retention.png", dpi=220)
plt.close(fig)


def rank_relative(rank: str):
    joined = feature.join(taxonomy[[rank]], how="left")
    grouped = joined.groupby(rank).sum(numeric_only=True).T
    rel = grouped.div(grouped.sum(axis=1), axis=0)
    rel.index.name = "sample-id"
    long = rel.reset_index().melt("sample-id", var_name=rank, value_name="relative_abundance")
    long = long.merge(metadata, on="sample-id")
    means = (long.groupby(["contamination_level", "day_after_contamination", rank], as_index=False)
                 ["relative_abundance"].mean())
    return rel, means


class_rel, class_means = rank_relative("Class")
genus_rel, genus_means = rank_relative("Genus")
class_means.to_csv(TABLES / "class_relative_abundance_group_means.csv", index=False)
genus_means.to_csv(TABLES / "genus_relative_abundance_group_means.csv", index=False)


def stacked_rank_figure(means: pd.DataFrame, rank: str, n_top: int, filename: str, title: str):
    top = (means.groupby(rank)["relative_abundance"].mean().nlargest(n_top).index.tolist())
    data = means.copy()
    data[rank] = data[rank].where(data[rank].isin(top), "Other")
    data = data.groupby(["contamination_level", "day_after_contamination", rank], as_index=False)["relative_abundance"].sum()
    data["group"] = data["contamination_level"] + "_" + data["day_after_contamination"]
    order = [f"{c}_{t}" for c in COND_ORDER for t in TIME_ORDER]
    pivot = data.pivot(index="group", columns=rank, values="relative_abundance").fillna(0).reindex(order)
    cols = [x for x in top if x in pivot.columns] + (["Other"] if "Other" in pivot.columns else [])
    pivot = pivot[cols]
    fig, ax = plt.subplots(figsize=(11.2, 5.2))
    palette = list(plt.get_cmap("tab20").colors)
    bottom = np.zeros(len(pivot))
    for i, col in enumerate(pivot.columns):
        vals = pivot[col].to_numpy() * 100
        ax.bar(np.arange(len(pivot)), vals, bottom=bottom, label=col, color=palette[i % len(palette)], width=0.82)
        bottom += vals
    ax.set_xticks(np.arange(len(pivot)), [x.replace("_d", "\nd") for x in pivot.index])
    ax.set_ylabel("Mean relative abundance, %")
    ax.set_title(title)
    ax.set_ylim(0, 100)
    ax.legend(title=rank, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, title_fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / filename, dpi=220, bbox_inches="tight")
    plt.close(fig)


stacked_rank_figure(class_means, "Class", 8, "task5_taxonomy_level3.png",
                    "Taxonomic composition at Level 3 (class)")
stacked_rank_figure(genus_means, "Genus", 10, "task5_taxonomy_level6.png",
                    "Taxonomic composition at Level 6 (genus)")

# Ranked taxonomic facts for the report.
tax_facts = {}
for time in ["d3", "d360"]:
    block = class_means[class_means.day_after_contamination == time]
    for cond in COND_ORDER:
        top = block[block.contamination_level == cond].nlargest(5, "relative_abundance")
        tax_facts[f"classes_{cond}_{time}"] = [
            {"taxon": r.Class, "percent": 100*r.relative_abundance} for r in top.itertuples()
        ]
contam = genus_means[genus_means.contamination_level.isin(["k5", "k25"])]
top_contam = contam.groupby("Genus")["relative_abundance"].mean().sort_values(ascending=False)
tax_facts["top_contaminated_genera"] = [
    {"taxon": idx, "percent": 100*val} for idx, val in top_contam.head(10).items()
]

# Weighted UniFrac PCoA: coordinates, centroids and recovery distances.
pcoa, eigvals, props = load_pcoa(TABLES / "weighted_unifrac_pcoa.txt")
pcoa = pcoa.merge(metadata, on="sample-id")
centroids = (pcoa.groupby(["contamination_level", "day_after_contamination"], as_index=False)
                 [["Axis1", "Axis2", "Axis3"]].mean())
centroids.to_csv(TABLES / "weighted_unifrac_centroids.csv", index=False)
dist_rows = []
for time in TIME_ORDER:
    ref = centroids[(centroids.contamination_level == "k0") & (centroids.day_after_contamination == time)].iloc[0]
    for cond in COND_ORDER:
        row = centroids[(centroids.contamination_level == cond) & (centroids.day_after_contamination == time)].iloc[0]
        distance = math.sqrt(sum((row[f"Axis{i}"] - ref[f"Axis{i}"])**2 for i in (1, 2, 3)))
        dist_rows.append({"contamination_level": cond, "day_after_contamination": time,
                          "distance_to_time_matched_k0_centroid_axes1_3": distance})
distances = pd.DataFrame(dist_rows)
distances.to_csv(TABLES / "weighted_unifrac_recovery_distances.csv", index=False)

fig, ax = plt.subplots(figsize=(7.6, 6.0))
markers = {"d3": "o", "d90": "s", "d180": "^", "d360": "D"}
for cond in COND_ORDER:
    for time in TIME_ORDER:
        block = pcoa[(pcoa.contamination_level == cond) & (pcoa.day_after_contamination == time)]
        ax.scatter(block.Axis1, block.Axis2, s=38, alpha=0.55, color=COLORS[cond], marker=markers[time])
    path = centroids[centroids.contamination_level == cond].set_index("day_after_contamination").loc[TIME_ORDER]
    ax.plot(path.Axis1, path.Axis2, color=COLORS[cond], linewidth=2.2, marker="o", label=cond)
    for time, row in path.iterrows():
        ax.annotate(time, (row.Axis1, row.Axis2), xytext=(4, 4), textcoords="offset points", fontsize=8)
ax.axhline(0, color="#dddddd", linewidth=0.8); ax.axvline(0, color="#dddddd", linewidth=0.8)
ax.set_xlabel(f"Axis 1 ({props[0]*100:.1f}%)")
ax.set_ylabel(f"Axis 2 ({props[1]*100:.1f}%)")
ax.set_title("Weighted UniFrac PCoA: samples and group trajectories")
ax.legend(title="Kerosene", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES / "task6_weighted_unifrac_pcoa.png", dpi=220)
plt.close(fig)

# Shannon summaries and figure.
shannon = shannon.merge(metadata, on="sample-id")
shannon_summary = (shannon.groupby(["contamination_level", "day_after_contamination"])
                   ["shannon_entropy"].agg(["mean", "std", "min", "max"]).reset_index())
shannon_summary.to_csv(TABLES / "shannon_group_summary.csv", index=False)
fig, ax = plt.subplots(figsize=(7.8, 4.8))
x = np.arange(len(TIME_ORDER))
for cond in COND_ORDER:
    block = shannon_summary[shannon_summary.contamination_level == cond].set_index("day_after_contamination").loc[TIME_ORDER]
    ax.errorbar(x, block["mean"], yerr=block["std"], color=COLORS[cond], marker="o",
                capsize=3, linewidth=2.2, label=cond)
ax.set_xticks(x, ["3", "90", "180", "360"])
ax.set_xlabel("Days after contamination")
ax.set_ylabel("Shannon entropy (mean ± SD, n=3)")
ax.set_title("Alpha diversity over time")
ax.legend(title="Kerosene", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES / "task7_shannon.png", dpi=220)
plt.close(fig)

summary = {
    "dada2": dada_summary,
    "taxonomy": tax_facts,
    "pcoa_proportion_explained_first3": [float(x) for x in props[:3]],
    "weighted_unifrac_recovery_distances": dist_rows,
    "shannon": shannon_summary.to_dict(orient="records"),
}
(TABLES / "analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
print(json.dumps(summary, indent=2, ensure_ascii=False))
