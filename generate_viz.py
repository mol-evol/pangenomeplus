#!/usr/bin/env python3
"""Generate visualizations for README from existing pipeline output."""

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Paths
output_dir = Path("readme_example_output")
viz_dir = output_dir / "visualizations"
viz_dir.mkdir(parents=True, exist_ok=True)

# Load data
with open(output_dir / "mappings/compact_id_mappings.json") as f:
    mappings_data = json.load(f)
    id_mappings = mappings_data.get("compact_to_full", {})

presence_absence = pd.read_csv(
    output_dir / "matrices/presence_absence_matrix.csv", index_col=0
)

# Count families by type (extract from family ID prefix in presence/absence matrix)
family_counts: dict[str, int] = defaultdict(int)
for family_id in presence_absence.columns:
    # Family IDs are like: FAM_P1, FAM_I5K, FAM_T42, FAM_R7, FAM_C12
    if family_id.startswith("FAM_"):
        prefix = family_id[4]  # Extract type: P, T, R, I, C
        feature_map = {
            "P": "Proteins",
            "T": "tRNAs",
            "R": "rRNAs",
            "I": "Intergenic",
            "C": "CRISPR",
        }
        feature_type = feature_map.get(prefix, "Unknown")
        family_counts[feature_type] += 1

# Filter out zero counts for cleaner visualization
family_counts_filtered: dict[str, int] = {
    k: v for k, v in family_counts.items() if v > 0
}

# Family classification counts
family_presence_counts = presence_absence.sum(axis=0)
n_genomes = len(presence_absence)
core_threshold = 0.95 * n_genomes
accessory_threshold = 0.15 * n_genomes

core_count = sum(family_presence_counts >= core_threshold)
accessory_count = sum(
    (family_presence_counts >= accessory_threshold)
    & (family_presence_counts < core_threshold)
)
cloud_count = sum(family_presence_counts < accessory_threshold)

print(f"Genomes: {n_genomes}")
print(f"Total families: {len(presence_absence.columns)}")
print(f"Core families: {core_count}")
print(f"Accessory families: {accessory_count}")
print(f"Cloud families: {cloud_count}")

# Figure 2: Feature Type Distribution
fig, ax = plt.subplots(figsize=(10, 6))
types = list(family_counts_filtered.keys())
counts = list(family_counts_filtered.values())
bars = ax.bar(types, counts, color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"])

for bar in bars:
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2.0,
        height,
        f"{int(height):,}",
        ha="center",
        va="bottom",
    )

ax.set_xlabel("Feature Type", fontsize=12)
ax.set_ylabel("Number of Families", fontsize=12)
ax.set_title("Gene Family Distribution by Feature Type", fontsize=14, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(viz_dir / "feature_type_distribution.png", dpi=300)
print(f"✓ Created {viz_dir}/feature_type_distribution.png")
plt.close()

# Figure 3: Family Classification Pie Chart
fig, ax = plt.subplots(figsize=(8, 6))
sizes = [core_count, accessory_count, cloud_count]
labels = ["Core", "Accessory", "Cloud"]
colors = ["#ff9999", "#66b3ff", "#99ff99"]
explode = (0.05, 0, 0)

ax.pie(
    sizes,
    explode=explode,
    labels=labels,
    colors=colors,
    autopct="%1.1f%%",
    shadow=True,
    startangle=90,
)
ax.set_title("Pangenome Family Classification", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(viz_dir / "family_classification.png", dpi=300)
print(f"✓ Created {viz_dir}/family_classification.png")
plt.close()

# Figure 4: Presence/Absence Heatmap (top 50 variable families)
family_variance = presence_absence.var(axis=0)
top_variable = family_variance.nlargest(50).index
heatmap_data = presence_absence[top_variable]

fig, ax = plt.subplots(figsize=(14, 8))
sns.heatmap(
    heatmap_data,
    cmap=["white", "#3182bd"],
    cbar_kws={"label": "Present"},
    ax=ax,
    xticklabels=[f[:15] + "..." if len(f) > 15 else f for f in top_variable],
    yticklabels=presence_absence.index,
)
ax.set_title(
    "Presence/Absence Pattern (Top 50 Variable Families)",
    fontsize=14,
    fontweight="bold",
)
ax.set_xlabel("Gene Families", fontsize=12)
ax.set_ylabel("Genomes", fontsize=12)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(viz_dir / "presence_absence_heatmap.png", dpi=300, bbox_inches="tight")
print(f"✓ Created {viz_dir}/presence_absence_heatmap.png")
plt.close()

# Figure 1: Rarefaction curves (simulated from data)
# Since we don't have rarefaction data, create a simple growth curve
genomes_range = range(1, n_genomes + 1)
cumulative_families = []
for i in genomes_range:
    subset = presence_absence.iloc[:i]
    n_families = (subset.sum(axis=0) > 0).sum()
    cumulative_families.append(n_families)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(
    genomes_range, cumulative_families, "b-", linewidth=2, marker="o", label="Pangenome"
)
ax.set_xlabel("Number of Genomes", fontsize=12)
ax.set_ylabel("Number of Gene Families", fontsize=12)
ax.set_title("Pangenome Growth Curve", fontsize=14, fontweight="bold")
ax.legend(loc="best", fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(viz_dir / "rarefaction_curves.png", dpi=300)
print(f"✓ Created {viz_dir}/rarefaction_curves.png")
plt.close()

print("\n✓ All visualizations generated successfully!")
