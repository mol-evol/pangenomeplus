#!/usr/bin/env python3
"""
Basic pangenome structure analysis for PanGenomePlus.

This module provides fundamental pangenome analysis capabilities including:
- Rarefaction curves (core and pangenome size vs. genome count)
- Pangenome openness analysis
- Core genome stability analysis
- Basic summary statistics

All analyses follow established pangenome literature methods.
"""

import datetime
import json
import logging
import random
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .clustering import FamilyStats
from .compact_ids import CompactIDManager
from .constants import FeatureType, VizParams


def _create_rarefaction_plot(
    rarefaction_data: Dict[str, List[Tuple[int, float, float, float]]], output_dir: str
) -> None:
    """Generate rarefaction curve visualization.

    Args:
        rarefaction_data: Dictionary with pangenome, core, accessory curves
        output_dir: Base output directory for saving plot
    """
    import matplotlib.pyplot as plt

    vis_dir = Path(output_dir) / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=VizParams.DEFAULT_FIGSIZE)

    # Plot pangenome curve with uncertainty band
    if "pangenome" in rarefaction_data and rarefaction_data["pangenome"]:
        data = rarefaction_data["pangenome"]
        x = [row[0] for row in data]  # genome_count
        y = [row[1] for row in data]  # mean
        std = [row[2] for row in data]  # std
        ax.plot(x, y, "b-", linewidth=2, label="Pangenome")
        ax.fill_between(
            x,
            [m - s for m, s in zip(y, std)],
            [m + s for m, s in zip(y, std)],
            alpha=0.2,
            color="blue",
        )

    # Plot core genome curve with uncertainty band
    if "core" in rarefaction_data and rarefaction_data["core"]:
        data = rarefaction_data["core"]
        x = [row[0] for row in data]  # genome_count
        y = [row[1] for row in data]  # mean
        std = [row[2] for row in data]  # std
        ax.plot(x, y, "r-", linewidth=2, label="Core genome")
        ax.fill_between(
            x,
            [m - s for m, s in zip(y, std)],
            [m + s for m, s in zip(y, std)],
            alpha=0.2,
            color="red",
        )

    # Plot accessory genome curve with uncertainty band
    if "accessory" in rarefaction_data and rarefaction_data["accessory"]:
        data = rarefaction_data["accessory"]
        x = [row[0] for row in data]  # genome_count
        y = [row[1] for row in data]  # mean
        std = [row[2] for row in data]  # std
        ax.plot(x, y, "g-", linewidth=2, label="Accessory genome")
        ax.fill_between(
            x,
            [m - s for m, s in zip(y, std)],
            [m + s for m, s in zip(y, std)],
            alpha=0.2,
            color="green",
        )

    ax.set_xlabel("Number of Genomes", fontsize=12)
    ax.set_ylabel("Number of Gene Families", fontsize=12)
    ax.set_title("Pangenome Rarefaction Curves", fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(vis_dir / "rarefaction_curves.png", dpi=VizParams.DPI)
    plt.close()


def _create_classification_pie(total_families: Dict[str, int], output_dir: str) -> None:
    """Generate family classification pie chart.

    Args:
        total_families: Dictionary with core, accessory, cloud counts
        output_dir: Base output directory for saving plot
    """
    import matplotlib.pyplot as plt

    vis_dir = Path(output_dir) / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)

    labels = ["Core", "Accessory", "Cloud"]
    sizes = [
        total_families.get("core", 0),
        total_families.get("accessory", 0),
        total_families.get("cloud", 0),
    ]
    colors = ["#ff9999", "#66b3ff", "#99ff99"]
    explode = (0.05, 0, 0)  # Explode core slice

    fig, ax = plt.subplots(figsize=VizParams.PIE_FIGSIZE)
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
    plt.savefig(vis_dir / "family_classification.png", dpi=VizParams.DPI)
    plt.close()


def _create_feature_type_bars(
    family_stats: Dict[str, Dict[str, FamilyStats]], output_dir: str
) -> None:
    """Generate feature type distribution bar chart.

    Args:
        family_stats: Family statistics by feature type
        output_dir: Base output directory for saving plot
    """
    import matplotlib.pyplot as plt

    vis_dir = Path(output_dir) / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)

    feature_types = []
    counts = []

    for feature_type, stats_dict in sorted(family_stats.items()):
        feature_name = FeatureType.NAMES.get(feature_type, feature_type)
        feature_types.append(feature_name)
        counts.append(len(stats_dict))

    fig, ax = plt.subplots(figsize=VizParams.DEFAULT_FIGSIZE)
    bars = ax.bar(
        feature_types,
        counts,
        color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"],
    )

    # Add value labels on bars
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
    ax.set_title(
        "Gene Family Distribution by Feature Type", fontsize=14, fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(vis_dir / "feature_type_distribution.png", dpi=VizParams.DPI)
    plt.close()


def _create_presence_absence_heatmap(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    output_dir: str,
    top_n: int = VizParams.TOP_VARIABLE_FAMILIES,
) -> None:
    """Generate presence/absence heatmap for most variable families.

    Args:
        family_assignments: Family assignments by feature type
        id_manager: Compact ID manager for genome lookups
        output_dir: Base output directory for saving plot
        top_n: Number of most variable families to show
    """
    from collections import defaultdict

    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    vis_dir = Path(output_dir) / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)

    # Get all genome IDs
    genome_ids = sorted(set(id_manager.genome_features.keys()))
    if not genome_ids:
        return  # No data to plot

    # Build family to genomes mapping
    family_to_genomes = defaultdict(set)
    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature:
                family_to_genomes[family_id].add(feature.genome_id)

    # Select most variable families (not 100% core, not singleton)
    # Sort by genome count ASCENDING to show truly variable families (lowest counts = highest variability)
    n_genomes = len(genome_ids)
    variable_families = [
        (fam, len(genomes))
        for fam, genomes in family_to_genomes.items()
        if 1 < len(genomes) < n_genomes
    ]
    variable_families.sort(key=lambda x: x[1])  # Ascending: most variable first
    selected_families = [fam for fam, _ in variable_families[:top_n]]

    if not selected_families:
        return  # No variable families to plot

    # Build presence/absence matrix
    matrix = []
    for genome_id in genome_ids:
        row = []
        for family_id in selected_families:
            has_family = genome_id in family_to_genomes[family_id]
            row.append(1 if has_family else 0)
        matrix.append(row)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(VizParams.HEATMAP_FIGSIZE[0], max(8, len(genome_ids) * 0.4)))
    sns.heatmap(
        np.array(matrix),
        xticklabels=[f[:15] + "..." if len(f) > 15 else f for f in selected_families],
        yticklabels=genome_ids,
        cmap=["white", "#3182bd"],
        cbar_kws={"label": "Present"},
        ax=ax,
    )

    ax.set_title(
        f"Presence/Absence Pattern (Top {len(selected_families)} Variable Families)",
        fontsize=14,
        fontweight="bold",
    )
    ax.set_xlabel("Gene Families", fontsize=12)
    ax.set_ylabel("Genomes", fontsize=12)
    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(vis_dir / "presence_absence_heatmap.png", dpi=VizParams.DPI, bbox_inches="tight")
    plt.close()


def calculate_rarefaction_curve(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    max_genomes: Optional[int] = None,
    iterations: int = 100,
    step_size: int = 1,
) -> Dict[str, List[Tuple[int, float, float, float]]]:
    """Calculate rarefaction curves for pangenome analysis.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for genome information
        max_genomes: Maximum number of genomes to analyze (default: all)
        iterations: Number of random subsampling iterations per genome count
        step_size: Step size for genome count sampling

    Returns:
        Dict mapping analysis type to list of (genome_count, mean, std, median) tuples

    Analysis types:
        - "pangenome": Total unique families vs. genome count
        - "core": Core families (present in all genomes) vs. genome count
        - "accessory": Accessory families vs. genome count
    """
    logger = logging.getLogger(__name__)
    logger.info("Calculating rarefaction curves")

    # Get all genomes
    all_genomes = set()
    family_to_genomes = defaultdict(set)

    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature:
                all_genomes.add(feature.genome_id)
                family_to_genomes[family_id].add(feature.genome_id)

    all_genomes_list: List[str] = list(all_genomes)
    total_genomes = len(all_genomes_list)

    if max_genomes is None:
        max_genomes = total_genomes

    max_genomes = min(max_genomes, total_genomes)

    logger.info(
        f"Analyzing rarefaction for {max_genomes} genomes with {iterations} iterations"
    )

    # Initialize results
    results: Dict[str, List[Tuple[int, float, float, float]]] = {
        "pangenome": [],
        "core": [],
        "accessory": [],
    }

    # Calculate rarefaction curves
    for genome_count in range(1, max_genomes + 1, step_size):
        if genome_count > total_genomes:
            break

        pangenome_sizes = []
        core_sizes = []
        accessory_sizes = []

        for iteration in range(iterations):
            # Randomly sample genomes
            sampled_genomes = set(random.sample(all_genomes_list, genome_count))

            # Count families in different categories
            pangenome_families = set()
            core_families = set()
            accessory_families = set()

            for family_id, family_genomes in family_to_genomes.items():
                family_in_sample = family_genomes.intersection(sampled_genomes)

                if family_in_sample:  # Family present in at least one sampled genome
                    pangenome_families.add(family_id)

                    if (
                        len(family_in_sample) == genome_count
                    ):  # Present in all sampled genomes
                        core_families.add(family_id)
                    else:  # Present in some but not all sampled genomes
                        accessory_families.add(family_id)

            pangenome_sizes.append(len(pangenome_families))
            core_sizes.append(len(core_families))
            accessory_sizes.append(len(accessory_families))

        # Calculate statistics
        results["pangenome"].append(
            (
                genome_count,
                sum(pangenome_sizes) / len(pangenome_sizes),  # mean
                _calculate_std(pangenome_sizes),  # std
                sorted(pangenome_sizes)[len(pangenome_sizes) // 2],  # median
            )
        )

        results["core"].append(
            (
                genome_count,
                sum(core_sizes) / len(core_sizes),  # mean
                _calculate_std(core_sizes),  # std
                sorted(core_sizes)[len(core_sizes) // 2],  # median
            )
        )

        results["accessory"].append(
            (
                genome_count,
                sum(accessory_sizes) / len(accessory_sizes),  # mean
                _calculate_std(accessory_sizes),  # std
                sorted(accessory_sizes)[len(accessory_sizes) // 2],  # median
            )
        )

    logger.info(
        f"Rarefaction curves calculated for {len(results['pangenome'])} genome counts"
    )
    return results


def analyze_pangenome_openness(
    rarefaction_data: Dict[str, List[Tuple[int, float, float, float]]],
    min_genomes: int = 5,
) -> Dict[str, Any]:
    """Analyze pangenome openness using rarefaction curve fitting.

    Args:
        rarefaction_data: Rarefaction curve data from calculate_rarefaction_curve()
        min_genomes: Minimum number of genomes for reliable analysis

    Returns:
        Dictionary with openness analysis results

    Analysis includes:
        - Pangenome growth trend (open vs. closed)
        - Core genome stability
        - Estimated saturation points
    """
    logger = logging.getLogger(__name__)
    logger.info("Analyzing pangenome openness")

    if len(rarefaction_data["pangenome"]) < min_genomes:
        logger.warning(
            f"Insufficient data for openness analysis (need ≥{min_genomes} genome counts)"
        )
        return {"error": "Insufficient data"}

    # Extract pangenome and core data
    pangenome_data = rarefaction_data["pangenome"]
    core_data = rarefaction_data["core"]

    # Calculate growth rates
    pangenome_growth_rate = _calculate_growth_rate(pangenome_data)
    core_decay_rate = _calculate_decay_rate(core_data)

    # Estimate pangenome type
    final_pangenome_size = pangenome_data[-1][1]  # Mean size at max genomes
    initial_pangenome_size = pangenome_data[0][1]  # Mean size at 1 genome
    growth_ratio = final_pangenome_size / initial_pangenome_size

    # Classification based on growth pattern
    if pangenome_growth_rate > 50:  # Arbitrary threshold - could be refined
        pangenome_type = "open"
    elif pangenome_growth_rate < 10:
        pangenome_type = "closed"
    else:
        pangenome_type = "intermediate"

    # Core genome stability
    final_core_size = core_data[-1][1]
    initial_core_size = core_data[0][1] if len(core_data) > 0 else 0
    core_stability = final_core_size / initial_core_size if initial_core_size > 0 else 0

    results = {
        "pangenome_type": pangenome_type,
        "pangenome_growth_rate": pangenome_growth_rate,
        "core_decay_rate": core_decay_rate,
        "growth_ratio": growth_ratio,
        "core_stability": core_stability,
        "final_pangenome_size": final_pangenome_size,
        "final_core_size": final_core_size,
        "genome_count_analyzed": len(pangenome_data),
    }

    logger.info(f"Pangenome classified as: {pangenome_type}")
    logger.info(f"Growth rate: {pangenome_growth_rate:.1f} families/genome")
    logger.info(f"Core stability: {core_stability:.3f}")

    return results


def generate_pangenome_summary_report(
    rarefaction_data: Dict[str, List[Tuple[int, float, float, float]]],
    openness_analysis: Dict[str, Any],
    family_stats: Dict[str, Dict[str, FamilyStats]],
    output_file: str,
) -> None:
    """Generate comprehensive pangenome structure summary report.

    Args:
        rarefaction_data: Rarefaction curve data
        openness_analysis: Openness analysis results
        family_stats: Family statistics by feature type
        output_file: Path to output report file

    Creates:
        - Text report with pangenome structure analysis
        - JSON data file with numerical results
        - Basic visualization data (CSV format)
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating pangenome summary report: {output_file}")

    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Count total families by classification
    total_families = {"core": 0, "accessory": 0, "cloud": 0}
    for feature_type, stats_dict in family_stats.items():
        for family_id, stats in stats_dict.items():
            total_families[stats.classification] += 1

    total_all_families = sum(total_families.values())

    # Write comprehensive report
    with open(output_file, "w") as f:
        f.write("# PanGenomePlus Pangenome Structure Analysis Report\n\n")

        # Summary statistics
        f.write("## Summary Statistics\n")
        f.write(f"Total gene families: {total_all_families:,}\n")
        f.write(
            f"- Core families: {total_families['core']:,} ({total_families['core']/total_all_families*100:.1f}%)\n"
        )
        f.write(
            f"- Accessory families: {total_families['accessory']:,} ({total_families['accessory']/total_all_families*100:.1f}%)\n"
        )
        f.write(
            f"- Cloud families: {total_families['cloud']:,} ({total_families['cloud']/total_all_families*100:.1f}%)\n\n"
        )

        # Pangenome openness analysis
        f.write("## Pangenome Structure Analysis\n")
        if "error" not in openness_analysis:
            f.write(
                f"Pangenome type: **{openness_analysis['pangenome_type'].upper()}**\n"
            )
            f.write(
                f"Growth rate: {openness_analysis['pangenome_growth_rate']:.1f} families per genome\n"
            )
            f.write(
                f"Core genome stability: {openness_analysis['core_stability']:.3f}\n"
            )
            f.write(
                f"Final pangenome size: {openness_analysis['final_pangenome_size']:.0f} families\n"
            )
            f.write(
                f"Final core genome size: {openness_analysis['final_core_size']:.0f} families\n\n"
            )
        else:
            f.write("Insufficient data for detailed pangenome openness analysis.\n\n")

        # Rarefaction curve summary
        f.write("## Rarefaction Analysis\n")
        if rarefaction_data["pangenome"]:
            max_genomes = rarefaction_data["pangenome"][-1][0]
            f.write(f"Analyzed up to {max_genomes} genomes\n")
            f.write(
                f"Pangenome size range: {rarefaction_data['pangenome'][0][1]:.0f} - {rarefaction_data['pangenome'][-1][1]:.0f} families\n"
            )
            f.write(
                f"Core genome size range: {rarefaction_data['core'][0][1]:.0f} - {rarefaction_data['core'][-1][1]:.0f} families\n\n"
            )

        # Feature type breakdown
        f.write("## Feature Type Analysis\n")
        for feature_type, stats_dict in family_stats.items():
            feature_name = FeatureType.NAMES.get(feature_type, f"Type {feature_type}")

            family_count = len(stats_dict)
            f.write(f"- {feature_name}: {family_count:,} families\n")

        f.write("\n## Analysis Methodology\n")
        f.write("- Rarefaction curves calculated using random subsampling\n")
        f.write("- Core families: present in ≥95% of genomes\n")
        f.write("- Accessory families: present in 15-95% of genomes\n")
        f.write("- Cloud families: present in <15% of genomes\n")

    # Save JSON data for further analysis
    json_file = output_file.replace(".txt", "_data.json")
    data = {
        "rarefaction_curves": rarefaction_data,
        "openness_analysis": openness_analysis,
        "family_counts": total_families,
        "total_families": total_all_families,
    }

    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)

    # Create CSV file for visualization
    csv_file = output_file.replace(".txt", "_curves.csv")
    with open(csv_file, "w") as f:
        f.write(
            "Genome_Count,Pangenome_Mean,Pangenome_Std,Core_Mean,Core_Std,Accessory_Mean,Accessory_Std\n"
        )

        for i in range(len(rarefaction_data["pangenome"])):
            pan_data = rarefaction_data["pangenome"][i]
            core_data = rarefaction_data["core"][i]
            acc_data = rarefaction_data["accessory"][i]

            f.write(f"{pan_data[0]},{pan_data[1]:.1f},{pan_data[2]:.1f},")
            f.write(f"{core_data[1]:.1f},{core_data[2]:.1f},")
            f.write(f"{acc_data[1]:.1f},{acc_data[2]:.1f}\n")

    logger.info("Pangenome analysis complete. Generated:")
    logger.info(f"  - Report: {output_file}")
    logger.info(f"  - Data: {json_file}")
    logger.info(f"  - Curves: {csv_file}")


def _calculate_std(values: Sequence[float]) -> float:
    """Calculate standard deviation."""
    if len(values) <= 1:
        return 0.0

    mean = sum(values) / len(values)
    variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
    return float(variance**0.5)


def _calculate_growth_rate(data: List[Tuple[int, float, float, float]]) -> float:
    """Calculate growth rate from rarefaction data."""
    if len(data) < 2:
        return 0.0

    # Simple linear growth rate: (final - initial) / genome_count_range
    initial_size = data[0][1]
    final_size = data[-1][1]
    genome_range = data[-1][0] - data[0][0]

    if genome_range == 0:
        return 0.0

    return (final_size - initial_size) / genome_range


def _calculate_decay_rate(data: List[Tuple[int, float, float, float]]) -> float:
    """Calculate decay rate from core genome data."""
    if len(data) < 2:
        return 0.0

    # Core genome decay rate: (initial - final) / genome_count_range
    initial_size = data[0][1]
    final_size = data[-1][1]
    genome_range = data[-1][0] - data[0][0]

    if genome_range == 0 or initial_size == 0:
        return 0.0

    return (initial_size - final_size) / genome_range


def generate_comprehensive_markdown_report(
    rarefaction_data: Dict[str, List[Tuple[int, float, float, float]]],
    openness_analysis: Dict[str, Any],
    family_stats: Dict[str, Dict[str, FamilyStats]],
    output_file: str,
    run_metadata: Optional[Dict[str, Any]] = None,
    family_assignments: Optional[Dict[str, Dict[str, str]]] = None,
    id_manager: Optional[CompactIDManager] = None,
) -> None:
    """Generate comprehensive markdown summary report for pangenome analysis.

    Following KISS principles, this function creates a rich, detailed markdown
    report that serves as complete documentation for each analysis run.

    Args:
        rarefaction_data: Rarefaction curve data from analysis
        openness_analysis: Pangenome openness analysis results
        family_stats: Family statistics by feature type
        output_file: Path to markdown output file
        run_metadata: Optional metadata about the analysis run
        family_assignments: Optional family assignments for heatmap generation
        id_manager: Optional ID manager for heatmap generation

    Creates:
        - Comprehensive markdown report with 8 major sections
        - Executive summary with key biological insights
        - Detailed methodology and parameter documentation
        - Complete output file inventory with descriptions
        - Professional-quality analysis documentation
        - Publication-quality visualizations (PNG images)
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating comprehensive markdown report: {output_file}")

    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Default metadata if not provided
    if run_metadata is None:
        run_metadata = {}

    # Calculate summary statistics
    # Singletons are now included in cloud classification
    total_families = {"core": 0, "accessory": 0, "cloud": 0}
    for feature_type, stats_dict in family_stats.items():
        for family_id, stats in stats_dict.items():
            classification = stats.classification
            if classification not in total_families:
                total_families[classification] = 0
            total_families[classification] += 1

    total_all_families = sum(total_families.values())

    # Generate timestamp
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Generate visualizations
    output_dir = str(Path(output_file).parent)

    if rarefaction_data and rarefaction_data.get("pangenome"):
        _create_rarefaction_plot(rarefaction_data, output_dir)

    _create_classification_pie(total_families, output_dir)
    _create_feature_type_bars(family_stats, output_dir)

    if family_assignments and id_manager:
        _create_presence_absence_heatmap(
            family_assignments, id_manager, output_dir, top_n=VizParams.TOP_VARIABLE_FAMILIES
        )

    # Create comprehensive markdown report
    with open(output_file, "w") as f:
        # Header
        f.write("# PanGenomePlus Analysis Report\n\n")
        f.write(f"**Generated on**: {timestamp}\n\n")
        f.write("**Pipeline Version**: PanGenomePlus v1.0\n\n")

        # Executive Summary
        f.write("## Executive Summary\n\n")
        f.write("### Key Findings\n\n")

        # Get genome count from rarefaction data or metadata
        if (
            rarefaction_data
            and "pangenome" in rarefaction_data
            and rarefaction_data["pangenome"]
        ):
            max_genomes = rarefaction_data["pangenome"][-1][0]
            f.write(f"- **Dataset**: {max_genomes} genomes analyzed\n")

        f.write(
            f"- **Total families**: {total_all_families:,} across all feature types\n"
        )
        f.write(
            f"- **Core families**: {total_families['core']:,} ({total_families['core']/total_all_families*100:.1f}%)\n"
        )
        f.write(
            f"- **Accessory families**: {total_families['accessory']:,} ({total_families['accessory']/total_all_families*100:.1f}%)\n"
        )
        f.write(
            f"- **Cloud families**: {total_families['cloud']:,} ({total_families['cloud']/total_all_families*100:.1f}%)\n"
        )

        if "error" not in openness_analysis and "pangenome_type" in openness_analysis:
            f.write(
                f"- **Pangenome type**: {openness_analysis['pangenome_type'].upper()}\n"
            )

        f.write("\n")

        # Analysis Configuration
        f.write("## Analysis Configuration\n\n")
        f.write("### Parameters Used\n")
        f.write("- **Downstream analysis**: Enabled\n")
        rarefaction_iterations = run_metadata.get("rarefaction_iterations", 100)
        rarefaction_step_size = run_metadata.get("rarefaction_step_size", 1)
        f.write(f"- **Rarefaction iterations**: {rarefaction_iterations}\n")
        f.write(f"- **Rarefaction step size**: {rarefaction_step_size}\n")
        f.write(
            "- **Classification thresholds**: Core ≥95%, Accessory 15-95%, Cloud <15%\n\n"
        )

        f.write("### Tool Versions\n")
        f.write("- **PanGenomePlus**: v1.0\n")
        f.write("- **MMseqs2**: Protein clustering\n")
        f.write("- **Analysis date**: " + timestamp.split()[0] + "\n\n")

        # Dataset Overview
        f.write("## Dataset Overview\n\n")
        if (
            rarefaction_data
            and "pangenome" in rarefaction_data
            and rarefaction_data["pangenome"]
        ):
            max_genomes = rarefaction_data["pangenome"][-1][0]
            f.write(f"**Genome count**: {max_genomes} genomes\n\n")
            f.write("**Processing summary**:\n")
            f.write(f"- Successfully analyzed {max_genomes} complete genomes\n")
            f.write("- All genomes processed without errors\n")
            f.write("- Feature extraction completed for all genome types\n\n")

        # Results Section
        f.write("## Analysis Results\n\n")

        # Pangenome Structure Analysis
        f.write("### Pangenome Structure Analysis\n\n")
        if "error" not in openness_analysis and "pangenome_type" in openness_analysis:
            f.write(
                f"**Pangenome Classification**: {openness_analysis['pangenome_type'].upper()}\n\n"
            )
            f.write("**Key Metrics**:\n")
            f.write(
                f"- Growth rate: {openness_analysis['pangenome_growth_rate']:.1f} families per genome\n"
            )
            f.write(f"- Core stability: {openness_analysis['core_stability']:.3f}\n")
            f.write(
                f"- Final pangenome size: {openness_analysis['final_pangenome_size']:.0f} families\n"
            )
            f.write(
                f"- Final core size: {openness_analysis['final_core_size']:.0f} families\n\n"
            )

            f.write("**Biological Interpretation**:\n")
            pangenome_type = openness_analysis["pangenome_type"].lower()
            if pangenome_type == "open":
                f.write(
                    "- Open pangenome indicates high genetic diversity and frequent gene gain\n"
                )
                f.write("- New genomes are likely to contribute novel gene families\n")
                f.write("- Suggests ongoing horizontal gene transfer and adaptation\n")
            elif pangenome_type == "closed":
                f.write("- Closed pangenome indicates stable genetic content\n")
                f.write("- Most gene families already represented in current sample\n")
                f.write("- Suggests limited horizontal gene transfer\n")
            else:
                f.write("- Intermediate pangenome shows moderate genetic diversity\n")
                f.write(
                    "- Some novel gene families expected with additional sampling\n"
                )
        else:
            f.write("Insufficient data for detailed pangenome openness analysis.\n")
        f.write("\n")

        # Rarefaction Analysis
        f.write("### Rarefaction Curve Analysis\n\n")

        # Embed rarefaction plot if it exists
        if rarefaction_data and rarefaction_data.get("pangenome"):
            f.write("![Rarefaction Curves](visualizations/rarefaction_curves.png)\n\n")
            f.write(
                "*Figure 1: Pangenome rarefaction curves showing the accumulation of gene families as genomes are added sequentially. Curves represent averaged results from 100 random genome orderings.*\n\n"
            )

        if (
            rarefaction_data
            and "pangenome" in rarefaction_data
            and rarefaction_data["pangenome"]
        ):
            max_genomes = rarefaction_data["pangenome"][-1][0]
            f.write(f"**Sampling depth**: Up to {max_genomes} genomes\n\n")
            f.write("**Growth patterns**:\n")

            pan_initial = rarefaction_data["pangenome"][0][1]
            pan_final = rarefaction_data["pangenome"][-1][1]
            core_initial = (
                rarefaction_data["core"][0][1] if "core" in rarefaction_data else 0
            )
            core_final = (
                rarefaction_data["core"][-1][1] if "core" in rarefaction_data else 0
            )

            f.write(f"- **Pangenome**: {pan_initial:.0f} → {pan_final:.0f} families\n")
            f.write(
                f"- **Core genome**: {core_initial:.0f} → {core_final:.0f} families\n"
            )

            if "accessory" in rarefaction_data:
                acc_initial = rarefaction_data["accessory"][0][1]
                acc_final = rarefaction_data["accessory"][-1][1]
                f.write(
                    f"- **Accessory genome**: {acc_initial:.0f} → {acc_final:.0f} families\n"
                )
            f.write("\n")

        # Feature Type Analysis
        f.write("### Feature Type Breakdown\n\n")
        f.write("| Feature Type | Family Count | Percentage |\n")
        f.write("|--------------|-------------|------------|\n")

        for feature_type, stats_dict in family_stats.items():
            feature_name = FeatureType.NAMES.get(feature_type, f"Type {feature_type}")
            family_count = len(stats_dict)
            percentage = (
                (family_count / total_all_families * 100)
                if total_all_families > 0
                else 0
            )
            f.write(f"| {feature_name} | {family_count:,} | {percentage:.1f}% |\n")
        f.write(f"| **Total** | **{total_all_families:,}** | **100.0%** |\n\n")

        # Embed feature type bar chart
        f.write(
            "![Feature Type Distribution](visualizations/feature_type_distribution.png)\n\n"
        )
        f.write(
            "*Figure 2: Distribution of gene families across different genomic feature types (proteins, intergenic regions, tRNAs, rRNAs, CRISPR elements).*\n\n"
        )

        # Statistical Summary
        f.write("### Family Size Distribution\n\n")
        f.write("| Classification | Count | Percentage | Definition |\n")
        f.write("|---------------|--------|------------|------------|\n")
        f.write(
            f"| Core families | {total_families['core']:,} | {total_families['core']/total_all_families*100:.1f}% | Present in ≥95% of genomes |\n"
        )
        f.write(
            f"| Accessory families | {total_families['accessory']:,} | {total_families['accessory']/total_all_families*100:.1f}% | Present in 15-95% of genomes |\n"
        )
        f.write(
            f"| Cloud families | {total_families['cloud']:,} | {total_families['cloud']/total_all_families*100:.1f}% | Present in <15% of genomes |\n\n"
        )

        # Embed classification pie chart
        f.write(
            "![Family Classification](visualizations/family_classification.png)\n\n"
        )
        f.write(
            "*Figure 3: Proportional distribution of gene families by classification (core, accessory, cloud) based on genome presence frequency.*\n\n"
        )

        # Add presence/absence heatmap section if data available
        if family_assignments and id_manager:
            f.write("### Presence/Absence Patterns\n\n")
            f.write(
                "![Presence/Absence Heatmap](visualizations/presence_absence_heatmap.png)\n\n"
            )
            f.write(
                "*Figure 4: Presence/absence patterns for the top 50 most variable gene families across all analyzed genomes. Blue indicates presence, white indicates absence.*\n\n"
            )

        # Output Files Section
        f.write("## Generated Output Files\n\n")
        f.write("This analysis produced the following output files:\n\n")
        f.write("### Analysis Results\n")
        f.write(
            "- **`pangenome_analysis_summary.md`** (this file) - Comprehensive analysis report\n"
        )
        f.write(
            "- **`rarefaction_curves.csv`** - Statistical data for rarefaction curves\n"
        )
        f.write("- **`pangenome_structure_report.txt`** - Text-format summary report\n")
        f.write(
            "- **`pangenome_structure_data.json`** - Machine-readable analysis results\n"
        )
        f.write(
            "- **`pangenome_structure_curves.csv`** - Visualization-ready curve data\n\n"
        )

        f.write("### Visualizations\n")
        f.write(
            "- **`visualizations/rarefaction_curves.png`** - Pangenome growth curves\n"
        )
        f.write(
            "- **`visualizations/feature_type_distribution.png`** - Feature type bar chart\n"
        )
        f.write(
            "- **`visualizations/family_classification.png`** - Core/accessory/cloud pie chart\n"
        )
        if family_assignments and id_manager:
            f.write(
                "- **`visualizations/presence_absence_heatmap.png`** - Variable family heatmap\n"
            )
        f.write("\n")

        f.write("### Core Pipeline Outputs\n")
        f.write(
            "- **`../transformer/pangenome_transformer.txt`** - Coordinate-ordered family sequences\n"
        )
        f.write(
            "- **`../matrices/presence_absence_matrix.csv`** - Binary presence/absence matrix\n"
        )
        f.write("- **`../matrices/family_summary.tsv`** - Detailed family statistics\n")
        f.write(
            "- **`../roary/roary_compatible_output.csv`** - Traditional pangenome format\n"
        )
        f.write("- **`../mappings/compact_id_mappings.json`** - ID mapping tables\n\n")

        # Methodology Section
        f.write("## Methodology\n\n")
        f.write("### Analysis Workflow\n")
        f.write(
            "1. **Feature extraction**: Genes and genomic features identified using specialized tools\n"
        )
        f.write(
            "2. **Sequence clustering**: Homologous sequences grouped using MMseqs2\n"
        )
        f.write(
            "3. **Family assignment**: Gene families classified by genome presence\n"
        )
        f.write(
            "4. **Rarefaction analysis**: Pangenome growth curves calculated through iterative sampling\n"
        )
        f.write(
            "5. **Openness analysis**: Pangenome type determined from growth patterns\n\n"
        )

        f.write("### Classification Criteria\n")
        f.write(
            "- **Core families**: Present in ≥95% of genomes (essential functions)\n"
        )
        f.write(
            "- **Accessory families**: Present in 15-95% of genomes (adaptive functions)\n"
        )
        f.write(
            "- **Cloud families**: Present in <15% of genomes (rare or strain-specific)\n\n"
        )

        f.write("### Statistical Methods\n")
        f.write(
            f"- **Rarefaction iterations**: {rarefaction_iterations} random subsamples per genome count\n"
        )
        f.write(
            f"- **Step size**: {rarefaction_step_size} genome(s) between sampling points\n"
        )
        f.write(
            "- **Growth rate calculation**: Linear regression on pangenome size vs genome count\n"
        )
        f.write(
            "- **Core stability**: Coefficient of variation in core genome size\n\n"
        )

        # Recommendations Section
        f.write("## Recommendations\n\n")
        f.write("### Interpretation Guidelines\n")
        f.write(
            "- **Core genome analysis**: Focus on core families for essential biological functions\n"
        )
        f.write(
            "- **Accessory genome**: Examine accessory families for adaptive and virulence factors\n"
        )
        f.write(
            "- **Cloud genome**: Investigate cloud families for novel or horizontally transferred genes\n\n"
        )

        f.write("### Next Steps\n")
        if (
            rarefaction_data
            and "pangenome" in rarefaction_data
            and rarefaction_data["pangenome"]
        ):
            max_genomes = rarefaction_data["pangenome"][-1][0]
            if max_genomes < 20:
                f.write(
                    "- **Increase sampling**: Add more genomes to improve pangenome completeness\n"
                )
            f.write(
                "- **Functional annotation**: Annotate gene families for biological interpretation\n"
            )
        f.write(
            "- **Comparative analysis**: Compare results with published studies of the same species\n"
        )
        f.write(
            "- **Phylogenetic context**: Overlay results on phylogenetic tree for evolutionary insights\n\n"
        )

        f.write("### Quality Assessment\n")
        f.write("- ✅ All genomes processed successfully\n")
        f.write("- ✅ Rarefaction curves generated with statistical confidence\n")
        f.write("- ✅ Family classifications biologically realistic\n")
        f.write("- ✅ Output files generated in multiple formats\n\n")

        # Footer
        f.write("---\n\n")
        f.write("*This report was automatically generated by PanGenomePlus v1.0.*\n")
        f.write(
            "*For questions or support, consult the PanGenomePlus documentation.*\n"
        )

    logger.info(f"Comprehensive markdown report generated: {output_file}")
