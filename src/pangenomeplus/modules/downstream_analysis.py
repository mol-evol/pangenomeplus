#!/usr/bin/env python3
"""
Downstream analysis module for pangenome data.

This module provides functionality for:
- Efficient binary matrix creation from gene family data
- Feature-type separation (protein-coding vs intergenic)
- Core/shell/cloud classification
- Statistical analysis and visualization
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import json
from typing import Dict, Tuple, List, Optional
from scipy import stats
from scipy.spatial.distance import pdist, squareform

logger = logging.getLogger(__name__)


def extract_genome_ids(gene_ids: pd.Series, pattern: str = r'(E_coli_\d+)') -> pd.Series:
    """
    Extract genome IDs from gene IDs using regex pattern.

    Args:
        gene_ids: Series of gene identifiers
        pattern: Regex pattern to extract genome ID

    Returns:
        Series of extracted genome IDs
    """
    return gene_ids.str.extract(pattern)[0]


def create_efficient_binary_matrix(
    gene_to_family: pd.DataFrame,
    family_summary: pd.DataFrame,
    genome_id_pattern: str = r'(E_coli_\d+)'
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """
    Create binary presence/absence matrix efficiently using pivot_table.

    Args:
        gene_to_family: DataFrame with gene_id and family_id columns
        family_summary: DataFrame with family information
        genome_id_pattern: Regex pattern to extract genome IDs

    Returns:
        Tuple of (binary_matrix, genomes_list, families_list)
    """
    logger.info("Creating binary presence/absence matrices efficiently...")

    # Extract genome IDs from gene IDs
    gene_to_family = gene_to_family.copy()
    gene_to_family['genome_id'] = extract_genome_ids(gene_to_family['gene_id'], genome_id_pattern)

    # Get unique genomes and families
    genomes = sorted(gene_to_family['genome_id'].unique())
    families = family_summary['family_id'].tolist()

    logger.info(f"Found {len(genomes)} genomes and {len(families)} families")

    # Add presence indicator
    gene_to_family['presence'] = 1

    # Create binary matrix efficiently using pivot_table
    binary_matrix = gene_to_family.pivot_table(
        index='family_id',
        columns='genome_id',
        values='presence',
        fill_value=0,
        aggfunc='max'  # In case of duplicates, take max (1)
    )

    # Reindex to ensure all families are included
    binary_matrix = binary_matrix.reindex(families, fill_value=0)
    binary_matrix = binary_matrix.reindex(columns=genomes, fill_value=0)

    logger.info(f"Binary matrix shape: {binary_matrix.shape}")

    return binary_matrix, genomes, families


def separate_by_feature_type(
    binary_matrix: pd.DataFrame,
    family_summary: pd.DataFrame,
    feature_type_column: str = 'feature_types'
) -> Dict[str, pd.DataFrame]:
    """
    Separate binary matrix by feature types.

    Args:
        binary_matrix: Binary presence/absence matrix
        family_summary: DataFrame with family information
        feature_type_column: Column name containing feature types

    Returns:
        Dictionary with feature types as keys and matrices as values
    """
    logger.info("Separating families by feature type...")

    # Get feature type information from family summary
    protein_families = family_summary[
        family_summary[feature_type_column] == 'protein_coding_gene'
    ]['family_id'].tolist()

    intergenic_families = family_summary[
        family_summary[feature_type_column] == 'intergenic_region'
    ]['family_id'].tolist()

    # Create separated matrices
    protein_matrix = binary_matrix.loc[protein_families]
    intergenic_matrix = binary_matrix.loc[intergenic_families]

    logger.info(f"Protein-coding families: {len(protein_families)}")
    logger.info(f"Intergenic families: {len(intergenic_families)}")

    return {
        'overall': binary_matrix,
        'protein': protein_matrix,
        'intergenic': intergenic_matrix
    }


def calculate_pangenome_statistics(
    matrix: pd.DataFrame,
    core_threshold: float = 0.95,
    cloud_threshold: float = 0.15
) -> Dict[str, float]:
    """
    Calculate core/shell/cloud statistics for a binary matrix.

    Args:
        matrix: Binary presence/absence matrix
        core_threshold: Threshold for core families (fraction of genomes)
        cloud_threshold: Threshold for cloud families (fraction of genomes)

    Returns:
        Dictionary with statistics
    """
    family_counts = matrix.sum(axis=1)  # Count genomes per family
    total_genomes = matrix.shape[1]

    core_families = (family_counts >= core_threshold * total_genomes).sum()
    cloud_families = (family_counts <= cloud_threshold * total_genomes).sum()
    shell_families = len(family_counts) - core_families - cloud_families

    stats = {
        'total_families': int(len(family_counts)),
        'core_families': int(core_families),
        'shell_families': int(shell_families),
        'cloud_families': int(cloud_families),
        'core_percentage': float((core_families / len(family_counts)) * 100),
        'shell_percentage': float((shell_families / len(family_counts)) * 100),
        'cloud_percentage': float((cloud_families / len(family_counts)) * 100),
        'total_genomes': int(total_genomes)
    }

    return stats


def analyze_feature_types(
    matrices: Dict[str, pd.DataFrame],
    names: Optional[Dict[str, str]] = None
) -> Dict[str, Dict[str, float]]:
    """
    Calculate statistics for each feature type.

    Args:
        matrices: Dictionary of matrices by feature type
        names: Optional mapping of keys to display names

    Returns:
        Dictionary of statistics by feature type
    """
    if names is None:
        names = {
            'overall': 'Overall',
            'protein': 'Protein-coding',
            'intergenic': 'Intergenic'
        }

    results = {}

    for key, matrix in matrices.items():
        if key in names:
            stats = calculate_pangenome_statistics(matrix)
            results[key] = stats

            name = names[key]
            logger.info(f"{name} statistics:")
            logger.info(f"  Total families: {stats['total_families']}")
            logger.info(
                f"  Core families: {stats['core_families']} ({stats['core_percentage']:.1f}%)"
            )
            logger.info(
                f"  Shell families: {stats['shell_families']} ({stats['shell_percentage']:.1f}%)"
            )
            logger.info(
                f"  Cloud families: {stats['cloud_families']} ({stats['cloud_percentage']:.1f}%)"
            )

    return results


def create_feature_type_comparison_plot(
    stats: Dict[str, Dict[str, float]],
    output_path: Path,
    title: str = "Feature Type Comparison"
) -> None:
    """
    Create a comparison plot for protein-coding vs intergenic features.

    Args:
        stats: Statistics dictionary from analyze_feature_types
        output_path: Path to save the plot
        title: Plot title
    """
    logger.info("Creating comparison plot...")

    categories = ['Core', 'Shell', 'Cloud']
    protein_data = [
        stats['protein']['core_percentage'],
        stats['protein']['shell_percentage'],
        stats['protein']['cloud_percentage']
    ]
    intergenic_data = [
        stats['intergenic']['core_percentage'],
        stats['intergenic']['shell_percentage'],
        stats['intergenic']['cloud_percentage']
    ]

    x = np.arange(len(categories))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(
        x - width/2, protein_data, width, label='Protein-coding', color='#FF6B6B', alpha=0.8
    )
    bars2 = ax.bar(x + width/2, intergenic_data, width, label='Intergenic', color='#4ECDC4', alpha=0.8)

    ax.set_xlabel('Genome Category')
    ax.set_ylabel('Percentage of Families')
    ax.set_title(title, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def run_feature_type_analysis(
    family_summary_file: Path,
    gene_to_family_file: Path,
    output_dir: Path,
    genome_id_pattern: str = r'(E_coli_\d+)'
) -> Dict[str, Dict[str, float]]:
    """
    Run complete feature-type separated analysis.

    Args:
        family_summary_file: Path to family summary TSV
        gene_to_family_file: Path to gene-to-family mapping TSV
        output_dir: Directory to save results
        genome_id_pattern: Regex pattern for genome ID extraction

    Returns:
        Dictionary of statistics by feature type
    """
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading gene family data...")
    family_summary = pd.read_csv(family_summary_file, sep='\t')
    gene_to_family = pd.read_csv(gene_to_family_file, sep='\t')

    logger.info(f"Loaded {len(family_summary)} gene families")
    logger.info(f"Loaded {len(gene_to_family)} gene-to-family mappings")

    # Create binary matrix
    binary_matrix, genomes, families = create_efficient_binary_matrix(
        gene_to_family, family_summary, genome_id_pattern
    )

    # Separate by feature type
    matrices = separate_by_feature_type(binary_matrix, family_summary)

    # Calculate statistics
    stats = analyze_feature_types(matrices)

    # Save results
    with open(output_dir / "pangenome_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)

    # Create plot
    create_feature_type_comparison_plot(
        stats,
        output_dir / "feature_type_comparison.png",
        f"Feature Type Comparison ({len(genomes)} Genomes)"
    )

    logger.info(f"Analysis completed! Results saved to {output_dir}")

    return stats


def calculate_pairwise_distances(matrix: pd.DataFrame, feature_type: str) -> Dict:
    """
    Calculate pairwise Jaccard distances for genomes.

    Args:
        matrix: Binary presence/absence matrix (families x genomes)
        feature_type: Type of features being analyzed

    Returns:
        Dictionary containing distance statistics
    """
    logger.info(f"Calculating {feature_type} pairwise distances...")

    # Transpose matrix for genome-by-gene format
    genome_matrix = matrix.T

    # Calculate Jaccard distances
    distances = pdist(genome_matrix, metric='jaccard')
    distance_matrix = squareform(distances)

    logger.info(f"Distance matrix: {distance_matrix.shape}")
    logger.info(f"Mean distance: {np.mean(distances):.3f}")
    logger.info(f"Distance range: {np.min(distances):.3f} - {np.max(distances):.3f}")

    return {
        'distance_vector': distances,
        'distance_matrix': distance_matrix,
        'mean_distance': np.mean(distances),
        'std_distance': np.std(distances),
        'min_distance': np.min(distances),
        'max_distance': np.max(distances)
    }


def analyze_distance_correlation(protein_distances: Dict, intergenic_distances: Dict) -> Dict:
    """
    Analyze correlation between protein and intergenic distances.

    Args:
        protein_distances: Distance results for protein-coding genes
        intergenic_distances: Distance results for intergenic regions

    Returns:
        Dictionary containing correlation statistics
    """
    logger.info("Analyzing distance correlation between protein-coding and intergenic regions...")

    # Calculate Pearson correlation
    correlation, p_value = stats.pearsonr(
        protein_distances['distance_vector'],
        intergenic_distances['distance_vector']
    )

    # Calculate Spearman correlation (rank-based)
    spearman_corr, spearman_p = stats.spearmanr(
        protein_distances['distance_vector'],
        intergenic_distances['distance_vector']
    )

    logger.info(f"Pearson correlation: r = {correlation:.3f} (p = {p_value:.2e})")
    logger.info(f"Spearman correlation: ρ = {spearman_corr:.3f} (p = {spearman_p:.2e})")

    return {
        'pearson_correlation': correlation,
        'pearson_p_value': p_value,
        'spearman_correlation': spearman_corr,
        'spearman_p_value': spearman_p,
        'is_significant': p_value < 0.05
    }


def run_pairwise_distance_analysis(
    family_summary_file: Path,
    gene_to_family_file: Path,
    output_dir: Path,
    genome_id_pattern: str = r'(E_coli_\d+)'
) -> Dict:
    """
    Run comprehensive pairwise distance analysis.

    Args:
        family_summary_file: Path to family summary TSV
        gene_to_family_file: Path to gene-to-family mapping TSV
        output_dir: Output directory for results
        genome_id_pattern: Regex pattern to extract genome IDs

    Returns:
        Dictionary containing all analysis results
    """
    logger.info("Starting pairwise distance analysis...")

    # Load data
    family_summary = pd.read_csv(family_summary_file, sep='\t')
    gene_to_family = pd.read_csv(gene_to_family_file, sep='\t')

    # Create binary matrix
    binary_matrix, genomes, families = create_efficient_binary_matrix(
        gene_to_family, family_summary, genome_id_pattern
    )

    # Separate by feature type
    matrices = separate_by_feature_type(binary_matrix, family_summary)

    # Calculate distances for each feature type
    results = {}
    for feature_type, matrix in matrices.items():
        if len(matrix) > 0:  # Only if matrix has data
            results[feature_type] = calculate_pairwise_distances(matrix, feature_type)

    # Analyze correlation if both feature types are present
    if 'protein_coding_gene' in results and 'intergenic_region' in results:
        correlation_results = analyze_distance_correlation(
            results['protein_coding_gene'],
            results['intergenic_region']
        )
        results['correlation'] = correlation_results

    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "pairwise_distance_analysis.json", 'w') as f:
        # Convert numpy arrays to lists for JSON serialization
        serializable_results = {}
        for key, value in results.items():
            if isinstance(value, dict):
                serializable_results[key] = {}
                for sub_key, sub_value in value.items():
                    if isinstance(sub_value, np.ndarray):
                        serializable_results[key][sub_key] = sub_value.tolist()
                    else:
                        serializable_results[key][sub_key] = sub_value
            else:
                serializable_results[key] = value

        json.dump(serializable_results, f, indent=2)

    logger.info(f"Pairwise distance analysis completed! Results saved to {output_dir}")

    return results