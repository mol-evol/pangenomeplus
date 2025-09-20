"""Output generation module for different formats."""

import logging
from pathlib import Path
from typing import Dict, List, Any, Literal, Optional
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from pangenomeplus.core.types import GenomeFeature, FeatureType
import gzip


logger = logging.getLogger(__name__)


def create_efficient_lookups(
    features: List[GenomeFeature],
    gene_families: Dict[str, str]
) -> Dict[str, Any]:
    """
    Create efficient lookup tables for fast post-processing operations.

    This eliminates O(n²) and O(n³) operations by pre-computing mappings.

    Args:
        features: All genome features
        gene_families: Gene family assignments

    Returns:
        Dictionary containing efficient lookup structures
    """
    lookups = {
        "genome_to_features": {},
        "feature_to_genome": {},
        "genome_to_families": {},
        "family_to_genomes": {},
        "genome_stats": {}
    }

    # Build genome-to-features mapping (O(n) instead of repeated searches)
    for feature in features:
        genome_id = feature.genome_id

        # Genome to features
        if genome_id not in lookups["genome_to_features"]:
            lookups["genome_to_features"][genome_id] = []
        lookups["genome_to_features"][genome_id].append(feature)

        # Feature to genome (for fast reverse lookup)
        lookups["feature_to_genome"][feature.pangenomeplus_id] = genome_id

        # Family to genomes mapping
        family_id = gene_families.get(feature.pangenomeplus_id)
        if family_id:
            if family_id not in lookups["family_to_genomes"]:
                lookups["family_to_genomes"][family_id] = set()
            lookups["family_to_genomes"][family_id].add(genome_id)

            # Genome to families
            if genome_id not in lookups["genome_to_families"]:
                lookups["genome_to_families"][genome_id] = set()
            lookups["genome_to_families"][genome_id].add(family_id)

    # Pre-compute genome statistics (avoid repeated calculations)
    for genome_id, genome_features in lookups["genome_to_features"].items():
        # Sort features by position within each genome (needed for transformer format)
        genome_features.sort(key=lambda f: (f.contig_id, f.start))

        stats = {
            "total_features": len(genome_features),
            "protein_count": sum(
                1 for f in genome_features if f.feature_type.value == 'protein_coding_gene'
            ),
            "intergenic_count": sum(
                1 for f in genome_features if f.feature_type.value == 'intergenic_region'
            ),
            "trna_count": sum(1 for f in genome_features if f.feature_type.value == 'trna'),
            "rrna_count": sum(1 for f in genome_features if f.feature_type.value == 'rrna')
        }
        lookups["genome_stats"][genome_id] = stats

    logger.info(f"Created efficient lookups for {len(lookups['genome_to_features'])} genomes")
    return lookups


def _process_genome_transformer_format(
    genome_data: tuple
) -> tuple:
    """
    Process a single genome for transformer format generation.
    Designed for multiprocessing workers.

    Args:
        genome_data: Tuple of (genome_id, features, gene_families,
                     include_singletons, feature_type_prefix)

    Returns:
        Tuple of (genome_id, family_sequence)
    """
    genome_id, genome_features, gene_families, include_singletons, feature_type_prefix = genome_data

    # Helper function to add feature type prefix and optimize family ID
    def add_feature_prefix(family_id: str, feature: GenomeFeature) -> str:
        """Optimize family ID by removing redundant prefixes and adding feature type prefix."""
        # Optimize family ID by removing redundant prefixes
        optimized_id = family_id

        # Remove PGP_ prefix from singleton IDs
        if family_id.startswith("S_PGP_"):
            optimized_id = family_id[4:]  # Remove "PGP_" keeping "S_"
        elif family_id.startswith("SINGLETON_PGP_"):
            optimized_id = family_id.replace("SINGLETON_PGP_", "S_", 1)

        if not feature_type_prefix:
            return optimized_id

        if feature.feature_type == FeatureType.PROTEIN_CODING_GENE:
            return f"P_{optimized_id}"
        elif feature.feature_type == FeatureType.INTERGENIC_REGION:
            return f"I_{optimized_id}"
        elif feature.feature_type == FeatureType.TRNA:
            return f"T_{optimized_id}"
        elif feature.feature_type == FeatureType.RRNA:
            return f"R_{optimized_id}"
        elif feature.feature_type == FeatureType.CRISPR_ARRAY:
            return f"C_{optimized_id}"
        else:
            # For any other feature types, use generic prefix
            return f"G_{optimized_id}"

    family_sequence = []
    for feature in genome_features:  # Already sorted by position
        family_id = gene_families.get(feature.pangenomeplus_id)

        if family_id:
            if family_id.startswith("SINGLETON_"):
                if include_singletons:
                    prefixed_family = add_feature_prefix(family_id, feature)
                    family_sequence.append(prefixed_family)
            else:
                prefixed_family = add_feature_prefix(family_id, feature)
                family_sequence.append(prefixed_family)

    return genome_id, family_sequence


def _process_genome_statistics(
    genome_data: tuple
) -> tuple:
    """
    Process a single genome's statistics for report generation.
    Designed for multiprocessing workers.

    Args:
        genome_data: Tuple of (genome_id, genome_stats)

    Returns:
        Tuple of (genome_id, formatted_stats_dict)
    """
    genome_id, stats = genome_data

    return genome_id, {
        "genome_id": genome_id,
        "total_features": stats["total_features"],
        "protein_count": stats["protein_count"],
        "intergenic_count": stats["intergenic_count"],
        "trna_count": stats["trna_count"],
        "rrna_count": stats["rrna_count"]
    }


def generate_transformer_format_parallel(
    gene_families: Dict[str, str],
    efficient_lookups: Dict[str, Any],
    output_file: Path,
    include_singletons: bool = True,
    singleton_prefix: str = "SINGLETON",
    feature_type_prefix: bool = True,
    compress_output: bool = False
) -> Path:
    """
    Generate AI-ready transformer format output using parallel processing.

    Format: genome_name FAM_001 FAM_045 FAM_002 FAM_089...

    Args:
        gene_families: Gene family assignments
        efficient_lookups: Pre-computed lookup tables
        output_file: Output file path
        include_singletons: Whether to include singleton genes
        singleton_prefix: Prefix for singleton gene identifiers

    Returns:
        Path to generated file
    """
    genome_to_features = efficient_lookups["genome_to_features"]

    # Determine number of parallel processes
    max_workers = min(len(genome_to_features), multiprocessing.cpu_count(), 8)

    # For small datasets, don't use parallel processing overhead
    if len(genome_to_features) <= 2:
        logger.info("Small dataset - processing transformer format sequentially")
        return generate_transformer_format(
            gene_families, [], output_file, include_singletons,
            singleton_prefix, feature_type_prefix, compress_output
        )

    logger.info(
        f"Processing {len(genome_to_features)} genomes in parallel using {max_workers} workers"
    )

    # Prepare data for parallel processing
    genome_data_list = [
        (genome_id, genome_features, gene_families, include_singletons, feature_type_prefix)
        for genome_id, genome_features in genome_to_features.items()
    ]

    # Process genomes in parallel
    genome_results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all genome processing jobs
        future_to_genome = {
            executor.submit(_process_genome_transformer_format, genome_data): genome_data[0]
            for genome_data in genome_data_list
        }

        # Collect results as they complete
        for future in as_completed(future_to_genome):
            genome_id = future_to_genome[future]
            try:
                result_genome_id, family_sequence = future.result()
                genome_results[result_genome_id] = family_sequence
                logger.debug(f"Completed transformer format for {result_genome_id}")
            except Exception as e:
                logger.error(f"Failed to process transformer format for {genome_id}: {e}")
                genome_results[genome_id] = []  # Empty sequence on error

    # Write results to file (ensure consistent ordering)
    open_func = gzip.open if compress_output else open
    mode = 'wt' if compress_output else 'w'

    with open_func(output_file, mode) as f:
        for genome_id in sorted(genome_results.keys()):
            family_sequence = genome_results[genome_id]
            if family_sequence:
                f.write(f"{genome_id} {' '.join(family_sequence)}\n")

    compression_note = " (compressed)" if compress_output else ""
    logger.info(f"Generated parallel transformer format{compression_note}: {output_file}")
    return output_file


def generate_outputs(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    features: List[GenomeFeature],
    id_mappings: Dict[str, Any],
    output_dir: Path,
    formats: List[str],
    config: Optional[Dict[str, Any]] = None,
    pipeline_results: Optional[Dict[str, Any]] = None,
    annotation_results: Optional[Dict[str, Dict[str, Any]]] = None,
    cluster_results: Optional[Dict[Any, List[Any]]] = None
) -> Dict[str, Path]:
    """
    Generate all requested output formats.

    Args:
        gene_families: Gene family assignments
        family_members: Detailed family membership
        features: All genome features
        id_mappings: ID mapping dictionaries
        output_dir: Output directory
        formats: List of output formats to generate
        pipeline_results: Full pipeline results for comprehensive reporting
        annotation_results: Raw annotation results per genome
        cluster_results: Clustering results per feature type

    Returns:
        Dictionary mapping format names to output file paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files = {}

    # Create efficient lookup tables for fast processing (eliminates O(n²) and O(n³) operations)
    logger.info("Creating efficient lookup tables for fast post-processing...")
    efficient_lookups = create_efficient_lookups(features, gene_families)
    logger.info("Lookup tables created successfully")

    for format_name in formats:
        logger.info(f"Generating {format_name} output format")

        try:
            if format_name == "transformer":
                # Get transformer configuration
                transformer_config = config.get("output", {}).get("transformer", {}) if config else {}

                # Extract parameters with defaults
                include_singletons = transformer_config.get("include_singletons", True)
                singleton_prefix = transformer_config.get("singleton_prefix", "SINGLETON")
                feature_type_prefix = transformer_config.get("feature_type_prefix", True)
                compress_output = transformer_config.get("compress", False)

                output_file = generate_transformer_format_parallel(
                    gene_families=gene_families,
                    efficient_lookups=efficient_lookups,
                    output_file=output_dir / "transformer_format.txt",
                    include_singletons=include_singletons,
                    singleton_prefix=singleton_prefix,
                    feature_type_prefix=feature_type_prefix,
                    compress_output=compress_output
                )
                output_files["transformer"] = output_file

            elif format_name == "presence_absence":
                output_file = generate_presence_absence_matrix(
                    gene_families=gene_families,
                    features=features,
                    output_file=output_dir / "presence_absence_matrix.tsv"
                )
                output_files["presence_absence"] = output_file

            elif format_name == "roary":
                roary_files = generate_roary_compatible_output_optimized(
                    gene_families=gene_families,
                    family_members=family_members,
                    efficient_lookups=efficient_lookups,
                    output_dir=output_dir / "roary"
                )
                output_files.update(roary_files)

            elif format_name == "fasta":
                fasta_files = generate_fasta_outputs(
                    gene_families=gene_families,
                    features=features,
                    output_dir=output_dir / "sequences"
                )
                output_files.update(fasta_files)

            elif format_name == "markdown_report":
                report_file = generate_markdown_report_optimized(
                    gene_families=gene_families,
                    family_members=family_members,
                    efficient_lookups=efficient_lookups,
                    id_mappings=id_mappings,
                    output_dir=output_dir / "reports",
                    pipeline_results=pipeline_results,
                    annotation_results=annotation_results,
                    cluster_results=cluster_results
                )
                output_files["markdown_report"] = report_file

            else:
                logger.warning(f"Unknown output format: {format_name}")

        except Exception as e:
            logger.error(f"Failed to generate {format_name} output: {e}")

    logger.info(f"Generated {len(output_files)} output files")
    return output_files


def generate_transformer_format(
    gene_families: Dict[str, str],
    features: List[GenomeFeature],
    output_file: Path,
    include_singletons: bool = True,
    singleton_prefix: str = "SINGLETON",
    feature_type_prefix: bool = True,
    compress_output: bool = False
) -> Path:
    """
    Generate AI-ready transformer format output with feature type distinction.

    Format: genome_name [P_|I_|T_|R_|C_]FAM_001 [P_|I_|T_|R_|C_]FAM_045...

    Args:
        gene_families: Gene family assignments
        features: All genome features with positional information
        output_file: Output file path
        include_singletons: Whether to include singleton genes
        singleton_prefix: Prefix for singleton gene identifiers
        feature_type_prefix: Whether to add feature type prefixes (P_/I_/T_/R_/C_)
        compress_output: Whether to compress output with gzip

    Feature Type Prefixes:
        P_ = Protein-coding genes
        I_ = Intergenic regions
        T_ = tRNA genes
        R_ = rRNA genes
        C_ = CRISPR arrays

    Example Output (with feature_type_prefix=True):
        genome_A P_FAM_001 T_FAM_002 I_INT_001 R_FAM_003
        genome_B P_FAM_001 I_INT_002 C_FAM_004 P_FAM_002
    """
    # Group features by genome and sort by position
    genome_features = {}
    for feature in features:
        genome_id = feature.genome_id
        if genome_id not in genome_features:
            genome_features[genome_id] = []
        genome_features[genome_id].append(feature)

    # Sort features within each genome by contig and position
    for genome_id in genome_features:
        genome_features[genome_id].sort(key=lambda f: (f.contig_id, f.start))

    # Helper function to add feature type prefix and optimize family ID
    def add_feature_prefix(family_id: str, feature: GenomeFeature) -> str:
        """Optimize family ID by removing redundant prefixes and adding feature type prefix."""
        # Optimize family ID by removing redundant prefixes
        optimized_id = family_id

        # Remove PGP_ prefix from singleton IDs
        if family_id.startswith("S_PGP_"):
            optimized_id = family_id[4:]  # Remove "PGP_" keeping "S_"
        elif family_id.startswith("SINGLETON_PGP_"):
            optimized_id = family_id.replace("SINGLETON_PGP_", "S_", 1)

        if not feature_type_prefix:
            return optimized_id

        if feature.feature_type == FeatureType.PROTEIN_CODING_GENE:
            return f"P_{optimized_id}"
        elif feature.feature_type == FeatureType.INTERGENIC_REGION:
            return f"I_{optimized_id}"
        elif feature.feature_type == FeatureType.TRNA:
            return f"T_{optimized_id}"
        elif feature.feature_type == FeatureType.RRNA:
            return f"R_{optimized_id}"
        elif feature.feature_type == FeatureType.CRISPR_ARRAY:
            return f"C_{optimized_id}"
        else:
            # For any other feature types, use generic prefix
            return f"G_{optimized_id}"

    # Determine output mode (compressed or regular)
    open_func = gzip.open if compress_output else open
    mode = 'wt' if compress_output else 'w'

    with open_func(output_file, mode) as f:
        for genome_id, genome_feature_list in genome_features.items():
            family_sequence = []

            for feature in genome_feature_list:
                family_id = gene_families.get(feature.pangenomeplus_id)

                if family_id:
                    if family_id.startswith("SINGLETON_"):
                        if include_singletons:
                            prefixed_family = add_feature_prefix(family_id, feature)
                            family_sequence.append(prefixed_family)
                    else:
                        prefixed_family = add_feature_prefix(family_id, feature)
                        family_sequence.append(prefixed_family)

            if family_sequence:
                f.write(f"{genome_id} {' '.join(family_sequence)}\n")

    compression_note = " (compressed)" if compress_output else ""
    logger.info(f"Generated transformer format{compression_note}: {output_file}")
    return output_file


def generate_presence_absence_matrix(
    gene_families: Dict[str, str],
    features: List[GenomeFeature],
    output_file: Path,
    format_type: Literal["csv", "tsv", "phylip"] = "tsv"
) -> Path:
    """
    Generate presence/absence matrix compatible with Roary/Panaroo.

    Args:
        gene_families: Gene family assignments
        features: All genome features
        output_file: Output file path
        format_type: Output format

    Returns:
        Path to generated file
    """
    # Create genome-family matrix
    genomes = sorted(set(f.genome_id for f in features))
    families = sorted(set(gene_families.values()))

    # Initialize matrix
    matrix = pd.DataFrame(0, index=families, columns=genomes)

    # Fill matrix
    for feature in features:
        family_id = gene_families.get(feature.pangenomeplus_id)
        if family_id:
            matrix.loc[family_id, feature.genome_id] = 1

    # Save matrix
    if format_type == "csv":
        matrix.to_csv(output_file, sep=',')
    elif format_type == "tsv":
        matrix.to_csv(output_file, sep='\t')
    elif format_type == "phylip":
        # PHYLIP format is more complex, simplified version here
        with open(output_file, 'w') as f:
            f.write(f"{len(genomes)} {len(families)}\n")
            for genome in genomes:
                f.write(f"{genome[:10]:10s} {''.join(map(str, matrix[genome].values))}\n")

    logger.info(f"Generated presence/absence matrix: {output_file}")
    return output_file


def generate_roary_compatible_output(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    features: List[GenomeFeature],
    output_dir: Path
) -> Dict[str, Path]:
    """
    Generate Roary-compatible output files.

    Args:
        gene_families: Gene family assignments
        family_members: Family membership details
        features: All genome features
        output_dir: Output directory

    Returns:
        Dictionary of generated file paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {}

    # Main gene presence/absence table
    gene_pa_file = output_dir / "gene_presence_absence.csv"
    files["gene_presence_absence"] = create_roary_gene_table(
        gene_families, family_members, features, gene_pa_file
    )

    # Summary statistics
    summary_file = output_dir / "summary_statistics.txt"
    files["summary_statistics"] = create_summary_statistics(
        gene_families, family_members, features, summary_file
    )

    # Pangenome statistics for R
    r_files = create_r_compatible_tables(gene_families, features, output_dir)
    files.update(r_files)

    return files


def generate_roary_compatible_output_optimized(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    efficient_lookups: Dict[str, Any],
    output_dir: Path
) -> Dict[str, Path]:
    """
    Generate Roary-compatible output files using optimized processing.

    Args:
        gene_families: Gene family assignments
        family_members: Family membership details
        efficient_lookups: Pre-computed lookup tables
        output_dir: Output directory

    Returns:
        Dictionary of generated file paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {}

    # Main gene presence/absence table (optimized version)
    gene_pa_file = output_dir / "gene_presence_absence.csv"
    files["gene_presence_absence"] = create_roary_gene_table_optimized(
        gene_families, family_members, efficient_lookups, gene_pa_file
    )

    # Summary statistics (use efficient lookups)
    summary_file = output_dir / "summary_statistics.txt"
    files["summary_statistics"] = create_summary_statistics_optimized(
        gene_families, family_members, efficient_lookups, summary_file
    )

    # Pangenome statistics for R (optimized)
    r_files = create_r_compatible_tables_optimized(efficient_lookups, output_dir)
    files.update(r_files)

    return files


def create_summary_statistics_optimized(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    efficient_lookups: Dict[str, Any],
    output_file: Path,
    core_threshold: float = 0.95,
    shell_threshold: float = 0.15
) -> Path:
    """Create pangenome summary statistics using efficient lookups."""
    family_to_genomes = efficient_lookups["family_to_genomes"]
    genome_to_features = efficient_lookups["genome_to_features"]
    n_genomes = len(genome_to_features)

    # Count family types using efficient lookups
    core_families = 0
    shell_families = 0
    cloud_families = 0

    for family_id, members in family_members.items():
        if family_id.startswith("SINGLETON_"):
            cloud_families += 1
            continue

        # Use efficient lookup to get genomes for this family
        family_genomes = family_to_genomes.get(family_id, set())
        prevalence = len(family_genomes) / n_genomes

        if prevalence >= core_threshold:  # Core
            core_families += 1
        elif prevalence >= shell_threshold:  # Shell
            shell_families += 1
        else:  # Cloud
            cloud_families += 1

    with open(output_file, 'w') as f:
        f.write("Pangenome Analysis Summary\n")
        f.write("=" * 30 + "\n\n")
        f.write(f"Number of genomes: {n_genomes}\n")
        f.write(f"Total gene families: {len(family_members)}\n\n")
        f.write(f"Core families (≥95%): {core_families}\n")
        f.write(f"Shell families (15-95%): {shell_families}\n")
        f.write(f"Cloud families (<15%): {cloud_families}\n\n")
        f.write(f"Pangenome size: {len(family_members)}\n")
        f.write(f"Core genome size: {core_families}\n")

    return output_file


def create_r_compatible_tables_optimized(
    efficient_lookups: Dict[str, Any],
    output_dir: Path
) -> Dict[str, Path]:
    """Create R-compatible data tables using efficient lookups."""
    files = {}

    # Use pre-computed statistics from efficient lookups
    genome_to_features = efficient_lookups["genome_to_features"]

    # This would create the various .Rtab files that Roary generates
    # Simplified implementation for now
    conserved_file = output_dir / "number_of_conserved_genes.Rtab"
    with open(conserved_file, 'w') as f:
        f.write("Genes\tGenomes\n")
        # Would calculate actual conserved gene statistics using efficient lookups
        f.write(f"{len(genome_to_features) * 1000}\t1\n")  # Simplified placeholder

    files["conserved_genes"] = conserved_file

    return files


def create_roary_gene_table_optimized(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    efficient_lookups: Dict[str, Any],
    output_file: Path
) -> Path:
    """Create Roary-style gene presence/absence table using efficient lookups."""
    feature_to_genome = efficient_lookups["feature_to_genome"]
    genome_to_features = efficient_lookups["genome_to_features"]
    genomes = sorted(genome_to_features.keys())

    # Pre-compute ID mappings for faster lookup
    id_mappings = efficient_lookups.get("id_mappings", {})

    with open(output_file, 'w') as f:
        # Header
        header = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
                 "Avg sequences per isolate", "Genome fragment", "Order within fragment",
                 "Accessory fragment", "Accessory order with fragment", "QC", "Min group size nuc",
                 "Max group size nuc", "Avg group size nuc"] + genomes
        f.write(','.join(header) + '\n')

        # Family rows - process families in parallel for large datasets
        max_workers = min(len(family_members), multiprocessing.cpu_count(), 8)

        if len(family_members) > 100:  # Use parallel processing for large datasets
            logger.info(f"Processing {len(family_members)} gene families in parallel using {max_workers} workers")

            # Prepare data for parallel processing
            family_data_list = [
                (family_id, members, feature_to_genome, genomes)
                for family_id, members in family_members.items()
                if not family_id.startswith("SINGLETON_")
            ]

            # Process families in parallel
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all family processing jobs
                future_to_family = {
                    executor.submit(_process_family_for_roary, family_data): family_data[0]
                    for family_data in family_data_list
                }

                # Collect and write results as they complete
                for future in as_completed(future_to_family):
                    family_id = future_to_family[future]
                    try:
                        row_data = future.result()
                        f.write(','.join(row_data) + '\n')
                        logger.debug(f"Completed Roary row for family {family_id}")
                    except Exception as e:
                        logger.error(f"Failed to process family {family_id} for Roary table: {e}")

        else:
            # Sequential processing for small datasets
            for family_id, members in family_members.items():
                if family_id.startswith("SINGLETON_"):
                    continue

                # Use efficient lookup instead of O(n²) search
                genome_presence = {genome: [] for genome in genomes}
                for member in members:
                    # O(1) lookup instead of O(n) search through all features
                    genome_id = feature_to_genome.get(member["pangenomeplus_id"])
                    if genome_id:
                        genome_presence[genome_id].append(member["original_id"])

                # Basic statistics
                n_isolates = sum(1 for genes in genome_presence.values() if genes)
                n_sequences = len(members)
                avg_per_isolate = n_sequences / max(n_isolates, 1)

                # Representative annotation
                annotation = "hypothetical protein"  # Simplified
                for member in members:
                    if member["is_representative"]:
                        # Would extract actual annotation here
                        break

                # Build row
                row = [
                    family_id,  # Gene
                    family_id,  # Non-unique Gene name
                    annotation,  # Annotation
                    str(n_isolates),  # No. isolates
                    str(n_sequences),  # No. sequences
                    f"{avg_per_isolate:.2f}",  # Avg sequences per isolate
                    "",  # Genome fragment
                    "",  # Order within fragment
                    "",  # Accessory fragment
                    "",  # Accessory order with fragment
                    "",  # QC
                    "",  # Min group size nuc
                    "",  # Max group size nuc
                    ""   # Avg group size nuc
                ]

                # Add genome presence
                for genome in genomes:
                    genes = genome_presence[genome]
                    if genes:
                        row.append('\t'.join(genes))
                    else:
                        row.append("")

                f.write(','.join(row) + '\n')

    return output_file


def _process_family_for_roary(family_data: tuple) -> List[str]:
    """
    Process a single gene family for Roary table generation.
    Designed for multiprocessing workers.

    Args:
        family_data: Tuple of (family_id, members, feature_to_genome, genomes)

    Returns:
        List of strings representing the Roary table row
    """
    family_id, members, feature_to_genome, genomes = family_data

    # Use efficient lookup instead of O(n²) search
    genome_presence = {genome: [] for genome in genomes}
    for member in members:
        # O(1) lookup instead of O(n) search through all features
        genome_id = feature_to_genome.get(member["pangenomeplus_id"])
        if genome_id:
            genome_presence[genome_id].append(member["original_id"])

    # Basic statistics
    n_isolates = sum(1 for genes in genome_presence.values() if genes)
    n_sequences = len(members)
    avg_per_isolate = n_sequences / max(n_isolates, 1)

    # Representative annotation
    annotation = "hypothetical protein"  # Simplified
    for member in members:
        if member["is_representative"]:
            # Would extract actual annotation here
            break

    # Build row
    row = [
        family_id,  # Gene
        family_id,  # Non-unique Gene name
        annotation,  # Annotation
        str(n_isolates),  # No. isolates
        str(n_sequences),  # No. sequences
        f"{avg_per_isolate:.2f}",  # Avg sequences per isolate
        "",  # Genome fragment
        "",  # Order within fragment
        "",  # Accessory fragment
        "",  # Accessory order with fragment
        "",  # QC
        "",  # Min group size nuc
        "",  # Max group size nuc
        ""   # Avg group size nuc
    ]

    # Add genome presence
    for genome in genomes:
        genes = genome_presence[genome]
        if genes:
            row.append('\t'.join(genes))
        else:
            row.append("")

    return row


def create_roary_gene_table(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    features: List[GenomeFeature],
    output_file: Path
) -> Path:
    """Create Roary-style gene presence/absence table (legacy version)."""
    # This is the original slow implementation - kept for compatibility
    # In practice, the optimized version will be used
    genomes = sorted(set(f.genome_id for f in features))

    with open(output_file, 'w') as f:
        # Header
        header = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
                 "Avg sequences per isolate", "Genome fragment", "Order within fragment",
                 "Accessory fragment", "Accessory order with fragment", "QC", "Min group size nuc",
                 "Max group size nuc", "Avg group size nuc"] + genomes
        f.write(','.join(header) + '\n')

        # Family rows
        for family_id, members in family_members.items():
            if family_id.startswith("SINGLETON_"):
                continue

            # Count presence per genome
            genome_presence = {genome: [] for genome in genomes}
            for member in members:
                # Find genome for this member
                for feature in features:
                    if feature.pangenomeplus_id == member["pangenomeplus_id"]:
                        genome_presence[feature.genome_id].append(member["original_id"])
                        break

            # Basic statistics
            n_isolates = sum(1 for genes in genome_presence.values() if genes)
            n_sequences = len(members)
            avg_per_isolate = n_sequences / max(n_isolates, 1)

            # Representative annotation
            annotation = "hypothetical protein"  # Simplified
            for member in members:
                if member["is_representative"]:
                    # Would extract actual annotation here
                    break

            # Build row
            row = [
                family_id,  # Gene
                family_id,  # Non-unique Gene name
                annotation,  # Annotation
                str(n_isolates),  # No. isolates
                str(n_sequences),  # No. sequences
                f"{avg_per_isolate:.2f}",  # Avg sequences per isolate
                "",  # Genome fragment
                "",  # Order within fragment
                "",  # Accessory fragment
                "",  # Accessory order with fragment
                "",  # QC
                "",  # Min group size nuc
                "",  # Max group size nuc
                ""   # Avg group size nuc
            ]

            # Add genome presence
            for genome in genomes:
                genes = genome_presence[genome]
                if genes:
                    row.append('\t'.join(genes))
                else:
                    row.append("")

            f.write(','.join(row) + '\n')

    return output_file


def create_summary_statistics(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    features: List[GenomeFeature],
    output_file: Path,
    core_threshold: float = 0.95,
    shell_threshold: float = 0.15
) -> Path:
    """Create pangenome summary statistics."""
    genomes = set(f.genome_id for f in features)
    n_genomes = len(genomes)

    # Count family types
    core_families = 0
    shell_families = 0
    cloud_families = 0

    for family_id, members in family_members.items():
        if family_id.startswith("SINGLETON_"):
            cloud_families += 1
            continue

        # Count genomes for this family
        family_genomes = set()
        for member in members:
            for feature in features:
                if feature.pangenomeplus_id == member["pangenomeplus_id"]:
                    family_genomes.add(feature.genome_id)
                    break

        prevalence = len(family_genomes) / n_genomes

        if prevalence >= core_threshold:  # Core
            core_families += 1
        elif prevalence >= shell_threshold:  # Shell
            shell_families += 1
        else:  # Cloud
            cloud_families += 1

    with open(output_file, 'w') as f:
        f.write("Pangenome Analysis Summary\n")
        f.write("=" * 30 + "\n\n")
        f.write(f"Number of genomes: {n_genomes}\n")
        f.write(f"Total gene families: {len(family_members)}\n\n")
        f.write(f"Core families (≥95%): {core_families}\n")
        f.write(f"Shell families (15-95%): {shell_families}\n")
        f.write(f"Cloud families (<15%): {cloud_families}\n\n")
        f.write(f"Pangenome size: {len(family_members)}\n")
        f.write(f"Core genome size: {core_families}\n")

    return output_file


def create_r_compatible_tables(
    gene_families: Dict[str, str],
    features: List[GenomeFeature],
    output_dir: Path
) -> Dict[str, Path]:
    """Create R-compatible data tables."""
    files = {}

    # This would create the various .Rtab files that Roary generates
    # Simplified implementation for now

    conserved_file = output_dir / "number_of_conserved_genes.Rtab"
    with open(conserved_file, 'w') as f:
        f.write("Genes\tGenomes\n")
        # Would calculate actual conserved gene statistics
        f.write("1000\t1\n")

    files["conserved_genes"] = conserved_file

    return files


def generate_fasta_outputs(
    gene_families: Dict[str, str],
    features: List[GenomeFeature],
    output_dir: Path
) -> Dict[str, Path]:
    """
    Generate FASTA format outputs.

    Args:
        gene_families: Gene family assignments
        features: All genome features
        output_dir: Output directory

    Returns:
        Dictionary of generated FASTA files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {}

    # Representative sequences (one per family)
    repr_file = output_dir / "representative_sequences.fasta"
    with open(repr_file, 'w') as f:
        # This would extract actual sequences
        # For now, create placeholder
        f.write(">placeholder\nMSTARTSEQUENCE\n")

    files["representative_sequences"] = repr_file

    logger.info(f"Generated FASTA outputs in {output_dir}")
    return files


def calculate_pangenome_statistics(
    presence_absence_matrix: pd.DataFrame,
    core_threshold: float = 0.95,
    shell_threshold: float = 0.15
) -> Dict[str, Any]:
    """
    Calculate comprehensive pangenome statistics.

    Args:
        presence_absence_matrix: Binary presence/absence matrix
        core_threshold: Threshold for core gene classification
        shell_threshold: Threshold for shell gene classification

    Returns:
        Dictionary with pangenome statistics
    """
    n_genomes = presence_absence_matrix.shape[1]
    n_families = presence_absence_matrix.shape[0]

    # Calculate prevalence for each family
    prevalence = presence_absence_matrix.sum(axis=1) / n_genomes

    # Classify families
    core_families = (prevalence >= core_threshold).sum()
    shell_families = ((prevalence >= shell_threshold) & (prevalence < core_threshold)).sum()
    cloud_families = (prevalence < shell_threshold).sum()

    statistics = {
        "n_genomes": n_genomes,
        "n_families": n_families,
        "core_families": core_families,
        "shell_families": shell_families,
        "cloud_families": cloud_families,
        "pangenome_size": n_families,
        "core_genome_size": core_families,
        "accessory_genome_size": shell_families + cloud_families
    }

    return statistics


def generate_markdown_report(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    features: List[GenomeFeature],
    id_mappings: Dict[str, Any],
    output_dir: Path,
    pipeline_results: Optional[Dict[str, Any]] = None,
    annotation_results: Optional[Dict[str, Dict[str, Any]]] = None,
    cluster_results: Optional[Dict[Any, List[Any]]] = None
) -> Path:
    """
    Generate comprehensive markdown report for pangenome analysis.

    Args:
        gene_families: Gene family assignments
        family_members: Detailed family membership
        features: All genome features
        id_mappings: ID mapping dictionaries
        output_dir: Output directory for report
        pipeline_results: Full pipeline results (runtime, config, etc.)
        annotation_results: Raw annotation results per genome
        cluster_results: Clustering results per feature type

    Returns:
        Path to generated markdown report
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    report_file = output_dir / "pangenome_analysis_report.md"

    # Collect basic statistics
    genomes = sorted(set(f.genome_id for f in features))
    n_genomes = len(genomes)
    n_features = len(features)
    n_families = len(set(gene_families.values()))

    # Group features by type and genome
    features_by_type = {}
    features_by_genome = {}

    for feature in features:
        # By type
        if feature.feature_type not in features_by_type:
            features_by_type[feature.feature_type] = []
        features_by_type[feature.feature_type].append(feature)

        # By genome
        if feature.genome_id not in features_by_genome:
            features_by_genome[feature.genome_id] = []
        features_by_genome[feature.genome_id].append(feature)

    # Calculate pangenome statistics
    presence_matrix = _create_presence_matrix(gene_families, features)
    pangenome_stats = calculate_pangenome_statistics(presence_matrix)

    with open(report_file, 'w') as f:
        # Title and metadata
        f.write("# PanGenome PLUS Analysis Report\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        if pipeline_results:
            f.write(f"**Pipeline Version:** {pipeline_results.get('pipeline_version', 'Unknown')}\n")
            f.write(f"**Analysis Strategy:** {pipeline_results.get('strategy', 'Unknown')}\n")
            f.write(f"**Runtime:** {pipeline_results.get('runtime_formatted', 'Unknown')}\n\n")

        f.write("---\n\n")

        # 1. Executive Summary
        f.write("## 📊 Executive Summary\n\n")
        f.write(f"- **Genomes Analyzed:** {n_genomes:,}\n")
        f.write(f"- **Total Features:** {n_features:,}\n")
        f.write(f"- **Gene Families:** {n_families:,}\n")
        f.write(f"- **Core Families:** {pangenome_stats['core_families']:,} ({pangenome_stats['core_families']/n_families*100:.1f}%)\n")
        f.write(f"- **Accessory Families:** {pangenome_stats['accessory_genome_size']:,} ({pangenome_stats['accessory_genome_size']/n_families*100:.1f}%)\n\n")

        # 2. Input Data Analysis
        f.write("## 🧬 Input Data Analysis\n\n")
        f.write("### Genome Information\n\n")
        f.write("| Genome ID | Features | Protein Genes | Intergenic Regions |\n")
        f.write("|-----------|----------|---------------|--------------------|\n")

        for genome_id in sorted(genomes):
            genome_features = features_by_genome.get(genome_id, [])
            protein_count = sum(1 for f in genome_features if f.feature_type.value == 'protein_coding_gene')
            intergenic_count = sum(1 for f in genome_features if f.feature_type.value == 'intergenic_region')
            f.write(f"| {genome_id} | {len(genome_features):,} | {protein_count:,} | {intergenic_count:,} |\n")

        f.write("\n")

        # Feature type distribution
        f.write("### Feature Type Distribution\n\n")
        for feature_type, type_features in features_by_type.items():
            percentage = len(type_features) / n_features * 100
            f.write(f"- **{feature_type.value.replace('_', ' ').title()}:** {len(type_features):,} features ({percentage:.1f}%)\n")
        f.write("\n")

        # 3. Annotation Results
        if annotation_results:
            f.write("## 🔧 Annotation Results\n\n")
            f.write("### Tools Used\n\n")

            # Extract tool information from annotation results
            tools_used = set()
            tool_success = {}

            for genome_id, genome_data in annotation_results.items():
                for tool, tool_data in genome_data.get('annotation_files', {}).items():
                    tools_used.add(tool)
                    if tool not in tool_success:
                        tool_success[tool] = {'success': 0, 'total': 0}
                    tool_success[tool]['total'] += 1
                    if 'error' not in tool_data:
                        tool_success[tool]['success'] += 1

            for tool in sorted(tools_used):
                success_rate = tool_success[tool]['success'] / tool_success[tool]['total'] * 100
                f.write(f"- **{tool.upper()}:** {tool_success[tool]['success']}/{tool_success[tool]['total']} genomes ({success_rate:.1f}% success)\n")
            f.write("\n")

        # 4. Clustering Analysis
        if cluster_results:
            f.write("## 🎯 Clustering Analysis\n\n")

            # Add clustering parameters if available from pipeline config
            if pipeline_results and 'config' in pipeline_results:
                config = pipeline_results['config']
                clustering_config = config.get('clustering', {})

                f.write("### Clustering Parameters\n\n")

                # Protein clustering parameters
                protein_config = clustering_config.get('protein', {})
                if protein_config:
                    f.write("**Protein Clustering:**\n")
                    f.write(f"- Identity threshold: {protein_config.get('identity', 0.8)*100:.0f}%\n")
                    f.write(f"- Coverage threshold: {protein_config.get('coverage', 0.8)*100:.0f}%\n")
                    f.write(f"- Sensitivity: {protein_config.get('sensitivity', 7.5)}\n")
                    f.write(f"- GPU acceleration: {'Yes' if protein_config.get('use_gpu', False) else 'No'}\n\n")

                # Intergenic clustering parameters
                intergenic_config = clustering_config.get('intergenic', {})
                if intergenic_config:
                    f.write("**Intergenic Clustering:**\n")
                    f.write(f"- Identity threshold: {intergenic_config.get('identity', 0.7)*100:.0f}%\n")
                    f.write(f"- Coverage threshold: {intergenic_config.get('coverage', 0.7)*100:.0f}%\n")
                    f.write(f"- Sensitivity: {intergenic_config.get('sensitivity', 6.0)}\n")
                    f.write(f"- GPU acceleration: {'Yes' if intergenic_config.get('use_gpu', False) else 'No'}\n\n")

            f.write("### Clustering Results by Feature Type\n\n")
            f.write("| Feature Type | Features | Clusters | Avg Cluster Size | Largest Cluster |\n")
            f.write("|--------------|----------|----------|------------------|------------------|\n")

            for feature_type, clusters in cluster_results.items():
                type_features = features_by_type.get(feature_type, [])
                n_clusters = len(clusters)
                avg_size = len(type_features) / n_clusters if n_clusters > 0 else 0

                # Find largest cluster size
                largest_cluster = max(len(members) for members in clusters.values()) if clusters else 0

                f.write(f"| {feature_type.value.replace('_', ' ').title()} | {len(type_features):,} | {n_clusters:,} | {avg_size:.1f} | {largest_cluster} |\n")

                # Add validation warning for unrealistic protein clustering
                if (feature_type.value == 'protein_coding_gene' and avg_size > 10):
                    f.write(f"\n⚠️ **Warning**: Average protein cluster size ({avg_size:.1f}) is unusually large. Expected range: 2-6 proteins per family for bacterial pangenomes.\n")

            f.write("\n")

        # 5. Gene Family Analysis
        f.write("## 👥 Gene Family Analysis\n\n")
        f.write("### Pangenome Classification\n\n")
        f.write(f"- **Core Genome (≥95% genomes):** {pangenome_stats['core_families']:,} families\n")
        f.write(f"- **Shell Genome (15-95% genomes):** {pangenome_stats['shell_families']:,} families\n")
        f.write(f"- **Cloud Genome (<15% genomes):** {pangenome_stats['cloud_families']:,} families\n\n")

        f.write(f"**Pangenome Size:** {pangenome_stats['pangenome_size']:,} families\n")
        f.write(f"**Core Genome Size:** {pangenome_stats['core_genome_size']:,} families\n\n")

        # Family size distribution
        family_sizes = []
        for family_id, members in family_members.items():
            if not family_id.startswith("SINGLETON_"):
                family_sizes.append(len(members))

        if family_sizes:
            f.write("### Family Size Statistics\n\n")
            f.write(f"- **Total Families:** {len(family_sizes):,}\n")
            f.write(f"- **Mean Family Size:** {np.mean(family_sizes):.1f}\n")
            f.write(f"- **Median Family Size:** {np.median(family_sizes):.1f}\n")
            f.write(f"- **Largest Family:** {max(family_sizes):,} members\n")
            f.write(f"- **Singleton Families:** {sum(1 for fid in family_members.keys() if fid.startswith('SINGLETON_')):,}\n\n")

        # 6. Performance Metrics
        if pipeline_results:
            f.write("## ⚡ Performance Metrics\n\n")
            runtime_seconds = pipeline_results.get('runtime_seconds', 0)
            f.write(f"- **Total Runtime:** {pipeline_results.get('runtime_formatted', 'Unknown')}\n")

            # Avoid division by zero errors
            if runtime_seconds > 0:
                runtime_minutes = runtime_seconds / 60
                if runtime_minutes > 0:
                    f.write(f"- **Processing Rate:** {n_genomes / runtime_minutes:.1f} genomes/minute\n")
                f.write(f"- **Feature Processing Rate:** {n_features / runtime_seconds:.0f} features/second\n")
            else:
                f.write("- **Processing Rate:** Unable to calculate (runtime too short)\n")

            if 'strategy' in pipeline_results:
                f.write(f"- **Strategy Used:** {pipeline_results['strategy']}\n")
            f.write("\n")

        # 7. Output Files Summary
        f.write("## 📁 Output Files\n\n")
        f.write("This analysis generated the following output files:\n\n")
        f.write("- **Transformer Format:** AI-ready sequence format for machine learning\n")
        f.write("- **Presence/Absence Matrix:** Binary matrix for phylogenetic analysis\n")
        f.write("- **Gene Family Assignments:** Detailed family membership information\n")
        f.write("- **Feature Summary:** Complete catalog of genomic features\n")
        f.write("- **This Report:** Comprehensive analysis summary\n\n")

        # 8. Methods Summary
        f.write("## 🔬 Methods Summary\n\n")
        f.write("### Annotation Pipeline\n")
        f.write("- **Gene Prediction:** Prodigal v2.6.3 with bacterial genetic code\n")
        f.write("- **tRNA Detection:** tRNAscan-SE with bacterial search mode\n")
        f.write("- **rRNA Detection:** Barrnap with bacterial HMM models\n")
        f.write("- **Intergenic Extraction:** Regions ≥50bp between annotated features\n\n")

        f.write("### Clustering Analysis\n")
        f.write("- **Algorithm:** MMseqs2 with GPU acceleration\n")
        f.write("- **Protein Identity:** 80% with 80% coverage\n")
        f.write("- **Intergenic Identity:** 70% with 70% coverage\n")
        f.write("- **Sensitivity:** High sensitivity mode (7.5 for proteins, 6.0 for intergenic)\n\n")

        f.write("### Gene Family Assignment\n")
        f.write("- **Core Genes:** Present in ≥95% of genomes\n")
        f.write("- **Shell Genes:** Present in 15-95% of genomes\n")
        f.write("- **Cloud Genes:** Present in <15% of genomes\n")
        f.write("- **Singletons:** Unique to individual genomes\n\n")

        # Footer
        f.write("---\n\n")
        f.write("*Report generated by PanGenome PLUS v1.0.0*\n")
        f.write("*For questions or support, please refer to the documentation*\n")

    logger.info(f"Generated comprehensive markdown report: {report_file}")
    return report_file


def _create_presence_matrix(gene_families: Dict[str, str], features: List[GenomeFeature]) -> pd.DataFrame:
    """Create presence/absence matrix for pangenome statistics."""
    genomes = sorted(set(f.genome_id for f in features))
    families = sorted(set(gene_families.values()))

    # Initialize matrix
    matrix = pd.DataFrame(0, index=families, columns=genomes)

    # Fill matrix
    for feature in features:
        family_id = gene_families.get(feature.pangenomeplus_id)
        if family_id:
            matrix.loc[family_id, feature.genome_id] = 1

    return matrix


def generate_markdown_report_optimized(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[Dict[str, Any]]],
    efficient_lookups: Dict[str, Any],
    id_mappings: Dict[str, Any],
    output_dir: Path,
    pipeline_results: Optional[Dict[str, Any]] = None,
    annotation_results: Optional[Dict[str, Dict[str, Any]]] = None,
    cluster_results: Optional[Dict[Any, List[Any]]] = None
) -> Path:
    """
    Generate comprehensive markdown report using parallel processing and efficient lookups.

    Args:
        gene_families: Gene family assignments
        family_members: Detailed family membership
        efficient_lookups: Pre-computed lookup tables
        id_mappings: ID mapping dictionaries
        output_dir: Output directory for report
        pipeline_results: Full pipeline results (runtime, config, etc.)
        annotation_results: Raw annotation results per genome
        cluster_results: Clustering results per feature type

    Returns:
        Path to generated markdown report
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    report_file = output_dir / "pangenome_analysis_report.md"

    # Get efficient lookups data
    genome_to_features = efficient_lookups["genome_to_features"]
    genome_stats = efficient_lookups["genome_stats"]
    family_to_genomes = efficient_lookups["family_to_genomes"]

    # Collect basic statistics using efficient lookups
    genomes = sorted(genome_to_features.keys())
    n_genomes = len(genomes)
    n_features = sum(len(features) for features in genome_to_features.values())
    n_families = len(set(gene_families.values()))

    # Group features by type using pre-computed stats
    features_by_type = {}
    for genome_id, stats in genome_stats.items():
        for feature_type in ['protein_coding_gene', 'intergenic_region', 'trna', 'rrna']:
            if feature_type not in features_by_type:
                features_by_type[feature_type] = 0
            count_key = feature_type.replace('_coding_gene', '').replace('_region', '') + '_count'
            features_by_type[feature_type] += stats.get(count_key, 0)

    # Calculate pangenome statistics using efficient lookups
    pangenome_stats = calculate_pangenome_statistics_optimized(family_to_genomes, family_members, n_genomes)

    # Generate genome statistics in parallel
    max_workers = min(len(genomes), multiprocessing.cpu_count(), 8)
    genome_stat_results = {}

    if len(genomes) > 4:  # Use parallel processing for larger datasets
        logger.info(f"Generating genome statistics in parallel using {max_workers} workers")

        genome_data_list = [
            (genome_id, genome_stats[genome_id])
            for genome_id in genomes
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all genome statistics jobs
            future_to_genome = {
                executor.submit(_process_genome_statistics, genome_data): genome_data[0]
                for genome_data in genome_data_list
            }

            # Collect results as they complete
            for future in as_completed(future_to_genome):
                genome_id = future_to_genome[future]
                try:
                    result_genome_id, stats = future.result()
                    genome_stat_results[result_genome_id] = stats
                except Exception as e:
                    logger.error(f"Failed to process statistics for {genome_id}: {e}")
                    # Use original stats as fallback
                    genome_stat_results[genome_id] = {
                        "genome_id": genome_id,
                        "total_features": genome_stats[genome_id]["total_features"],
                        "protein_count": genome_stats[genome_id]["protein_count"],
                        "intergenic_count": genome_stats[genome_id]["intergenic_count"],
                        "trna_count": genome_stats[genome_id].get("trna_count", 0),
                        "rrna_count": genome_stats[genome_id].get("rrna_count", 0)
                    }
    else:
        # Sequential processing for small datasets
        for genome_id in genomes:
            genome_stat_results[genome_id] = {
                "genome_id": genome_id,
                "total_features": genome_stats[genome_id]["total_features"],
                "protein_count": genome_stats[genome_id]["protein_count"],
                "intergenic_count": genome_stats[genome_id]["intergenic_count"],
                "trna_count": genome_stats[genome_id].get("trna_count", 0),
                "rrna_count": genome_stats[genome_id].get("rrna_count", 0)
            }

    with open(report_file, 'w') as f:
        # Title and metadata
        f.write("# PanGenome PLUS Analysis Report\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        if pipeline_results:
            f.write(f"**Pipeline Version:** {pipeline_results.get('pipeline_version', 'Unknown')}\n")
            f.write(f"**Analysis Strategy:** {pipeline_results.get('strategy', 'Unknown')}\n")
            f.write(f"**Runtime:** {pipeline_results.get('runtime_formatted', 'Unknown')}\n\n")

        f.write("---\n\n")

        # 1. Executive Summary
        f.write("## 📊 Executive Summary\n\n")
        f.write(f"- **Genomes Analyzed:** {n_genomes:,}\n")
        f.write(f"- **Total Features:** {n_features:,}\n")
        f.write(f"- **Gene Families:** {n_families:,}\n")
        f.write(f"- **Core Families:** {pangenome_stats['core_families']:,} ({pangenome_stats['core_families']/n_families*100:.1f}%)\n")
        f.write(f"- **Accessory Families:** {pangenome_stats['accessory_genome_size']:,} ({pangenome_stats['accessory_genome_size']/n_families*100:.1f}%)\n\n")

        # 2. Input Data Analysis (using parallel-generated statistics)
        f.write("## 🧬 Input Data Analysis\n\n")
        f.write("### Genome Information\n\n")
        f.write("| Genome ID | Features | Protein Genes | Intergenic Regions |\n")
        f.write("|-----------|----------|---------------|--------------------|\n")

        for genome_id in sorted(genomes):
            stats = genome_stat_results[genome_id]
            f.write(f"| {genome_id} | {stats['total_features']:,} | {stats['protein_count']:,} | {stats['intergenic_count']:,} |\n")

        f.write("\n")

        # Feature type distribution (using efficient calculations)
        f.write("### Feature Type Distribution\n\n")
        for feature_type, count in features_by_type.items():
            percentage = count / n_features * 100 if n_features > 0 else 0
            f.write(f"- **{feature_type.replace('_', ' ').title()}:** {count:,} features ({percentage:.1f}%)\n")
        f.write("\n")

        # 3. Annotation Results
        if annotation_results:
            f.write("## 🔧 Annotation Results\n\n")
            f.write("### Tools Used\n\n")

            # Extract tool information from annotation results
            tools_used = set()
            tool_success = {}

            for genome_id, genome_data in annotation_results.items():
                for tool, tool_data in genome_data.get('annotation_files', {}).items():
                    tools_used.add(tool)
                    if tool not in tool_success:
                        tool_success[tool] = {'success': 0, 'total': 0}
                    tool_success[tool]['total'] += 1
                    if 'error' not in tool_data:
                        tool_success[tool]['success'] += 1

            for tool in sorted(tools_used):
                success_rate = tool_success[tool]['success'] / tool_success[tool]['total'] * 100
                f.write(f"- **{tool.upper()}:** {tool_success[tool]['success']}/{tool_success[tool]['total']} genomes ({success_rate:.1f}% success)\n")
            f.write("\n")

        # 4. Gene Family Analysis
        f.write("## 👥 Gene Family Analysis\n\n")
        f.write("### Pangenome Classification\n\n")
        f.write(f"- **Core Genome (≥95% genomes):** {pangenome_stats['core_families']:,} families\n")
        f.write(f"- **Shell Genome (15-95% genomes):** {pangenome_stats['shell_families']:,} families\n")
        f.write(f"- **Cloud Genome (<15% genomes):** {pangenome_stats['cloud_families']:,} families\n\n")

        f.write(f"**Pangenome Size:** {pangenome_stats['pangenome_size']:,} families\n")
        f.write(f"**Core Genome Size:** {pangenome_stats['core_genome_size']:,} families\n\n")

        # Family size distribution
        family_sizes = []
        for family_id, members in family_members.items():
            if not family_id.startswith("SINGLETON_"):
                family_sizes.append(len(members))

        if family_sizes:
            f.write("### Family Size Statistics\n\n")
            f.write(f"- **Total Families:** {len(family_sizes):,}\n")
            f.write(f"- **Mean Family Size:** {np.mean(family_sizes):.1f}\n")
            f.write(f"- **Median Family Size:** {np.median(family_sizes):.1f}\n")
            f.write(f"- **Largest Family:** {max(family_sizes):,} members\n")
            f.write(f"- **Singleton Families:** {sum(1 for fid in family_members.keys() if fid.startswith('SINGLETON_')):,}\n\n")

        # 5. Performance Metrics
        if pipeline_results:
            f.write("## ⚡ Performance Metrics\n\n")
            runtime_seconds = pipeline_results.get('runtime_seconds', 0)
            f.write(f"- **Total Runtime:** {pipeline_results.get('runtime_formatted', 'Unknown')}\n")

            # Avoid division by zero errors
            if runtime_seconds > 0:
                runtime_minutes = runtime_seconds / 60
                if runtime_minutes > 0:
                    f.write(f"- **Processing Rate:** {n_genomes / runtime_minutes:.1f} genomes/minute\n")
                f.write(f"- **Feature Processing Rate:** {n_features / runtime_seconds:.0f} features/second\n")
            else:
                f.write("- **Processing Rate:** Unable to calculate (runtime too short)\n")

            if 'strategy' in pipeline_results:
                f.write(f"- **Strategy Used:** {pipeline_results['strategy']}\n")
            f.write("\n")

        # Footer
        f.write("---\n\n")
        f.write("*Report generated by PanGenome PLUS v1.0.0 with parallel processing optimizations*\n")
        f.write("*For questions or support, please refer to the documentation*\n")

    logger.info(f"Generated optimized markdown report: {report_file}")
    return report_file


def calculate_pangenome_statistics_optimized(
    family_to_genomes: Dict[str, set],
    family_members: Dict[str, List[Dict[str, Any]]],
    n_genomes: int,
    core_threshold: float = 0.95,
    shell_threshold: float = 0.15
) -> Dict[str, Any]:
    """
    Calculate comprehensive pangenome statistics using efficient lookups.

    Args:
        family_to_genomes: Pre-computed mapping of families to genomes
        family_members: Family member details
        n_genomes: Total number of genomes
        core_threshold: Threshold for core gene classification
        shell_threshold: Threshold for shell gene classification

    Returns:
        Dictionary with pangenome statistics
    """
    # Count family types using efficient lookups
    core_families = 0
    shell_families = 0
    cloud_families = 0

    for family_id, members in family_members.items():
        if family_id.startswith("SINGLETON_"):
            cloud_families += 1
            continue

        # Use efficient lookup to get genomes for this family
        family_genomes = family_to_genomes.get(family_id, set())
        prevalence = len(family_genomes) / n_genomes

        if prevalence >= core_threshold:  # Core
            core_families += 1
        elif prevalence >= shell_threshold:  # Shell
            shell_families += 1
        else:  # Cloud
            cloud_families += 1

    statistics = {
        "n_genomes": n_genomes,
        "n_families": len(family_members),
        "core_families": core_families,
        "shell_families": shell_families,
        "cloud_families": cloud_families,
        "pangenome_size": len(family_members),
        "core_genome_size": core_families,
        "accessory_genome_size": shell_families + cloud_families
    }

    return statistics