"""Output format generation for PanGenomePlus.

This module generates multiple output formats from clustering results:
- Transformer format: coordinate-ordered family sequences per genome
- Presence/absence matrices: binary CSV/TSV matrices
- Family summaries: comprehensive statistics and classifications
- Roary-compatible format: traditional pangenome analysis format

All outputs maintain genomic coordinate ordering and biological accuracy.
"""

import csv
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from .clustering import FamilyStats
from .compact_ids import CompactIDManager
from .constants import CLOUD_GENOME_THRESHOLD, CORE_GENOME_THRESHOLD, FeatureType


def reconstruct_genomic_coordinates(
    genome_id: str,
    family_ids: List[str],
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
) -> List[Tuple[str, int]]:
    """Reconstruct genomic coordinates for families in a genome.

    Args:
        genome_id: Target genome identifier
        family_ids: List of family IDs present in genome
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for coordinate lookups

    Returns:
        List of (family_id, start_coordinate) tuples sorted by coordinate

    Raises:
        ValueError: If family mapping cannot be reconstructed
    """
    # Create reverse lookup: family_id -> list of compact_ids
    family_to_compact_ids: Dict[str, List[str]] = {}
    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            if family_id not in family_to_compact_ids:
                family_to_compact_ids[family_id] = []
            family_to_compact_ids[family_id].append(compact_id)

    family_coordinates = []

    for family_id in family_ids:
        # Get compact IDs for this family
        compact_ids = family_to_compact_ids.get(family_id, [])

        # Find the first occurrence in the target genome
        for compact_id in compact_ids:
            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature and feature.genome_id == genome_id:
                start_coord = feature.start
                family_coordinates.append((family_id, start_coord))
                break  # Use first occurrence for coordinate ordering

    # Sort by genomic coordinate (NOT alphabetical)
    family_coordinates.sort(key=lambda x: x[1])

    return family_coordinates


def generate_transformer_format(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    output_dir: str,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Generate coordinate-ordered family sequences per genome.

    Creates transformer format: 'genome_name FAM_P1 FAM_I5 FAM_P2 FAM_T42'
    Families are ordered by genomic coordinate, NOT alphabetical order.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for coordinate reconstruction
        output_dir: Output directory path
        logger: Optional logger instance

    Returns:
        Path to generated transformer format file

    Raises:
        IOError: If output file cannot be written
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating transformer format with coordinate ordering")

    # Create output directory
    transformer_dir = Path(output_dir) / "transformer"
    transformer_dir.mkdir(parents=True, exist_ok=True)

    # Collect all families per genome
    genome_families: Dict[str, Set[str]] = {}  # genome_id -> set of family_ids

    # Process all feature types
    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature:
                genome_id = feature.genome_id
                if genome_id not in genome_families:
                    genome_families[genome_id] = set()
                genome_families[genome_id].add(family_id)

    # Generate coordinate-ordered output
    output_file = transformer_dir / "pangenome_transformer.txt"

    with open(output_file, "w") as f:
        for genome_id in sorted(genome_families.keys()):
            family_list = list(genome_families[genome_id])

            # Reconstruct coordinates and sort by position
            family_coords = reconstruct_genomic_coordinates(
                genome_id, family_list, family_assignments, id_manager
            )

            # Create coordinate-ordered token sequence
            ordered_families = [family_id for family_id, _ in family_coords]

            # Write transformer format line
            line = f"{genome_id} " + " ".join(ordered_families) + "\n"
            f.write(line)

    logger.info(f"Transformer format written to {output_file}")
    logger.info(f"Generated {len(genome_families)} genome sequences")

    return str(output_file)


def generate_presence_absence_matrix(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    output_dir: str,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Generate presence/absence matrix in CSV format.

    Creates binary matrix with genomes as rows and families as columns.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for genome identification
        output_dir: Output directory path
        logger: Optional logger instance

    Returns:
        Path to generated presence/absence matrix file
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating presence/absence matrix")

    # Create output directory
    matrices_dir = Path(output_dir) / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    # Collect all families and genomes
    all_families: Set[str] = set()
    genome_families: Dict[str, Set[str]] = {}  # genome_id -> set of family_ids

    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            all_families.add(family_id)

            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature:
                genome_id = feature.genome_id
                if genome_id not in genome_families:
                    genome_families[genome_id] = set()
                genome_families[genome_id].add(family_id)

    # Sort families and genomes for consistent output
    sorted_families = sorted(all_families)
    sorted_genomes = sorted(genome_families.keys())

    # Generate matrix
    output_file = matrices_dir / "presence_absence_matrix.csv"

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)

        # Write header
        header = ["Genome"] + sorted_families
        writer.writerow(header)

        # Write matrix rows
        for genome_id in sorted_genomes:
            genome_set = genome_families.get(genome_id, set())
            row: List[str] = [genome_id]

            for family_id in sorted_families:
                presence = "1" if family_id in genome_set else "0"
                row.append(presence)

            writer.writerow(row)

    logger.info(f"Presence/absence matrix written to {output_file}")
    logger.info(
        f"Matrix dimensions: {len(sorted_genomes)} genomes × {len(sorted_families)} families"
    )

    return str(output_file)


def generate_family_summary(
    family_stats: Dict[str, Dict[str, FamilyStats]],
    output_dir: str,
    total_genomes: int,
    id_manager: Optional[CompactIDManager] = None,
    family_assignments: Optional[Dict[str, Dict[str, str]]] = None,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Generate comprehensive family summary statistics with sequence information.

    Creates TSV file with family classifications, statistics, and sequence metrics.

    Args:
        family_stats: Dict mapping feature types to family_id -> FamilyStats
        output_dir: Output directory path
        total_genomes: Total number of genomes for percentage calculations
        id_manager: Optional CompactIDManager for sequence statistics
        family_assignments: Optional Dict mapping feature types to compact_id -> family_id
        logger: Optional logger instance

    Returns:
        Path to generated family summary file
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating family summary statistics")

    # Create output directory
    matrices_dir = Path(output_dir) / "matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    # Generate summary
    output_file = matrices_dir / "family_summary.tsv"

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Enhanced header with sequence statistics
        if id_manager is not None:
            header = [
                "Family_ID",
                "Feature_Type",
                "Size",
                "Genome_Count",
                "Genome_Percentage",
                "Classification",
                "Representative",
                "Rep_Sequence_Length",
                "Min_Length",
                "Max_Length",
                "Mean_Length",
                "Sequence_Diversity",
                "Length_Variance",
            ]
        else:
            header = [
                "Family_ID",
                "Feature_Type",
                "Size",
                "Genome_Count",
                "Genome_Percentage",
                "Classification",
                "Representative",
            ]
        writer.writerow(header)

        # Process all feature types
        total_families = 0
        classification_counts = {"core": 0, "accessory": 0, "cloud": 0}
        classification_by_type: Dict[str, Dict[str, int]] = {}

        for feature_type, stats_dict in family_stats.items():
            # Initialize per-type classification counts
            if feature_type not in classification_by_type:
                classification_by_type[feature_type] = {
                    "core": 0,
                    "accessory": 0,
                    "cloud": 0,
                }

            for family_id, stats in stats_dict.items():
                genome_percentage = (stats.genome_count / total_genomes) * 100

                # Classification logic using constants
                if genome_percentage >= (CORE_GENOME_THRESHOLD * 100):
                    classification = "core"
                elif genome_percentage >= (CLOUD_GENOME_THRESHOLD * 100):
                    classification = "accessory"
                else:
                    classification = "cloud"

                classification_counts[classification] += 1
                classification_by_type[feature_type][classification] += 1
                total_families += 1

                # Basic family row
                row = [
                    family_id,
                    feature_type,
                    stats.member_count,
                    stats.genome_count,
                    f"{genome_percentage:.1f}%",
                    classification,
                    stats.representative_id,
                ]

                # Add sequence statistics if id_manager and family_assignments are available
                if id_manager is not None and family_assignments is not None:
                    # Find members of this family
                    family_members = []
                    for compact_id, assigned_family_id in family_assignments.get(
                        feature_type, {}
                    ).items():
                        if assigned_family_id == family_id:
                            family_members.append(compact_id)

                    seq_stats = _calculate_family_sequence_stats(
                        family_id, family_members, id_manager, logger
                    )
                    row.extend(
                        [
                            seq_stats.get("rep_length", "N/A"),
                            seq_stats.get("min_length", "N/A"),
                            seq_stats.get("max_length", "N/A"),
                            seq_stats.get("mean_length", "N/A"),
                            seq_stats.get("diversity", "N/A"),
                            seq_stats.get("variance", "N/A"),
                        ]
                    )

                writer.writerow(row)

    logger.info(f"Family summary written to {output_file}")
    logger.info(f"Total families: {total_families}")
    logger.info(
        f"Total - Core: {classification_counts['core']}, "
        f"Accessory: {classification_counts['accessory']}, "
        f"Cloud: {classification_counts['cloud']}"
    )

    # Log per-feature-type breakdown
    logger.info("Classification by feature type:")
    for feature_type, counts in sorted(classification_by_type.items()):
        feature_name = FeatureType.NAMES.get(feature_type, feature_type)
        logger.info(
            f"  {feature_name}: "
            f"Core={counts['core']}, "
            f"Accessory={counts['accessory']}, "
            f"Cloud={counts['cloud']}"
        )

    return str(output_file)


def generate_roary_compatible_output(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    output_dir: str,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Generate Roary-compatible output format.

    Creates CSV format compatible with existing Roary analysis tools.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for original gene name lookup
        output_dir: Output directory path
        logger: Optional logger instance

    Returns:
        Path to generated Roary-compatible file
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating Roary-compatible output")

    # Create output directory
    roary_dir = Path(output_dir) / "roary"
    roary_dir.mkdir(parents=True, exist_ok=True)

    # This is a simplified Roary format - full implementation would require
    # more detailed gene information and clustering statistics
    output_file = roary_dir / "roary_compatible_output.csv"

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)

        # Write simplified header
        header = [
            "Gene",
            "Non-unique Gene name",
            "Annotation",
            "No. isolates",
            "No. sequences",
        ]
        writer.writerow(header)

        # Process families and write simplified entries
        family_count = 0
        for feature_type, assignments in family_assignments.items():
            family_genes: Dict[str, List[str]] = {}  # family_id -> list of compact_ids

            for compact_id, family_id in assignments.items():
                if family_id not in family_genes:
                    family_genes[family_id] = []
                family_genes[family_id].append(compact_id)

            for family_id, compact_ids in family_genes.items():
                # Get representative gene info
                rep_id = compact_ids[0]
                feature = id_manager.get_feature_by_compact_id(rep_id)

                if feature:
                    original_id = feature.original_id
                    product = feature.metadata.get("product", "hypothetical protein")

                    # Count unique genomes for this family
                    unique_genomes = set()
                    for cid in compact_ids:
                        feat = id_manager.get_feature_by_compact_id(cid)
                        if feat:
                            unique_genomes.add(feat.genome_id)

                    row = [
                        family_id,
                        original_id,
                        product,
                        len(unique_genomes),
                        len(compact_ids),
                    ]
                    writer.writerow(row)
                    family_count += 1

    logger.info(f"Roary-compatible output written to {output_file}")
    logger.info(f"Generated {family_count} family entries")

    return str(output_file)


def generate_representative_sequences_fasta(
    family_assignments: Dict[str, Dict[str, str]],
    family_stats: Dict[str, Dict[str, FamilyStats]],
    id_manager: CompactIDManager,
    output_dir: str,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, str]:
    """Generate FASTA files with representative sequences for each family.

    Since sequences are no longer stored in JSON mappings for space optimization,
    this function generates FASTA files containing representative sequences
    for each gene family, organized by feature type.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        family_stats: Dict mapping feature types to family_id -> FamilyStats
        id_manager: CompactIDManager for sequence access
        output_dir: Output directory for FASTA files
        logger: Optional logger instance

    Returns:
        Dict mapping feature types to output FASTA file paths

    Raises:
        FileNotFoundError: If sequences cannot be accessed
        IOError: If file writing fails
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating representative sequence FASTA files")

    # Create sequences output directory
    sequences_dir = os.path.join(output_dir, "representative_sequences")
    Path(sequences_dir).mkdir(parents=True, exist_ok=True)

    output_files = {}

    # Process each feature type separately
    for feature_type, assignments in family_assignments.items():
        if not assignments:
            logger.debug(f"No families found for feature type {feature_type}")
            continue

        # Get feature type name for filename
        feature_name = FeatureType.NAMES.get(
            feature_type, f"type_{feature_type}"
        ).lower()

        output_file = os.path.join(
            sequences_dir, f"representatives_{feature_name}.fasta"
        )

        # Collect family representatives
        family_representatives = {}
        for compact_id, family_id in assignments.items():
            if family_id not in family_representatives:
                # Use first member as representative (could be improved to use actual cluster representative)
                family_representatives[family_id] = compact_id

        # Write FASTA file
        sequences_written = 0
        with open(output_file, "w") as f:
            for family_id, representative_compact_id in sorted(
                family_representatives.items()
            ):
                # Get feature with sequence
                feature = id_manager.get_feature_by_compact_id(
                    representative_compact_id
                )
                if feature and feature.sequence:
                    f.write(f">{family_id}\n")
                    f.write(f"{feature.sequence}\n")
                    sequences_written += 1
                else:
                    logger.warning(
                        f"No sequence found for representative {representative_compact_id} in family {family_id}"
                    )

        output_files[feature_type] = output_file
        logger.info(
            f"Generated {sequences_written} representative sequences for {feature_name}: {output_file}"
        )

    logger.info(
        f"Representative sequence FASTA generation complete: {len(output_files)} files created"
    )
    return output_files


def generate_all_outputs(
    family_assignments: Dict[str, Dict[str, str]],
    family_stats: Dict[str, Dict[str, FamilyStats]],
    id_manager: CompactIDManager,
    output_dir: str,
    total_genomes: int,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """Generate all output formats from clustering results.

    Orchestrates generation of transformer, matrix, summary, and Roary formats.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        family_stats: Dict mapping feature types to family_id -> FamilyStats
        id_manager: CompactIDManager for coordinate and metadata lookup
        output_dir: Output directory path
        total_genomes: Total number of genomes processed
        logger: Optional logger instance

    Returns:
        Dict mapping format names to output file paths
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Generating all output formats")

    # Create main output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    output_files: Dict[str, Any] = {}

    try:
        # Generate transformer format
        output_files["transformer"] = generate_transformer_format(
            family_assignments, id_manager, output_dir, logger
        )

        # Generate presence/absence matrix
        output_files["matrix"] = generate_presence_absence_matrix(
            family_assignments, id_manager, output_dir, logger
        )

        # Generate family summary
        output_files["summary"] = generate_family_summary(
            family_stats,
            output_dir,
            total_genomes,
            id_manager,
            family_assignments,
            logger,
        )

        # Generate Roary-compatible output
        output_files["roary"] = generate_roary_compatible_output(
            family_assignments, id_manager, output_dir, logger
        )

        # Generate representative sequence FASTA files
        rep_seq_files = generate_representative_sequences_fasta(
            family_assignments, family_stats, id_manager, output_dir, logger
        )
        output_files["representative_sequences"] = rep_seq_files

        logger.info("All output formats generated successfully")

    except Exception as e:
        logger.error(f"Output generation failed: {e}")
        raise

    return output_files


def _calculate_family_sequence_stats(
    family_id: str,
    family_members: List[str],
    id_manager: CompactIDManager,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """Calculate sequence statistics for a gene family.

    Args:
        family_id: Family identifier
        family_members: List of compact IDs belonging to this family
        id_manager: CompactIDManager for sequence access
        logger: Optional logger instance

    Returns:
        Dictionary with sequence statistics
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Collect sequences only for members of this specific family
    family_sequences = []

    for compact_id in family_members:
        feature = id_manager.get_feature_by_compact_id(compact_id)
        if feature and hasattr(feature, "sequence") and feature.sequence:
            family_sequences.append(feature.sequence)

    if not family_sequences:
        return {
            "rep_length": "N/A",
            "min_length": "N/A",
            "max_length": "N/A",
            "mean_length": "N/A",
            "diversity": "N/A",
            "variance": "N/A",
        }

    # Calculate basic statistics
    lengths = [len(seq) for seq in family_sequences]

    if not lengths:
        return {
            "rep_length": "N/A",
            "min_length": "N/A",
            "max_length": "N/A",
            "mean_length": "N/A",
            "diversity": "N/A",
            "variance": "N/A",
        }

    # Get representative sequence length (use first sequence as representative)
    rep_length = lengths[0]

    # Calculate diversity (unique sequences / total sequences)
    unique_sequences = len(set(family_sequences))
    diversity = unique_sequences / len(family_sequences)

    # Calculate variance
    if len(lengths) > 1:
        mean_length = sum(lengths) / len(lengths)
        variance = sum((x - mean_length) ** 2 for x in lengths) / (len(lengths) - 1)
    else:
        mean_length = lengths[0]
        variance = 0

    return {
        "rep_length": rep_length,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": f"{mean_length:.1f}",
        "diversity": f"{diversity:.3f}",
        "variance": f"{variance:.1f}",
    }
