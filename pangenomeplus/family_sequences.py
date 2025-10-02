#!/usr/bin/env python3
"""
Family sequence utilities for PanGenomePlus.

This module provides utilities for generating representative sequences
for gene families, handling sequence access optimization, and family
sequence analysis.
"""

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from .clustering import FamilyStats
from .compact_ids import CompactIDManager
from .core import Feature


def generate_family_representatives(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    selection_method: str = "first",
) -> Dict[str, str]:
    """Generate representative sequences for each gene family.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for sequence access
        selection_method: Method for selecting representatives ("first", "longest", "shortest")

    Returns:
        Dict mapping family_id to representative compact_id

    Raises:
        ValueError: If selection method is invalid
    """
    valid_methods = ["first", "longest", "shortest"]
    if selection_method not in valid_methods:
        raise ValueError(f"selection_method must be one of {valid_methods}")

    logger = logging.getLogger(__name__)
    logger.info(f"Generating family representatives using '{selection_method}' method")

    family_representatives = {}

    # Collect all members for each family
    family_members: Dict[str, List[str]] = {}
    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            if family_id not in family_members:
                family_members[family_id] = []
            family_members[family_id].append(compact_id)

    # Select representative for each family
    for family_id, members in family_members.items():
        if selection_method == "first":
            # Use first member (alphabetically sorted for consistency)
            representative = sorted(members)[0]

        elif selection_method == "longest":
            # Use member with longest sequence
            longest_id = None
            longest_length = -1

            for member_id in members:
                feature = id_manager.get_feature_by_compact_id(member_id)
                if feature and feature.sequence:
                    seq_length = len(feature.sequence)
                    if seq_length > longest_length:
                        longest_length = seq_length
                        longest_id = member_id

            representative = longest_id if longest_id else sorted(members)[0]

        elif selection_method == "shortest":
            # Use member with shortest sequence
            shortest_id = None
            shortest_length = float("inf")

            for member_id in members:
                feature = id_manager.get_feature_by_compact_id(member_id)
                if feature and feature.sequence:
                    seq_length = len(feature.sequence)
                    if seq_length < shortest_length:
                        shortest_length = seq_length
                        shortest_id = member_id

            representative = shortest_id if shortest_id else sorted(members)[0]

        family_representatives[family_id] = representative

    logger.info(f"Generated representatives for {len(family_representatives)} families")
    return family_representatives


def write_family_fasta(
    family_representatives: Dict[str, str],
    id_manager: CompactIDManager,
    output_file: str,
    include_metadata: bool = True,
) -> int:
    """Write family representative sequences to FASTA file.

    Args:
        family_representatives: Dict mapping family_id to representative compact_id
        id_manager: CompactIDManager for sequence access
        output_file: Path to output FASTA file
        include_metadata: Whether to include metadata in FASTA headers

    Returns:
        Number of sequences written

    Raises:
        IOError: If file writing fails
    """
    logger = logging.getLogger(__name__)
    sequences_written = 0

    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        for family_id, representative_id in sorted(family_representatives.items()):
            feature = id_manager.get_feature_by_compact_id(representative_id)

            if feature and feature.sequence:
                if include_metadata:
                    # Rich header with metadata
                    header = f">{family_id} representative={representative_id} genome={feature.genome_id} pos={feature.contig}:{feature.start}-{feature.end} strand={feature.strand}"
                else:
                    # Simple header
                    header = f">{family_id}"

                f.write(f"{header}\n")
                f.write(f"{feature.sequence}\n")
                sequences_written += 1
            else:
                logger.warning(
                    f"No sequence found for representative {representative_id} in family {family_id}"
                )

    logger.info(f"Wrote {sequences_written} representative sequences to {output_file}")
    return sequences_written


def analyze_family_sequence_stats(
    family_assignments: Dict[str, Dict[str, str]], id_manager: CompactIDManager
) -> Dict[str, Dict[str, Any]]:
    """Analyze sequence statistics for gene families.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for sequence access

    Returns:
        Dict mapping family_id to sequence statistics

    Statistics include:
        - member_count: Number of sequences in family
        - min_length, max_length, mean_length: Sequence length statistics
        - feature_type: Type of genomic feature
        - sequence_diversity: Basic diversity metrics
    """
    logger = logging.getLogger(__name__)
    logger.info("Analyzing family sequence statistics")

    family_stats = {}

    # Collect sequences for each family
    family_sequences = {}
    for feature_type, assignments in family_assignments.items():
        for compact_id, family_id in assignments.items():
            if family_id not in family_sequences:
                family_sequences[family_id] = {
                    "sequences": [],
                    "feature_type": feature_type,
                }

            feature = id_manager.get_feature_by_compact_id(compact_id)
            if feature and feature.sequence:
                seq_list: List[str] = family_sequences[family_id]["sequences"]  # type: ignore[assignment]
                seq_list.append(feature.sequence)

    # Calculate statistics for each family
    for family_id, data in family_sequences.items():
        sequences = data["sequences"]
        if not sequences:
            continue

        # Basic statistics
        lengths = [len(seq) for seq in sequences]
        member_count = len(sequences)

        # Sequence diversity (simple metric: unique sequences / total sequences)
        unique_sequences = len(set(sequences))
        diversity = unique_sequences / member_count if member_count > 0 else 0

        family_stats[family_id] = {
            "member_count": member_count,
            "feature_type": data["feature_type"],
            "min_length": min(lengths),
            "max_length": max(lengths),
            "mean_length": sum(lengths) / len(lengths),
            "unique_sequences": unique_sequences,
            "sequence_diversity": diversity,
            "length_variance": _calculate_variance(lengths),
        }

    logger.info(f"Analyzed sequence statistics for {len(family_stats)} families")
    return family_stats


def _calculate_variance(values: Sequence[float]) -> float:
    """Calculate variance of a list of values."""
    if len(values) <= 1:
        return 0.0

    mean = sum(values) / len(values)
    squared_diffs = [(x - mean) ** 2 for x in values]
    variance = sum(squared_diffs) / (len(values) - 1)
    return variance


def create_family_sequence_report(
    family_stats: Dict[str, Dict[str, Any]], output_file: str
) -> None:
    """Create a comprehensive report of family sequence statistics.

    Args:
        family_stats: Family statistics from analyze_family_sequence_stats()
        output_file: Path to output report file

    Raises:
        IOError: If file writing fails
    """
    logger = logging.getLogger(__name__)

    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        # Header
        f.write("# Family Sequence Analysis Report\n")
        f.write("# Generated by PanGenomePlus\n\n")

        # Summary statistics
        total_families = len(family_stats)
        feature_types = set(stats["feature_type"] for stats in family_stats.values())

        f.write("## Summary\n")
        f.write(f"Total families analyzed: {total_families}\n")
        f.write(f"Feature types: {', '.join(sorted(feature_types))}\n\n")

        # Per-feature-type summary
        f.write("## Feature Type Summary\n")
        for feature_type in sorted(feature_types):
            type_families = [
                fid
                for fid, stats in family_stats.items()
                if stats["feature_type"] == feature_type
            ]
            f.write(f"- {feature_type}: {len(type_families)} families\n")

        f.write("\n")

        # Detailed family statistics
        f.write("## Detailed Family Statistics\n")
        f.write(
            "Family_ID\tFeature_Type\tMembers\tMin_Length\tMax_Length\tMean_Length\tUnique_Seqs\tDiversity\tLength_Variance\n"
        )

        for family_id in sorted(family_stats.keys()):
            stats = family_stats[family_id]
            f.write(f"{family_id}\t{stats['feature_type']}\t{stats['member_count']}\t")
            f.write(
                f"{stats['min_length']}\t{stats['max_length']}\t{stats['mean_length']:.1f}\t"
            )
            f.write(f"{stats['unique_sequences']}\t{stats['sequence_diversity']:.3f}\t")
            f.write(f"{stats['length_variance']:.1f}\n")

    logger.info(f"Created family sequence report: {output_file}")


def generate_family_sequence_bundle(
    family_assignments: Dict[str, Dict[str, str]],
    id_manager: CompactIDManager,
    output_dir: str,
    selection_method: str = "first",
) -> Dict[str, str]:
    """Generate complete family sequence bundle with FASTA and reports.

    This is a convenience function that generates representative sequences,
    analyzes statistics, and creates comprehensive output files.

    Args:
        family_assignments: Dict mapping feature types to compact_id -> family_id
        id_manager: CompactIDManager for sequence access
        output_dir: Output directory for all files
        selection_method: Method for selecting representatives

    Returns:
        Dict mapping output type to file path

    Output files created:
        - representatives.fasta: Representative sequences
        - family_statistics.tsv: Detailed statistics report
        - sequence_summary.txt: Human-readable summary
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating family sequence bundle in {output_dir}")

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    output_files = {}

    # Generate representatives
    representatives = generate_family_representatives(
        family_assignments, id_manager, selection_method
    )

    # Write representative FASTA
    fasta_file = os.path.join(output_dir, "representatives.fasta")
    sequences_written = write_family_fasta(
        representatives, id_manager, fasta_file, include_metadata=True
    )
    output_files["fasta"] = fasta_file

    # Analyze sequence statistics
    sequence_stats = analyze_family_sequence_stats(family_assignments, id_manager)

    # Create detailed report
    report_file = os.path.join(output_dir, "family_statistics.tsv")
    create_family_sequence_report(sequence_stats, report_file)
    output_files["report"] = report_file

    # Create summary
    summary_file = os.path.join(output_dir, "sequence_summary.txt")
    with open(summary_file, "w") as f:
        f.write("# Family Sequence Bundle Summary\n")
        f.write("# Generated by PanGenomePlus\n\n")
        f.write(f"Representative selection method: {selection_method}\n")
        f.write(f"Total families: {len(representatives)}\n")
        f.write(f"Representative sequences written: {sequences_written}\n")
        f.write(f"Families with statistics: {len(sequence_stats)}\n\n")
        f.write("Files created:\n")
        for output_type, filepath in output_files.items():
            f.write(f"- {output_type}: {os.path.basename(filepath)}\n")

    output_files["summary"] = summary_file

    logger.info(f"Family sequence bundle complete: {len(output_files)} files created")
    return output_files
