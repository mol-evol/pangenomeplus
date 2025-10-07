"""MMseqs2-based clustering system for PanGenomePlus.

This module implements sequence clustering using MMseqs2 with feature-type-specific
parameters and bidirectional search capabilities for intergenic regions and CRISPR
spacers.

Pipeline Flow:
Features with Compact IDs → FASTA Generation → MMseqs2 Clustering → Family Assignment

Key Features:
- Unidirectional clustering for proteins, tRNAs, rRNAs
- Bidirectional clustering for intergenic regions and CRISPR spacers
- Feature-type-specific parameters for optimal clustering
- Family naming with feature type prefixes (FAM_P1, FAM_I42, etc.)
"""

import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .compact_ids import CompactIDManager
from .constants import (
    CLOUD_GENOME_THRESHOLD,
    CORE_GENOME_THRESHOLD,
    ClusteringDefaults,
)
from .core import Feature


class ClusteringError(Exception):
    """Exception raised during clustering operations."""

    pass


@dataclass
class ClusteringParams:
    """Parameters for MMseqs2 clustering."""

    min_seq_id: float = 0.8
    coverage: float = 0.8
    coverage_mode: int = 1  # 0: query, 1: target, 2: bidirectional
    cluster_mode: int = 0  # 0: set-cover, 1: connected-component, 2: greedy
    sensitivity: float = ClusteringDefaults.SENSITIVITY
    max_seqs: int = ClusteringDefaults.MAX_SEQS


@dataclass
class FamilyStats:
    """Statistics for a gene family."""

    family_id: str
    feature_type: str
    member_count: int
    genome_count: int
    representative_id: str
    classification: str  # "core", "accessory", "cloud", "singleton"


def get_clustering_params(
    feature_type: str, overrides: Optional[Dict[str, float]] = None
) -> ClusteringParams:
    """Get optimized clustering parameters for feature type.

    Args:
        feature_type: Feature type ("P", "I", "T", "R", "C")
        overrides: Optional dict with custom parameter values
                   (e.g., {"identity": 0.9, "coverage": 0.85})

    Returns:
        ClusteringParams with optimized settings

    Raises:
        ClusteringError: If feature type is invalid
    """
    params_map = {
        "P": ClusteringParams(
            min_seq_id=ClusteringDefaults.PROTEIN["identity"],
            coverage=ClusteringDefaults.PROTEIN["coverage"],
        ),  # Proteins
        "I": ClusteringParams(
            min_seq_id=ClusteringDefaults.INTERGENIC["identity"],
            coverage=ClusteringDefaults.INTERGENIC["coverage"],
        ),  # Intergenic (relaxed)
        "T": ClusteringParams(
            min_seq_id=ClusteringDefaults.TRNA["identity"],
            coverage=ClusteringDefaults.TRNA["coverage"],
        ),  # tRNAs (strict)
        "R": ClusteringParams(
            min_seq_id=ClusteringDefaults.RRNA["identity"],
            coverage=ClusteringDefaults.RRNA["coverage"],
        ),  # rRNAs (very strict)
        "C": ClusteringParams(
            min_seq_id=ClusteringDefaults.CRISPR["identity"],
            coverage=ClusteringDefaults.CRISPR["coverage"],
        ),  # CRISPR spacers
    }

    if feature_type not in params_map:
        raise ClusteringError(f"Invalid feature type: {feature_type}")

    params = params_map[feature_type]

    # Apply overrides if provided
    if overrides:
        if "identity" in overrides:
            params.min_seq_id = overrides["identity"]
        if "coverage" in overrides:
            params.coverage = overrides["coverage"]
        if "sensitivity" in overrides:
            params.sensitivity = overrides["sensitivity"]
        if "max_seqs" in overrides:
            params.max_seqs = int(overrides["max_seqs"])

    return params


def create_feature_fasta(
    features: List[Feature], output_file: str, bidirectional: bool = False
) -> str:
    """Create FASTA file from features using compact IDs as headers.

    Args:
        features: List of Feature objects
        output_file: Path to output FASTA file
        bidirectional: If True, include reverse complement sequences

    Returns:
        Path to created FASTA file

    Raises:
        ClusteringError: If file creation fails
    """
    if not features:
        raise ClusteringError("No features provided for FASTA creation")

    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        with open(output_file, "w") as f:
            for feature in features:
                if not feature.sequence:
                    continue

                # Forward sequence
                f.write(f">{feature.compact_id}\n")
                f.write(f"{feature.sequence}\n")

                # Reverse complement for bidirectional clustering
                if bidirectional:
                    reverse_complement = _get_reverse_complement(feature.sequence)
                    f.write(f">{feature.compact_id}_RC\n")
                    f.write(f"{reverse_complement}\n")

        return output_file

    except IOError as e:
        raise ClusteringError(f"Failed to create FASTA file {output_file}: {e}")


def _get_reverse_complement(sequence: str) -> str:
    """Get reverse complement of DNA sequence.

    Args:
        sequence: DNA sequence

    Returns:
        Reverse complement sequence
    """
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    # Handle mixed case and unknown nucleotides
    complement_map.update(
        {
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
            "n": "n",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "r": "y",
            "y": "r",
            "s": "s",
            "w": "w",
            "k": "m",
            "m": "k",
            "b": "v",
            "v": "b",
            "d": "h",
            "h": "d",
        }
    )

    reverse_seq = sequence[::-1]
    complement = "".join(complement_map.get(base, base) for base in reverse_seq)

    return complement


def run_mmseqs2_clustering(
    fasta_file: str,
    output_dir: str,
    params: ClusteringParams,
    temp_dir: Optional[str] = None,
) -> str:
    """Run MMseqs2 clustering on FASTA file.

    Args:
        fasta_file: Input FASTA file
        output_dir: Output directory for results
        params: Clustering parameters
        temp_dir: Temporary directory (optional)

    Returns:
        Path to cluster results TSV file

    Raises:
        ClusteringError: If clustering fails
    """
    if not os.path.exists(fasta_file):
        raise ClusteringError(f"FASTA file not found: {fasta_file}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Create temp directory if not provided
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="mmseqs2_")

    try:
        # MMseqs2 database files
        query_db = os.path.join(output_dir, "query_db")
        cluster_db = os.path.join(output_dir, "cluster_db")
        cluster_tsv = os.path.join(output_dir, "cluster_results.tsv")

        # Step 1: Create database from FASTA
        cmd_createdb = ["mmseqs", "createdb", fasta_file, query_db]

        result = subprocess.run(
            cmd_createdb, capture_output=True, text=True, check=False
        )

        if result.returncode != 0:
            raise ClusteringError(f"mmseqs createdb failed: {result.stderr}")

        # Step 2: Cluster sequences
        cmd_cluster = [
            "mmseqs",
            "cluster",
            query_db,
            cluster_db,
            temp_dir,
            "--min-seq-id",
            str(params.min_seq_id),
            "-c",
            str(params.coverage),
            "--cov-mode",
            str(params.coverage_mode),
            "--cluster-mode",
            str(params.cluster_mode),
            "-s",
            str(params.sensitivity),
            "--max-seqs",
            str(params.max_seqs),
        ]

        result = subprocess.run(
            cmd_cluster, capture_output=True, text=True, check=False
        )

        if result.returncode != 0:
            raise ClusteringError(f"mmseqs cluster failed: {result.stderr}")

        # Step 3: Convert to TSV
        cmd_createtsv = [
            "mmseqs",
            "createtsv",
            query_db,
            query_db,
            cluster_db,
            cluster_tsv,
        ]

        result = subprocess.run(
            cmd_createtsv, capture_output=True, text=True, check=False
        )

        if result.returncode != 0:
            raise ClusteringError(f"mmseqs createtsv failed: {result.stderr}")

        return cluster_tsv

    except subprocess.SubprocessError as e:
        raise ClusteringError(f"MMseqs2 execution failed: {e}")
    finally:
        # Clean up temporary files
        if temp_dir and os.path.exists(temp_dir):
            import shutil

            shutil.rmtree(temp_dir, ignore_errors=True)


def parse_cluster_results(
    cluster_tsv: str, bidirectional: bool = False
) -> Dict[str, List[str]]:
    """Parse MMseqs2 cluster results into representative -> members mapping.

    Args:
        cluster_tsv: Path to MMseqs2 TSV output
        bidirectional: If True, handle reverse complement mappings

    Returns:
        Dictionary mapping representative ID to list of member IDs

    Raises:
        ClusteringError: If parsing fails
    """
    if not os.path.exists(cluster_tsv):
        raise ClusteringError(f"Cluster results file not found: {cluster_tsv}")

    clusters = defaultdict(set)

    try:
        with open(cluster_tsv, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 2:
                    continue

                representative = parts[0]
                member = parts[1]

                # Handle bidirectional clustering
                if bidirectional:
                    # Remove _RC suffix for reverse complement sequences
                    if representative.endswith("_RC"):
                        representative = representative[:-3]
                    if member.endswith("_RC"):
                        member = member[:-3]

                clusters[representative].add(member)

        # Convert sets to sorted lists for consistency
        return {rep: sorted(list(members)) for rep, members in clusters.items()}

    except IOError as e:
        raise ClusteringError(f"Failed to parse cluster results: {e}")


def assign_family_ids(
    clusters: Dict[str, List[str]], feature_type: str
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """Assign family IDs to clusters and handle singletons.

    Args:
        clusters: Dictionary mapping representative to member list
        feature_type: Feature type for family naming

    Returns:
        Tuple of (compact_id -> family_id, family_id -> member_list)

    Raises:
        ClusteringError: If feature type is invalid
    """
    if feature_type not in {"P", "I", "T", "R", "C"}:
        raise ClusteringError(f"Invalid feature type: {feature_type}")

    compact_to_family = {}
    family_to_members = {}
    family_counter = 1

    for representative, members in clusters.items():
        if len(members) > 1:
            # Multi-member family
            family_id = f"FAM_{feature_type}{family_counter}"

            for member in members:
                compact_to_family[member] = family_id

            family_to_members[family_id] = members
            family_counter += 1

        else:
            # Singleton
            singleton_id = f"SING_{members[0]}"
            compact_to_family[members[0]] = singleton_id
            family_to_members[singleton_id] = members

    return compact_to_family, family_to_members


def classify_families(
    family_to_members: Dict[str, List[str]],
    id_manager: CompactIDManager,
    total_genomes: int,
    core_threshold: float = CORE_GENOME_THRESHOLD,
    cloud_threshold: float = CLOUD_GENOME_THRESHOLD,
) -> Dict[str, FamilyStats]:
    """Classify families as core, accessory, cloud, or singleton.

    Args:
        family_to_members: Mapping of family ID to member list
        id_manager: CompactIDManager for genome lookups
        total_genomes: Total number of genomes
        core_threshold: Threshold for core families (default: 95%)
        cloud_threshold: Threshold for cloud families (default: 15%)

    Returns:
        Dictionary mapping family ID to FamilyStats
    """
    family_stats = {}

    for family_id, members in family_to_members.items():
        # Count genomes containing this family
        genomes_with_family = set()

        for member_id in members:
            feature = id_manager.get_feature_by_compact_id(member_id)
            if feature:
                genomes_with_family.add(feature.genome_id)

        genome_count = len(genomes_with_family)
        genome_fraction = genome_count / total_genomes if total_genomes > 0 else 0

        # Determine classification based on genome frequency
        # Singletons are included - they're cloud families (present in 1 genome)
        if genome_fraction >= core_threshold:
            classification = "core"
        elif genome_fraction >= cloud_threshold:
            classification = "accessory"
        else:
            classification = "cloud"

        # Extract feature type and representative
        if family_id.startswith("FAM_"):
            feature_type = family_id[4]  # FAM_P1 -> P
            representative_id = members[0]  # Use first member as representative
        elif family_id.startswith("SING_"):
            feature_type = family_id[5]  # SING_P1 -> P
            representative_id = members[0]
        else:
            feature_type = "Unknown"
            representative_id = members[0] if members else ""

        family_stats[family_id] = FamilyStats(
            family_id=family_id,
            feature_type=feature_type,
            member_count=len(members),
            genome_count=genome_count,
            representative_id=representative_id,
            classification=classification,
        )

    return family_stats


def cluster_features_by_type(
    features: List[Feature],
    feature_type: str,
    output_dir: str,
    id_manager: CompactIDManager,
    total_genomes: int,
    cloud_threshold: float = 0.15,
    core_threshold: float = 0.95,
    clustering_overrides: Optional[Dict[str, float]] = None,
) -> Tuple[Dict[str, str], Dict[str, FamilyStats]]:
    """Cluster features of specific type and assign families.

    Args:
        features: List of features to cluster
        feature_type: Feature type ("P", "I", "T", "R", "C")
        output_dir: Output directory for intermediate files
        id_manager: CompactIDManager instance
        total_genomes: Total number of genomes for classification
        cloud_threshold: Threshold for cloud families (default: 0.15)
        core_threshold: Threshold for core families (default: 0.95)
        clustering_overrides: Optional custom clustering parameters

    Returns:
        Tuple of (compact_id -> family_id, family_id -> FamilyStats)

    Raises:
        ClusteringError: If clustering fails
    """
    if not features:
        return {}, {}

    # Create feature-type specific output directory
    type_output_dir = os.path.join(output_dir, f"clustering_{feature_type}")
    os.makedirs(type_output_dir, exist_ok=True)

    # Determine if bidirectional clustering is needed
    bidirectional = feature_type in {"I", "C"}  # Intergenic and CRISPR

    # Get clustering parameters with optional overrides
    params = get_clustering_params(feature_type, overrides=clustering_overrides)

    # Create FASTA file
    fasta_file = os.path.join(type_output_dir, f"sequences_{feature_type}.fasta")
    create_feature_fasta(features, fasta_file, bidirectional=bidirectional)

    # Run clustering
    cluster_results = run_mmseqs2_clustering(fasta_file, type_output_dir, params)

    # Parse results
    clusters = parse_cluster_results(cluster_results, bidirectional=bidirectional)

    # Assign family IDs
    compact_to_family, family_to_members = assign_family_ids(clusters, feature_type)

    # Classify families
    family_stats = classify_families(
        family_to_members,
        id_manager,
        total_genomes,
        cloud_threshold=cloud_threshold,
        core_threshold=core_threshold,
    )

    return compact_to_family, family_stats
