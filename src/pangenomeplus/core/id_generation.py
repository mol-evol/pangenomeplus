"""Functions for generating standardized PanGenomePlus IDs."""

import re
from pangenomeplus.core.types import FeatureType


def sanitize_id(original_id: str) -> str:
    """
    Sanitize ID by removing invalid characters and ensuring safe format.

    Args:
        original_id: Original identifier string

    Returns:
        Sanitized identifier safe for use in file systems and databases
    """
    # Replace problematic characters with underscores
    sanitized = re.sub(r'[^a-zA-Z0-9_.-]', '_', original_id)
    # Remove consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    # Ensure it starts with alphanumeric
    if sanitized and not sanitized[0].isalnum():
        sanitized = 'id_' + sanitized
    return sanitized or 'unknown'


def generate_pangenomeplus_id(
    genome_id: str,
    feature_type: FeatureType,
    original_id: str,
    contig_id: str,
    start: int,
    end: int
) -> str:
    """
    Generate standardized PanGenomePlus feature ID.

    Args:
        genome_id: Sanitized genome identifier
        feature_type: Type of genomic feature
        original_id: Original feature ID from annotation
        contig_id: Contig/scaffold identifier
        start: Start coordinate (1-based)
        end: End coordinate (1-based)

    Returns:
        Standardized PanGenomePlus feature ID

    Example:
        >>> generate_pangenomeplus_id(
        ...     "genome1", FeatureType.PROTEIN_CODING_GENE,
        ...     "gene_001", "contig_1", 1000, 2500
        ... )
        'PGP_genome1_protein_coding_gene_gene001_contig1_1000_2500'
    """
    # Sanitize all components
    genome_clean = sanitize_id(genome_id)
    feature_clean = feature_type.value
    original_clean = sanitize_id(original_id)
    contig_clean = sanitize_id(contig_id)

    # Build the ID
    pgp_id = f"PGP_{genome_clean}_{feature_clean}_{original_clean}_{contig_clean}_{start}_{end}"

    return pgp_id