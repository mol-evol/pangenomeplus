"""Core data types and structures for the pangenome pipeline."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Tuple, Any, Literal
from pathlib import Path
from enum import Enum


class FeatureType(Enum):
    """Types of genomic features that can be analyzed."""
    PROTEIN_CODING_GENE = "protein_coding_gene"
    TRNA = "tRNA"
    RRNA = "rRNA"
    CRISPR_ARRAY = "CRISPR_array"
    INTERGENIC_REGION = "intergenic_region"


@dataclass
class GenomeFeature:
    """Complete genome feature representation with validation."""
    pangenomeplus_id: str
    original_id: str
    genome_id: str
    feature_type: FeatureType
    contig_id: str
    start: int  # 1-based coordinates
    end: int    # 1-based coordinates
    strand: Literal['+', '-', '.']
    product: str
    phase: Optional[int] = None
    attributes: Dict[str, Any] = field(default_factory=dict)
    source_tool: str = "unknown"


def validate_genome_feature(feature: GenomeFeature) -> None:
    """
    Validate genome feature data integrity.

    Args:
        feature: GenomeFeature to validate

    Raises:
        ValueError: If validation fails
    """
    if feature.start < 1:
        raise ValueError(f"Start coordinate must be >= 1, got {feature.start}")
    if feature.end < feature.start:
        raise ValueError(f"End coordinate {feature.end} must be >= start {feature.start}")
    if feature.strand not in ['+', '-', '.']:
        raise ValueError(f"Invalid strand: {feature.strand}")
    if not feature.pangenomeplus_id.startswith("PGP_"):
        raise ValueError(f"Invalid PanGenomePlus ID format: {feature.pangenomeplus_id}")


@dataclass
class ValidationResult:
    """Result of validation operation."""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)