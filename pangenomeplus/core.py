"""Core data structures for PanGenomePlus."""

from dataclasses import dataclass
from typing import Any, Dict

from pangenomeplus.constants import FeatureType


@dataclass(frozen=True)
class Feature:
    """
    Represents a genomic feature with compact ID and metadata.

    This is the core data structure for all genomic features in
    PanGenomePlus,
    including protein-coding genes, intergenic regions, tRNAs, rRNAs,
    and CRISPR elements.

    The dataclass is frozen (immutable) to ensure data integrity throughout
    the analysis pipeline.
    """

    compact_id: str  # P1, I5, T42, R7, C1A (Base36 format)
    genome_id: str  # Genome identifier (e.g., "E_coli_001")
    contig: str  # Contig/chromosome identifier
    start: int  # Start coordinate (1-based)
    end: int  # End coordinate (inclusive)
    strand: str  # Strand ("+", "-", or ".")
    sequence: str  # DNA/protein sequence
    feature_type: str  # Feature type ("P", "I", "T", "R", "C")
    original_id: str  # Original ID from annotation tool
    metadata: Dict[str, Any]  # Additional metadata from annotation tools

    def __hash__(self) -> int:
        """
        Custom hash function that handles unhashable metadata dict.

        Only includes hashable fields in the hash to avoid TypeError.
        """
        return hash(
            (
                self.compact_id,
                self.genome_id,
                self.contig,
                self.start,
                self.end,
                self.strand,
                self.sequence,
                self.feature_type,
                self.original_id,
                # Convert metadata dict to sorted tuple of items for hashing
                tuple(sorted(self.metadata.items())) if self.metadata else (),
            )
        )

    def __post_init__(self) -> None:
        """Validate feature data after initialization."""
        # Basic validation - could be extended later
        if not self.compact_id:
            raise ValueError("compact_id cannot be empty")

        if not self.feature_type:
            raise ValueError("feature_type cannot be empty")

        # Validate feature type is one of the supported types
        if self.feature_type not in FeatureType.ALL_TYPES:
            raise ValueError(f"feature_type must be one of {FeatureType.ALL_TYPES}")

        # Validate strand
        valid_strands = {"+", "-", "."}
        if self.strand not in valid_strands:
            raise ValueError(f"strand must be one of {valid_strands}")

        # Validate coordinates
        if self.start < 1:
            raise ValueError("start coordinate must be >= 1")

        if self.end < self.start:
            raise ValueError("end coordinate must be >= start coordinate")
