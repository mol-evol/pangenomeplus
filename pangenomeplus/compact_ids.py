"""Compact ID management system for PanGenomePlus."""

import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .core import Feature


class CompactIDError(Exception):
    """Exception raised for compact ID management errors."""

    pass


def generate_base36_id(number: int) -> str:
    """
    Convert integer to Base36 string representation.

    Args:
        number: Positive integer to convert

    Returns:
        Base36 string representation

    Raises:
        ValueError: If number is <= 0
    """
    if number <= 0:
        raise ValueError("Number must be positive")

    if number == 0:
        return "0"

    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    result = ""

    while number > 0:
        result = digits[number % 36] + result
        number //= 36

    return result


def validate_compact_id(compact_id: str) -> bool:
    """
    Validate compact ID format.

    Valid format: {FEATURE_TYPE}{BASE36_NUMBER}
    where FEATURE_TYPE is one of P, I, T, R, C
    and BASE36_NUMBER contains only 0-9, A-Z

    Args:
        compact_id: Compact ID string to validate

    Returns:
        True if valid, False otherwise
    """
    if not compact_id:
        return False

    # Must be at least 2 characters (type + number)
    if len(compact_id) < 2:
        return False

    # First character must be valid feature type
    feature_type = compact_id[0]
    if feature_type not in {"P", "I", "T", "R", "C"}:
        return False

    # Second character must NOT be a feature type (avoid double prefix)
    if compact_id[1] in {"P", "I", "T", "R", "C"}:
        return False

    # Remaining characters must be valid Base36
    number_part = compact_id[1:]
    if not number_part:
        return False

    # Check if all characters in number part are valid Base36
    valid_base36_pattern = re.compile(r"^[0-9A-Z]+$")
    return bool(valid_base36_pattern.match(number_part))


class CompactIDManager:
    """
    Manages compact ID generation and bidirectional mapping.

    Provides universal Base36 ID assignment across all feature types
    with efficient lookup capabilities.
    """

    def __init__(self) -> None:
        """Initialize CompactIDManager with empty mappings."""
        # Forward mapping: compact_id -> Feature
        self.compact_to_full: Dict[str, Feature] = {}

        # Reverse mapping: (genome_id, contig, start, end) -> compact_id
        self.location_to_compact: Dict[Tuple[str, str, int, int], str] = {}

        # Genome-based mapping: genome_id -> {feature_type: [compact_ids]}
        self.genome_features: Dict[str, Dict[str, List[str]]] = defaultdict(
            lambda: defaultdict(list)
        )

        # Counters for each feature type
        self._counters: Dict[str, int] = {
            "P": 0,  # Proteins
            "I": 0,  # Intergenic regions
            "T": 0,  # tRNAs
            "R": 0,  # rRNAs
            "C": 0,  # CRISPR elements
        }

    def generate_compact_id(self, feature_type: str) -> str:
        """
        Generate next sequential compact ID for feature type.

        Args:
            feature_type: Feature type ("P", "I", "T", "R", "C")

        Returns:
            New compact ID (e.g., "P1", "I5", "T42")

        Raises:
            CompactIDError: If feature_type is invalid
        """
        if feature_type not in self._counters:
            raise CompactIDError(f"Invalid feature type: {feature_type}")

        # Increment counter and generate ID
        self._counters[feature_type] += 1
        number = self._counters[feature_type]

        base36_number = generate_base36_id(number)
        return f"{feature_type}{base36_number}"

    def register_feature(self, feature: Feature) -> None:
        """
        Register a feature in all mapping tables.

        Args:
            feature: Feature instance to register

        Raises:
            CompactIDError: If compact_id or location already exists
        """
        compact_id = feature.compact_id
        location_key = (feature.genome_id, feature.contig, feature.start, feature.end)

        # Check for duplicate compact ID
        if compact_id in self.compact_to_full:
            raise CompactIDError(f"Compact ID {compact_id} already exists")

        # Check for duplicate location
        if location_key in self.location_to_compact:
            existing_id = self.location_to_compact[location_key]
            raise CompactIDError(
                f"Location {location_key} already occupied by {existing_id}"
            )

        # Register in all mappings
        self.compact_to_full[compact_id] = feature
        self.location_to_compact[location_key] = compact_id
        self.genome_features[feature.genome_id][feature.feature_type].append(compact_id)

    def get_feature_by_compact_id(self, compact_id: str) -> Optional[Feature]:
        """
        Retrieve feature by compact ID.

        Args:
            compact_id: Compact ID to look up

        Returns:
            Feature instance or None if not found
        """
        return self.compact_to_full.get(compact_id)

    def get_compact_id_by_location(
        self, genome_id: str, contig: str, start: int, end: int
    ) -> Optional[str]:
        """
        Retrieve compact ID by genomic location.

        Args:
            genome_id: Genome identifier
            contig: Contig identifier
            start: Start coordinate
            end: End coordinate

        Returns:
            Compact ID or None if not found
        """
        location_key = (genome_id, contig, start, end)
        return self.location_to_compact.get(location_key)

    def get_features_for_genome(self, genome_id: str) -> Optional[Dict[str, List[str]]]:
        """
        Retrieve all features for a genome organized by type.

        Args:
            genome_id: Genome identifier

        Returns:
            Dictionary mapping feature type to list of compact IDs,
            or None if genome not found
        """
        if genome_id not in self.genome_features:
            return None

        # Return regular dict instead of defaultdict for cleaner interface
        return dict(self.genome_features[genome_id])
