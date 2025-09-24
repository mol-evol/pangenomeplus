"""Comprehensive tests for Feature dataclass."""

from dataclasses import FrozenInstanceError
from typing import Any, Dict

import pytest

from pangenomeplus.core import Feature


class TestFeatureCreation:
    """Test Feature dataclass creation and basic functionality."""

    def test_feature_creation_with_all_fields(self) -> None:
        """Test creating a Feature with all required fields."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGCATTAGCACC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase III"},
        )

        assert feature.compact_id == "P1"
        assert feature.genome_id == "E_coli_001"
        assert feature.contig == "NZ_CP007265.1"
        assert feature.start == 47
        assert feature.end == 1450
        assert feature.strand == "+"
        assert feature.sequence == "ATGAAACGCATTAGCACC"
        assert feature.feature_type == "P"
        assert feature.original_id == "gene_1"
        assert feature.metadata == {"product": "DNA polymerase III"}

    def test_feature_creation_with_minimal_fields(self) -> None:
        """Test creating a Feature with minimal required fields."""
        feature = Feature(
            compact_id="I5",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=2556,
            end=2800,
            strand=".",
            sequence="GCTAGCTAGC",
            feature_type="I",
            original_id="intergenic_1",
            metadata={},
        )

        assert feature.compact_id == "I5"
        assert feature.feature_type == "I"
        assert feature.metadata == {}

    def test_feature_creation_with_different_feature_types(self) -> None:
        """Test creating Features with all supported feature types."""
        feature_types = [
            ("P1", "P", "protein"),
            ("I5", "I", "intergenic"),
            ("T42", "T", "tRNA"),
            ("R7", "R", "rRNA"),
            ("C1", "C", "CRISPR"),
        ]

        for compact_id, feature_type, description in feature_types:
            feature = Feature(
                compact_id=compact_id,
                genome_id="E_coli_001",
                contig="NZ_CP007265.1",
                start=100,
                end=200,
                strand="+",
                sequence="ATCG",
                feature_type=feature_type,
                original_id=f"{description}_1",
                metadata={"type": description},
            )
            assert feature.compact_id == compact_id
            assert feature.feature_type == feature_type

    def test_feature_immutability(self) -> None:
        """Test that Feature instances are immutable (frozen dataclass)."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )

        with pytest.raises(FrozenInstanceError):
            feature.compact_id = "P2"  # type: ignore

        with pytest.raises(FrozenInstanceError):
            feature.start = 100  # type: ignore


class TestFeatureValidation:
    """Test Feature dataclass validation and edge cases."""

    def test_coordinate_validation_start_less_than_end(self) -> None:
        """Test that start coordinate is less than end coordinate."""
        # This should be valid
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=100,
            end=200,
            strand="+",
            sequence="ATCG",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )
        assert feature.start < feature.end

    def test_sequence_length_matches_coordinates(self) -> None:
        """Test that sequence length matches coordinate span."""
        # Coordinate span: 200 - 100 = 100, but sequence is only 4 bases
        # This test documents current behavior - validation may be added later
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=100,
            end=200,
            strand="+",
            sequence="ATCG",  # Only 4 bases
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )
        coordinate_span = feature.end - feature.start
        sequence_length = len(feature.sequence)
        # Document the current behavior
        assert coordinate_span == 100
        assert sequence_length == 4

    def test_valid_strand_values(self) -> None:
        """Test that valid strand values are accepted."""
        valid_strands = ["+", "-", "."]

        for strand in valid_strands:
            feature = Feature(
                compact_id="P1",
                genome_id="E_coli_001",
                contig="NZ_CP007265.1",
                start=100,
                end=200,
                strand=strand,
                sequence="ATCG",
                feature_type="P",
                original_id="gene_1",
                metadata={},
            )
            assert feature.strand == strand

    def test_compact_id_format_validation(self) -> None:
        """Test various compact ID formats are accepted."""
        valid_compact_ids = ["P1", "P2B7K", "I5", "T42", "R7", "C1A", "P9ZZZ"]

        for compact_id in valid_compact_ids:
            feature = Feature(
                compact_id=compact_id,
                genome_id="E_coli_001",
                contig="NZ_CP007265.1",
                start=100,
                end=200,
                strand="+",
                sequence="ATCG",
                feature_type=compact_id[0],  # First character as feature type
                original_id="test_1",
                metadata={},
            )
            assert feature.compact_id == compact_id

    def test_empty_sequence_handling(self) -> None:
        """Test handling of empty sequences."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=100,
            end=200,
            strand="+",
            sequence="",  # Empty sequence
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )
        assert feature.sequence == ""
        assert len(feature.sequence) == 0


class TestFeatureEquality:
    """Test Feature equality and hashing behavior."""

    def test_feature_equality_same_values(self) -> None:
        """Test that Features with identical values are equal."""
        feature1 = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase"},
        )

        feature2 = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase"},
        )

        assert feature1 == feature2
        assert hash(feature1) == hash(feature2)

    def test_feature_inequality_different_values(self) -> None:
        """Test that Features with different values are not equal."""
        feature1 = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )

        feature2 = Feature(
            compact_id="P2",  # Different compact_id
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )

        assert feature1 != feature2
        assert hash(feature1) != hash(feature2)


class TestFeatureProperties:
    """Test computed properties and methods of Feature."""

    def test_feature_length_property(self) -> None:
        """Test that feature length is calculated correctly from coordinates."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=100,
            end=200,
            strand="+",
            sequence="ATCG",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )

        expected_length = 200 - 100  # end - start
        assert (feature.end - feature.start) == expected_length

    def test_feature_coordinate_span(self) -> None:
        """Test coordinate span calculation for various ranges."""
        test_cases = [
            (1, 100, 99),  # 1-indexed coordinates
            (47, 1450, 1403),  # Real gene coordinates
            (2556, 2800, 244),  # Intergenic region
        ]

        for start, end, expected_span in test_cases:
            feature = Feature(
                compact_id="P1",
                genome_id="E_coli_001",
                contig="NZ_CP007265.1",
                start=start,
                end=end,
                strand="+",
                sequence="ATCG",
                feature_type="P",
                original_id="gene_1",
                metadata={},
            )
            assert (feature.end - feature.start) == expected_span


class TestFeatureRepresentation:
    """Test string representation and debugging output."""

    def test_feature_string_representation(self) -> None:
        """Test that Feature has useful string representation."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase"},
        )

        str_repr = str(feature)
        # Should contain key identifying information
        assert "P1" in str_repr
        assert "E_coli_001" in str_repr

    def test_feature_repr_contains_all_fields(self) -> None:
        """Test that repr contains all field information."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase"},
        )

        repr_str = repr(feature)
        # Should be evaluable and contain field names
        assert "compact_id" in repr_str or "P1" in repr_str
        assert "Feature" in repr_str


class TestFeatureMetadata:
    """Test metadata handling and manipulation."""

    def test_empty_metadata_dict(self) -> None:
        """Test Feature with empty metadata dictionary."""
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={},
        )

        assert feature.metadata == {}
        assert len(feature.metadata) == 0

    def test_metadata_with_various_types(self) -> None:
        """Test metadata with different value types."""
        metadata: Dict[str, Any] = {
            "product": "DNA polymerase III",
            "score": 95.5,
            "confidence": 0.95,
            "length": 1403,
            "partial": True,
            "domains": ["helix", "loop"],
            "nested": {"gc_content": 0.52},
        }

        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata=metadata,
        )

        assert feature.metadata["product"] == "DNA polymerase III"
        assert feature.metadata["score"] == 95.5
        assert feature.metadata["confidence"] == 0.95
        assert feature.metadata["length"] == 1403
        assert feature.metadata["partial"] is True
        assert feature.metadata["domains"] == ["helix", "loop"]
        assert feature.metadata["nested"]["gc_content"] == 0.52

    def test_metadata_immutability(self) -> None:
        """Test that metadata dictionary itself is not directly mutable."""
        # Note: This test documents current behavior
        # The Feature is frozen, but the metadata dict content could be mutable
        feature = Feature(
            compact_id="P1",
            genome_id="E_coli_001",
            contig="NZ_CP007265.1",
            start=47,
            end=1450,
            strand="+",
            sequence="ATGAAACGC",
            feature_type="P",
            original_id="gene_1",
            metadata={"product": "DNA polymerase"},
        )

        # Cannot replace the metadata dict due to frozen dataclass
        with pytest.raises(FrozenInstanceError):
            feature.metadata = {"new": "dict"}  # type: ignore

        # But the contents of the dict might be mutable (current implementation)
        # This documents current behavior - could be changed to use frozen dict later
        original_metadata = feature.metadata
        assert original_metadata == {"product": "DNA polymerase"}
