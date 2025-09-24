"""Comprehensive tests for compact ID system."""

import pytest

from pangenomeplus.compact_ids import (
    CompactIDError,
    CompactIDManager,
    generate_base36_id,
    validate_compact_id,
)
from pangenomeplus.core import Feature


class TestBase36Encoding:
    """Test Base36 ID encoding functionality."""

    def test_base36_encoding_small_numbers(self) -> None:
        """Test Base36 encoding for small sequential numbers."""
        test_cases = [
            (1, "1"),
            (9, "9"),
            (10, "A"),
            (35, "Z"),
            (36, "10"),
            (37, "11"),
        ]

        for number, expected in test_cases:
            result = generate_base36_id(number)
            assert result == expected

    def test_base36_encoding_large_numbers(self) -> None:
        """Test Base36 encoding for large numbers."""
        test_cases = [
            (1296, "100"),  # 36^2 = 1296
            (46656, "1000"),  # 36^3 = 46656
            (1000000, "LFLS"),  # Large number example
        ]

        for number, expected in test_cases:
            result = generate_base36_id(number)
            assert result == expected

    def test_base36_encoding_edge_cases(self) -> None:
        """Test Base36 encoding edge cases."""
        # Zero should raise error or handle appropriately
        with pytest.raises((ValueError, CompactIDError)):
            generate_base36_id(0)

        # Negative numbers should raise error
        with pytest.raises((ValueError, CompactIDError)):
            generate_base36_id(-1)

    def test_base36_encoding_consistency(self) -> None:
        """Test that Base36 encoding is consistent and reversible."""
        test_numbers = [1, 100, 1000, 10000, 100000, 1000000]

        for number in test_numbers:
            encoded = generate_base36_id(number)
            # Verify it's a valid Base36 string
            assert all(c in "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ" for c in encoded)
            # Verify it can be decoded back
            decoded = int(encoded, 36)
            assert decoded == number


class TestCompactIDValidation:
    """Test compact ID format validation."""

    def test_valid_compact_id_formats(self) -> None:
        """Test that valid compact ID formats pass validation."""
        valid_ids = [
            "P1",
            "P2",
            "P10",
            "P36",
            "P100",
            "I5",
            "I1A",
            "I2B7K",
            "T42",
            "TAB",
            "T1Z",
            "R7",
            "RZ",
            "R10A",
            "C1",
            "C999",
            "CZZZZ",
        ]

        for compact_id in valid_ids:
            result = validate_compact_id(compact_id)
            assert result is True

    def test_invalid_compact_id_formats(self) -> None:
        """Test that invalid compact ID formats fail validation."""
        invalid_ids = [
            "",  # Empty string
            "1",  # No feature type prefix
            "PP1",  # Double prefix
            "X1",  # Invalid feature type
            "P",  # No ID number
            "p1",  # Lowercase prefix
            "P1a",  # Lowercase in number
            "P-1",  # Invalid character
            "P 1",  # Space character
            "P1.2",  # Decimal point
        ]

        for compact_id in invalid_ids:
            result = validate_compact_id(compact_id)
            assert result is False

    def test_compact_id_feature_type_extraction(self) -> None:
        """Test extracting feature type from compact ID."""
        test_cases = [
            ("P1", "P"),
            ("I5", "I"),
            ("T42", "T"),
            ("R7", "R"),
            ("C1A", "C"),
        ]

        for compact_id, expected_type in test_cases:
            assert compact_id[0] == expected_type

    def test_compact_id_number_extraction(self) -> None:
        """Test extracting Base36 number from compact ID."""
        test_cases = [
            ("P1", "1"),
            ("I5", "5"),
            ("T42", "42"),
            ("R7", "7"),
            ("C1A", "1A"),
            ("P2B7K", "2B7K"),
        ]

        for compact_id, expected_number in test_cases:
            number_part = compact_id[1:]
            assert number_part == expected_number
            # Should be valid Base36
            int(number_part, 36)  # Should not raise exception


class TestCompactIDManager:
    """Test CompactIDManager functionality."""

    def test_compact_id_manager_initialization(self) -> None:
        """Test CompactIDManager initialization."""
        manager = CompactIDManager()

        # Should start with empty mappings
        assert len(manager.compact_to_full) == 0
        assert len(manager.location_to_compact) == 0
        assert len(manager.genome_features) == 0

        # Should initialize counters for all feature types
        expected_types = ["P", "I", "T", "R", "C"]
        for feature_type in expected_types:
            assert manager._counters[feature_type] == 0

    def test_generate_new_compact_id(self) -> None:
        """Test generating new compact IDs."""
        manager = CompactIDManager()

        # Test sequential ID generation
        id1 = manager.generate_compact_id("P")
        assert id1 == "P1"

        id2 = manager.generate_compact_id("P")
        assert id2 == "P2"

        id3 = manager.generate_compact_id("I")
        assert id3 == "I1"

        id4 = manager.generate_compact_id("I")
        assert id4 == "I2"

        # Different feature types should have independent counters
        assert manager._counters["P"] == 2
        assert manager._counters["I"] == 2
        assert manager._counters["T"] == 0

    def test_register_feature_creates_mappings(self) -> None:
        """Test that registering a feature creates all necessary mappings."""
        manager = CompactIDManager()

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

        manager.register_feature(feature)

        # Forward mapping should exist
        assert "P1" in manager.compact_to_full
        assert manager.compact_to_full["P1"] == feature

        # Reverse mapping should exist
        location_key = ("E_coli_001", 47, 1450)
        assert location_key in manager.location_to_compact
        assert manager.location_to_compact[location_key] == "P1"

        # Genome features should be updated
        assert "E_coli_001" in manager.genome_features
        assert "P" in manager.genome_features["E_coli_001"]
        assert "P1" in manager.genome_features["E_coli_001"]["P"]

    def test_register_multiple_features_same_genome(self) -> None:
        """Test registering multiple features from the same genome."""
        manager = CompactIDManager()

        features = [
            Feature(
                "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                "P2", "E_coli_001", "contig1", 300, 400, "+", "GCTA", "P", "gene2", {}
            ),
            Feature(
                "I1", "E_coli_001", "contig1", 250, 280, ".", "AAAA", "I", "inter1", {}
            ),
            Feature(
                "T1", "E_coli_001", "contig1", 500, 580, "+", "TGCA", "T", "trna1", {}
            ),
        ]

        for feature in features:
            manager.register_feature(feature)

        # Check genome features organization
        genome_features = manager.genome_features["E_coli_001"]
        assert len(genome_features["P"]) == 2
        assert len(genome_features["I"]) == 1
        assert len(genome_features["T"]) == 1
        assert "R" not in genome_features or len(genome_features["R"]) == 0

        # Check all features are in mappings
        assert len(manager.compact_to_full) == 4
        assert len(manager.location_to_compact) == 4

    def test_register_features_multiple_genomes(self) -> None:
        """Test registering features from multiple genomes."""
        manager = CompactIDManager()

        features = [
            Feature(
                "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                "P2", "E_coli_002", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                "P3", "E_coli_003", "contig1", 150, 250, "+", "GCTA", "P", "gene1", {}
            ),
        ]

        for feature in features:
            manager.register_feature(feature)

        # Should have features for all genomes
        assert "E_coli_001" in manager.genome_features
        assert "E_coli_002" in manager.genome_features
        assert "E_coli_003" in manager.genome_features

        # Each genome should have one protein
        for genome_id in ["E_coli_001", "E_coli_002", "E_coli_003"]:
            assert len(manager.genome_features[genome_id]["P"]) == 1

    def test_lookup_by_compact_id(self) -> None:
        """Test looking up features by compact ID."""
        manager = CompactIDManager()

        feature = Feature(
            "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
        )
        manager.register_feature(feature)

        # Successful lookup
        result = manager.get_feature_by_compact_id("P1")
        assert result == feature

        # Failed lookup
        result = manager.get_feature_by_compact_id("P999")
        assert result is None

    def test_lookup_by_location(self) -> None:
        """Test looking up features by genomic location."""
        manager = CompactIDManager()

        feature = Feature(
            "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
        )
        manager.register_feature(feature)

        # Successful lookup
        result = manager.get_compact_id_by_location("E_coli_001", 100, 200)
        assert result == "P1"

        # Failed lookup - wrong coordinates
        result = manager.get_compact_id_by_location("E_coli_001", 100, 199)
        assert result is None

        # Failed lookup - wrong genome
        result = manager.get_compact_id_by_location("E_coli_002", 100, 200)
        assert result is None

    def test_get_genome_features(self) -> None:
        """Test retrieving all features for a genome."""
        manager = CompactIDManager()

        features = [
            Feature(
                "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                "I1", "E_coli_001", "contig1", 250, 280, ".", "AAAA", "I", "inter1", {}
            ),
            Feature(
                "T1", "E_coli_001", "contig1", 500, 580, "+", "TGCA", "T", "trna1", {}
            ),
        ]

        for feature in features:
            manager.register_feature(feature)

        # Get features for existing genome
        genome_features = manager.get_features_for_genome("E_coli_001")
        assert genome_features is not None
        assert len(genome_features["P"]) == 1
        assert len(genome_features["I"]) == 1
        assert len(genome_features["T"]) == 1

        # Get features for non-existing genome
        genome_features = manager.get_features_for_genome("E_coli_999")
        assert genome_features is None

    def test_duplicate_compact_id_handling(self) -> None:
        """Test handling of duplicate compact IDs."""
        manager = CompactIDManager()

        feature1 = Feature(
            "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
        )
        manager.register_feature(feature1)

        # Attempt to register different feature with same compact ID
        feature2 = Feature(
            "P1", "E_coli_002", "contig1", 300, 400, "+", "GCTA", "P", "gene2", {}
        )

        with pytest.raises(CompactIDError):
            manager.register_feature(feature2)

    def test_duplicate_location_handling(self) -> None:
        """Test handling of features at the same genomic location."""
        manager = CompactIDManager()

        feature1 = Feature(
            "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
        )
        manager.register_feature(feature1)

        # Attempt to register different feature at same location
        feature2 = Feature(
            "P2", "E_coli_001", "contig1", 100, 200, "+", "GCTA", "P", "gene2", {}
        )

        with pytest.raises(CompactIDError):
            manager.register_feature(feature2)

    def test_counter_synchronization(self) -> None:
        """Test that counters stay synchronized with registered features."""
        manager = CompactIDManager()

        # Generate some IDs
        id1 = manager.generate_compact_id("P")  # P1
        id2 = manager.generate_compact_id("P")  # P2
        id3 = manager.generate_compact_id("I")  # I1

        assert manager._counters["P"] == 2
        assert manager._counters["I"] == 1

        # Register features with those IDs
        features = [
            Feature(
                id1, "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                id2, "E_coli_001", "contig1", 300, 400, "+", "GCTA", "P", "gene2", {}
            ),
            Feature(
                id3, "E_coli_001", "contig1", 250, 280, ".", "AAAA", "I", "inter1", {}
            ),
        ]

        for feature in features:
            manager.register_feature(feature)

        # Counters should remain synchronized
        assert manager._counters["P"] == 2
        assert manager._counters["I"] == 1

        # Next generated ID should continue sequence
        next_id = manager.generate_compact_id("P")
        assert next_id == "P3"


class TestCompactIDManagerScalability:
    """Test CompactIDManager scalability and performance characteristics."""

    def test_large_number_of_features(self) -> None:
        """Test manager performance with large numbers of features."""
        manager = CompactIDManager()

        # Create a large number of features (simulate real dataset)
        num_features = 10000
        features = []

        for i in range(1, num_features + 1):
            compact_id = manager.generate_compact_id("P")
            feature = Feature(
                compact_id=compact_id,
                genome_id=f"genome_{i % 100}",  # 100 genomes
                contig="contig1",
                start=i * 100,
                end=(i * 100) + 99,
                strand="+",
                sequence="ATCG",
                feature_type="P",
                original_id=f"gene_{i}",
                metadata={},
            )
            features.append(feature)

        # Register all features
        for feature in features:
            manager.register_feature(feature)

        # Verify all mappings are correct
        assert len(manager.compact_to_full) == num_features
        assert len(manager.location_to_compact) == num_features
        assert len(manager.genome_features) == 100  # 100 genomes

        # Verify lookups work correctly
        test_feature = features[5000]  # Middle feature
        lookup_result = manager.get_feature_by_compact_id(test_feature.compact_id)
        assert lookup_result == test_feature

    def test_base36_id_space_capacity(self) -> None:
        """Test Base36 ID space capacity for very large numbers."""
        # Test that we can generate IDs for very large counters
        large_numbers = [1000000, 10000000, 100000000, 2176782335]  # Up to ZZZ999

        for number in large_numbers:
            base36_id = generate_base36_id(number)
            assert len(base36_id) <= 6  # Should fit in reasonable string length
            # Verify it's valid Base36
            decoded = int(base36_id, 36)
            assert decoded == number

    def test_memory_efficiency_validation(self) -> None:
        """Test that data structures are memory efficient."""
        manager = CompactIDManager()

        # Create features with various sizes of metadata
        features = [
            Feature(
                "P1", "E_coli_001", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {}
            ),
            Feature(
                "P2",
                "E_coli_001",
                "contig1",
                300,
                400,
                "+",
                "GCTA",
                "P",
                "gene2",
                {
                    "product": "Very long protein name with detailed description",
                    "score": 95.5,
                    "domains": ["domain1", "domain2", "domain3"],
                },
            ),
        ]

        for feature in features:
            manager.register_feature(feature)

        # Verify manager doesn't duplicate large data structures
        # Features should be stored by reference, not copied
        stored_feature = manager.compact_to_full["P2"]
        assert stored_feature is features[1]  # Same object reference


class TestCompactIDManagerEdgeCases:
    """Test CompactIDManager edge cases and error conditions."""

    def test_empty_genome_id_handling(self) -> None:
        """Test handling of features with empty genome IDs."""
        manager = CompactIDManager()

        feature = Feature("P1", "", "contig1", 100, 200, "+", "ATCG", "P", "gene1", {})

        # Should handle empty genome ID gracefully
        manager.register_feature(feature)
        assert "" in manager.genome_features

    def test_special_characters_in_ids(self) -> None:
        """Test handling of special characters in various ID fields."""
        manager = CompactIDManager()

        feature = Feature(
            compact_id="P1",
            genome_id="E.coli-K12_MG1655",
            contig="NZ_CP007265.1",
            start=100,
            end=200,
            strand="+",
            sequence="ATCG",
            feature_type="P",
            original_id="locus_tag_001",
            metadata={},
        )

        manager.register_feature(feature)

        # Should handle special characters in genome ID
        assert "E.coli-K12_MG1655" in manager.genome_features

    def test_zero_length_features(self) -> None:
        """Test handling of zero-length features."""
        manager = CompactIDManager()

        # Feature with same start and end (zero length)
        feature = Feature(
            "P1", "E_coli_001", "contig1", 100, 100, "+", "", "P", "gene1", {}
        )

        manager.register_feature(feature)

        # Should be stored normally
        result = manager.get_feature_by_compact_id("P1")
        assert result == feature
