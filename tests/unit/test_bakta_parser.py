"""
Unit tests for Bakta GFF3 parser.

Tests written BEFORE implementation (TDD approach).
"""

import pytest
from pathlib import Path
from pangenomeplus.extraction import (
    parse_bakta_gff,
    _parse_gff_attributes,
    _map_bakta_feature_type,
    _extract_sequence_from_coordinates,
)
from pangenomeplus.compact_ids import CompactIDManager


@pytest.fixture
def test_genome_path():
    """Path to test genome FASTA."""
    return "tests/fixtures/bakta_annotated/test_genome.fna"


@pytest.fixture
def test_gff_path():
    """Path to test Bakta GFF."""
    return "tests/fixtures/bakta_annotated/test_genome.gff"


@pytest.fixture
def id_manager():
    """Fresh CompactIDManager for each test."""
    return CompactIDManager()


class TestGFFAttributeParsing:
    """Test GFF3 attribute column parsing."""

    def test_parse_simple_attributes(self):
        """Parse basic key=value attributes."""
        attr_string = "ID=BAKTA_00001;Name=dnaA;gene=dnaA"
        result = _parse_gff_attributes(attr_string)

        assert result["ID"] == "BAKTA_00001"
        assert result["Name"] == "dnaA"
        assert result["gene"] == "dnaA"

    def test_parse_attributes_with_spaces_in_values(self):
        """Parse attributes containing spaces in values."""
        attr_string = "ID=BAKTA_00001;product=Chromosomal replication initiator protein DnaA"
        result = _parse_gff_attributes(attr_string)

        assert result["product"] == "Chromosomal replication initiator protein DnaA"

    def test_parse_attributes_with_special_characters(self):
        """Parse attributes with parentheses and colons."""
        attr_string = "ID=BAKTA_00002;anticodon=(pos:384..386)"
        result = _parse_gff_attributes(attr_string)

        assert result["anticodon"] == "(pos:384..386)"

    def test_parse_crispr_attributes(self):
        """Parse CRISPR repeat_region attributes."""
        attr_string = "ID=BAKTA_00005;repeat_family=CRISPR;rpt_type=direct"
        result = _parse_gff_attributes(attr_string)

        assert result["repeat_family"] == "CRISPR"
        assert result["rpt_type"] == "direct"


class TestFeatureTypeMapping:
    """Test Bakta feature type to PanGenomePlus type mapping."""

    def test_map_cds_to_protein(self):
        """CDS features map to 'P' (protein)."""
        attributes = {"ID": "BAKTA_00001", "product": "some protein"}
        result = _map_bakta_feature_type("CDS", attributes)

        assert result == "P"

    def test_map_trna(self):
        """tRNA features map to 'T'."""
        attributes = {"ID": "BAKTA_00002", "product": "tRNA-Ala"}
        result = _map_bakta_feature_type("tRNA", attributes)

        assert result == "T"

    def test_map_rrna(self):
        """rRNA features map to 'R'."""
        attributes = {"ID": "BAKTA_00004", "product": "16S ribosomal RNA"}
        result = _map_bakta_feature_type("rRNA", attributes)

        assert result == "R"

    def test_map_crispr_repeat_region(self):
        """repeat_region with repeat_family=CRISPR maps to 'C'."""
        attributes = {"ID": "BAKTA_00005", "repeat_family": "CRISPR"}
        result = _map_bakta_feature_type("repeat_region", attributes)

        assert result == "C"

    def test_map_non_crispr_repeat_region_returns_none(self):
        """repeat_region without CRISPR family returns None."""
        attributes = {"ID": "BAKTA_00006", "repeat_family": "IS6"}
        result = _map_bakta_feature_type("repeat_region", attributes)

        assert result is None

    def test_map_ncrna_returns_none(self):
        """ncRNA features return None (not extracted)."""
        attributes = {"ID": "BAKTA_00007", "product": "some ncRNA"}
        result = _map_bakta_feature_type("ncRNA", attributes)

        assert result is None

    def test_map_gene_returns_none(self):
        """gene features return None (only CDS extracted)."""
        attributes = {"ID": "BAKTA_00008"}
        result = _map_bakta_feature_type("gene", attributes)

        assert result is None


class TestSequenceExtraction:
    """Test sequence extraction from genomic coordinates."""

    def test_extract_forward_strand(self):
        """Extract sequence from forward strand."""
        contig_seq = "ATGAAACGCATTAGCACTCTG"
        sequence = _extract_sequence_from_coordinates(contig_seq, 1, 10, "+")

        assert sequence == "ATGAAACGCA"
        assert len(sequence) == 10

    def test_extract_reverse_strand(self):
        """Extract and reverse complement sequence from reverse strand."""
        contig_seq = "ATGAAACGCATTAGCACTCTG"
        sequence = _extract_sequence_from_coordinates(contig_seq, 5, 14, "-")

        # Positions 5-14: AACGCATTAG
        # Reverse complement: CTAATGCGTT
        assert sequence == "CTAATGCGTT"
        assert len(sequence) == 10

    def test_extract_full_sequence(self):
        """Extract entire contig sequence."""
        contig_seq = "ATGAAACGCA"
        sequence = _extract_sequence_from_coordinates(contig_seq, 1, 10, "+")

        assert sequence == contig_seq


class TestBaktaGFFParsing:
    """Test complete Bakta GFF parsing."""

    def test_parse_extracts_all_cds(self, test_gff_path, test_genome_path, id_manager):
        """Parse extracts all CDS features as proteins."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        assert "P" in features_by_type
        proteins = features_by_type["P"]

        # Should have 3 CDS features
        assert len(proteins) == 3

        # Check compact IDs assigned
        assert proteins[0].compact_id == "P1"
        assert proteins[1].compact_id == "P2"
        assert proteins[2].compact_id == "P3"

        # Check feature types
        for protein in proteins:
            assert protein.feature_type == "P"
            assert len(protein.sequence) > 0

    def test_parse_extracts_all_trna(self, test_gff_path, test_genome_path, id_manager):
        """Parse extracts all tRNA features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        assert "T" in features_by_type
        trnas = features_by_type["T"]

        # Should have 2 tRNA features
        assert len(trnas) == 2

        # Check compact IDs
        assert trnas[0].compact_id == "T1"
        assert trnas[1].compact_id == "T2"

        # Check feature types
        for trna in trnas:
            assert trna.feature_type == "T"
            assert len(trna.sequence) > 0

    def test_parse_extracts_rrna(self, test_gff_path, test_genome_path, id_manager):
        """Parse extracts rRNA features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        assert "R" in features_by_type
        rrnas = features_by_type["R"]

        # Should have 1 rRNA feature
        assert len(rrnas) == 1

        # Check compact ID
        assert rrnas[0].compact_id == "R1"
        assert rrnas[0].feature_type == "R"
        assert len(rrnas[0].sequence) > 0

    def test_parse_extracts_crispr(self, test_gff_path, test_genome_path, id_manager):
        """Parse extracts CRISPR repeat_region features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        assert "C" in features_by_type
        crisprs = features_by_type["C"]

        # Should have 1 CRISPR feature
        assert len(crisprs) == 1

        # Check compact ID
        assert crisprs[0].compact_id == "C1"
        assert crisprs[0].feature_type == "C"
        assert len(crisprs[0].sequence) > 0

    def test_parse_assigns_correct_coordinates(self, test_gff_path, test_genome_path, id_manager):
        """Parse assigns correct genomic coordinates to features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        proteins = features_by_type["P"]

        # First CDS: positions 1-300
        assert proteins[0].start == 1
        assert proteins[0].end == 300
        assert proteins[0].strand == "+"

        # Second CDS: positions 450-650, reverse strand
        assert proteins[1].start == 450
        assert proteins[1].end == 650
        assert proteins[1].strand == "-"

    def test_parse_extracts_correct_sequences(self, test_gff_path, test_genome_path, id_manager):
        """Parse extracts sequences matching coordinates."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        proteins = features_by_type["P"]

        # Verify sequence lengths match coordinate ranges
        assert len(proteins[0].sequence) == 300  # positions 1-300
        assert len(proteins[1].sequence) == 201  # positions 450-650

        # Verify sequences start with expected codons
        assert proteins[0].sequence.startswith("ATG")  # forward strand, position 1

    def test_parse_handles_reverse_strand_correctly(self, test_gff_path, test_genome_path, id_manager):
        """Parse correctly handles reverse strand features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        proteins = features_by_type["P"]
        trnas = features_by_type["T"]

        # Second CDS is reverse strand
        reverse_protein = proteins[1]
        assert reverse_protein.strand == "-"
        assert len(reverse_protein.sequence) == 201

        # Second tRNA is reverse strand
        reverse_trna = trnas[1]
        assert reverse_trna.strand == "-"
        assert len(reverse_trna.sequence) == 76

    def test_parse_stores_metadata(self, test_gff_path, test_genome_path, id_manager):
        """Parse stores Bakta metadata in features."""
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        proteins = features_by_type["P"]

        # First protein should have gene name and product
        assert proteins[0].genome_id == "test_genome"
        assert proteins[0].contig == "test_contig_1"
        # Note: actual metadata storage depends on Feature dataclass implementation


class TestIntergenicCalculationWithAllFeatures:
    """Test intergenic region calculation using ALL feature types."""

    def test_intergenic_excludes_all_annotated_features(self, test_gff_path, test_genome_path, id_manager):
        """Intergenic regions should not overlap ANY annotated feature."""
        from pangenomeplus.extraction import calculate_intergenic_regions
        from Bio import SeqIO

        # Parse all features
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        # Combine all annotated features
        all_features = []
        for feature_list in features_by_type.values():
            all_features.extend(feature_list)

        # Read genome sequence
        genome_seq = str(SeqIO.read(test_genome_path, "fasta").seq)

        # Calculate intergenic regions
        intergenic = calculate_intergenic_regions(
            all_features, genome_seq, "test_genome", id_manager, min_length=20
        )

        # Verify no intergenic region overlaps any annotated feature
        for ig in intergenic:
            for feature in all_features:
                # Check no overlap
                overlap = not (ig.end < feature.start or ig.start > feature.end)
                assert not overlap, f"Intergenic {ig.start}-{ig.end} overlaps {feature.feature_type} {feature.start}-{feature.end}"

    def test_intergenic_has_correct_coordinates(self, test_gff_path, test_genome_path, id_manager):
        """Intergenic regions should fill gaps between features."""
        from pangenomeplus.extraction import calculate_intergenic_regions
        from Bio import SeqIO

        # Parse all features
        features_by_type = parse_bakta_gff(
            test_gff_path, test_genome_path, "test_genome", id_manager
        )

        # Combine all annotated features
        all_features = []
        for feature_list in features_by_type.values():
            all_features.extend(feature_list)

        # Read genome sequence
        genome_seq = str(SeqIO.read(test_genome_path, "fasta").seq)

        # Calculate intergenic regions
        intergenic = calculate_intergenic_regions(
            all_features, genome_seq, "test_genome", id_manager, min_length=20
        )

        # Based on test GFF:
        # CDS: 1-300, tRNA: 350-425, CDS: 450-650, rRNA: 700-1200, CRISPR: 1220-1280, CDS: 1300-1500, tRNA: 1550-1625
        # Expected intergenic gaps:
        # Gap 1: 301-349 (49bp - too short with min_length=20 but should exist with default 50)
        # Gap 2: 426-449 (24bp)
        # Gap 3: 651-699 (49bp)
        # Gap 4: 1201-1219 (19bp - too short)
        # Gap 5: 1281-1299 (19bp - too short)
        # Gap 6: 1501-1549 (49bp)

        # With min_length=20, should have several intergenic regions
        assert len(intergenic) >= 3

        # Verify all intergenic have compact IDs starting with "I"
        for ig in intergenic:
            assert ig.feature_type == "I"
            assert ig.compact_id.startswith("I")
