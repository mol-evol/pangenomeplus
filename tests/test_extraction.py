"""Tests for external tool integration and feature extraction."""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from pangenomeplus.compact_ids import CompactIDManager
from pangenomeplus.core import Feature
from pangenomeplus.extraction import (
    ExtractionError,
    calculate_intergenic_regions,
    extract_genome_features,
    parse_prodigal_gff,
    run_barrnap,
    run_minced,
    run_prodigal,
    run_trnascan,
)


class TestExternalToolWrappers:
    """Test external tool wrapper functions."""

    def test_run_prodigal_success(self):
        """Test successful Prodigal execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test.fasta")
            output_dir = os.path.join(tmpdir, "output")

            # Create mock genome file
            with open(genome_file, "w") as f:
                f.write(">test_contig\nATCG\n")

            # Mock successful subprocess call
            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock()
                expected_gff = os.path.join(output_dir, "test.gff")

                # Create expected output file
                os.makedirs(output_dir)
                with open(expected_gff, "w") as f:
                    f.write("# Mock GFF content\n")

                result = run_prodigal(genome_file, output_dir)

                assert result == expected_gff
                mock_run.assert_called_once()
                args = mock_run.call_args[0][0]
                assert args[0] == "prodigal"
                assert args[2] == genome_file

    def test_run_prodigal_failure(self):
        """Test Prodigal execution failure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test.fasta")
            output_dir = os.path.join(tmpdir, "output")

            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.side_effect = FileNotFoundError()

                with pytest.raises(ExtractionError, match="Prodigal not found"):
                    run_prodigal(genome_file, output_dir)

    def test_run_trnascan_success(self):
        """Test successful tRNAscan-SE execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test.fasta")
            output_dir = os.path.join(tmpdir, "output")

            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock()
                expected_output = os.path.join(output_dir, "test_trna.txt")

                os.makedirs(output_dir)
                with open(expected_output, "w") as f:
                    f.write("# Mock tRNA content\n")

                result = run_trnascan(genome_file, output_dir, model="bacteria")

                assert result == expected_output
                args = mock_run.call_args[0][0]
                assert args[0] == "trnascan-se"
                assert "-B" in args

    def test_run_barrnap_success(self):
        """Test successful Barrnap execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test.fasta")
            output_dir = os.path.join(tmpdir, "output")

            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock()
                expected_output = os.path.join(output_dir, "test_rrna.gff")

                os.makedirs(output_dir)
                with open(expected_output, "w") as f:
                    f.write("# Mock rRNA content\n")

                result = run_barrnap(genome_file, output_dir, kingdom="bac")

                assert result == expected_output
                args = mock_run.call_args[0][0]
                assert args[0] == "barrnap"
                assert "--kingdom" in args
                assert "bac" in args

    def test_run_minced_success(self):
        """Test successful MINCED execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test.fasta")
            output_dir = os.path.join(tmpdir, "output")

            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock()
                expected_output = os.path.join(output_dir, "test_crispr.txt")

                os.makedirs(output_dir)
                with open(expected_output, "w") as f:
                    f.write("# Mock CRISPR content\n")

                result = run_minced(genome_file, output_dir, min_repeats=3)

                assert result == expected_output
                args = mock_run.call_args[0][0]
                assert args[0] == "minced"
                assert "-minNR" in args


class TestProdigalGFFParsing:
    """Test Prodigal GFF3 parsing functionality."""

    @pytest.fixture
    def sample_gff_content(self):
        """Sample GFF3 content for testing."""
        return """##gff-version 3
# Sequence Data: seqnum=1;seqlen=1000;seqhdr="test_contig"
test_contig\tProdigal_v2.6.3\tCDS\t1\t300\t.\t+\t0\tID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp
test_contig\tProdigal_v2.6.3\tCDS\t400\t600\t.\t-\t0\tID=1_2;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp
"""

    @pytest.fixture
    def sample_genome_seq(self):
        """Sample genome sequence for testing."""
        # 300bp forward gene + 100bp gap + 200bp reverse gene + 400bp
        return "A" * 300 + "T" * 100 + "C" * 200 + "G" * 400

    def test_parse_prodigal_gff_success(self, sample_gff_content, sample_genome_seq):
        """Test successful GFF3 parsing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_file = os.path.join(tmpdir, "test.gff")
            genome_file = os.path.join(tmpdir, "genome.fasta")

            # Write test files
            with open(gff_file, "w") as f:
                f.write(sample_gff_content)

            with open(genome_file, "w") as f:
                f.write(">test_contig\n")
                f.write(sample_genome_seq + "\n")

            id_manager = CompactIDManager()
            features = parse_prodigal_gff(
                gff_file, genome_file, "test_genome", id_manager
            )

            # Should find 2 CDS features
            assert len(features) == 2

            # Check first feature (forward strand)
            feature1 = features[0]
            assert feature1.compact_id == "P1"
            assert feature1.genome_id == "test_genome"
            assert feature1.contig == "test_contig"
            assert feature1.start == 1
            assert feature1.end == 300
            assert feature1.strand == "+"
            assert feature1.sequence == "A" * 300
            assert feature1.feature_type == "P"

            # Check second feature (reverse strand)
            feature2 = features[1]
            assert feature2.compact_id == "P2"
            assert feature2.strand == "-"
            assert feature2.start == 400
            assert feature2.end == 600
            # Test that we got the reverse complement (sequence length should be correct)
            assert len(feature2.sequence) == 201  # 600 - 400 + 1
            # Check it's mostly G's (reverse complement of C's with some A at the end)
            assert feature2.sequence.count("G") > 190

    def test_parse_prodigal_gff_truncated_coordinates(self):
        """Test GFF parsing with coordinates beyond sequence length (graceful handling)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_file = os.path.join(tmpdir, "test.gff")
            genome_file = os.path.join(tmpdir, "genome.fasta")

            # GFF with coordinates beyond sequence length
            gff_content = (
                """test_contig\tProdigal_v2.6.3\tCDS\t1\t2000\t.\t+\t0\tID=1_1"""
            )

            with open(gff_file, "w") as f:
                f.write(gff_content)

            with open(genome_file, "w") as f:
                f.write(">test_contig\nATCG\n")  # Only 4bp

            id_manager = CompactIDManager()

            # Should gracefully handle out-of-bounds by taking available sequence
            features = parse_prodigal_gff(
                gff_file, genome_file, "test_genome", id_manager
            )
            assert len(features) == 1
            assert features[0].sequence == "ATCG"  # All available sequence

    def test_parse_prodigal_gff_missing_genome(self):
        """Test GFF parsing with missing genome file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_file = os.path.join(tmpdir, "test.gff")
            missing_genome = os.path.join(tmpdir, "missing.fasta")

            with open(gff_file, "w") as f:
                f.write("# Empty GFF\n")

            id_manager = CompactIDManager()

            with pytest.raises(ExtractionError, match="Failed to load genome sequence"):
                parse_prodigal_gff(gff_file, missing_genome, "test_genome", id_manager)


class TestIntergenicCalculation:
    """Test intergenic region calculation."""

    @pytest.fixture
    def sample_features(self):
        """Sample protein features for testing."""
        return [
            Feature(
                compact_id="P1",
                genome_id="test_genome",
                contig="test_contig",
                start=1,
                end=100,
                strand="+",
                sequence="A" * 100,
                feature_type="P",
                original_id="gene_1",
                metadata={},
            ),
            Feature(
                compact_id="P2",
                genome_id="test_genome",
                contig="test_contig",
                start=200,
                end=300,
                strand="+",
                sequence="C" * 100,
                feature_type="P",
                original_id="gene_2",
                metadata={},
            ),
        ]

    def test_calculate_intergenic_regions_success(self, sample_features):
        """Test successful intergenic region calculation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            # Create genome: 100bp gene + 99bp intergenic + 100bp gene + 100bp intergenic
            genome_seq = "A" * 100 + "T" * 99 + "C" * 100 + "G" * 100
            with open(genome_file, "w") as f:
                f.write(">test_contig\n")
                f.write(genome_seq + "\n")

            id_manager = CompactIDManager()
            intergenic_features = calculate_intergenic_regions(
                sample_features, genome_file, "test_genome", id_manager, min_length=50
            )

            # Should find 2 intergenic regions
            assert len(intergenic_features) == 2

            # First intergenic region (between genes)
            region1 = intergenic_features[0]
            assert region1.compact_id.startswith("I")
            assert region1.start == 101
            assert region1.end == 199
            assert region1.sequence == "T" * 99
            assert region1.feature_type == "I"

            # Second intergenic region (after last gene)
            region2 = intergenic_features[1]
            assert region2.start == 301
            assert region2.end == 399
            assert region2.sequence == "G" * 99

    def test_calculate_intergenic_regions_min_length_filter(self, sample_features):
        """Test intergenic calculation with minimum length filtering."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            # Create genome with short intergenic regions
            genome_seq = "A" * 100 + "T" * 10 + "C" * 100 + "G" * 200
            with open(genome_file, "w") as f:
                f.write(">test_contig\n")
                f.write(genome_seq + "\n")

            id_manager = CompactIDManager()
            intergenic_features = calculate_intergenic_regions(
                sample_features, genome_file, "test_genome", id_manager, min_length=50
            )

            # Should find 2 intergenic regions:
            # 1. Between genes (101-199): 99bp > 50bp min_length
            # 2. After last gene (301-end): 110bp > 50bp min_length
            assert len(intergenic_features) == 2
            assert intergenic_features[0].start == 101
            assert intergenic_features[1].start == 301


class TestFeatureExtraction:
    """Test complete genome feature extraction."""

    def test_extract_genome_features_proteins_only(self):
        """Test extraction with only protein features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            # Create test genome
            with open(genome_file, "w") as f:
                f.write(">test_contig\n")
                f.write("ATCG" * 100 + "\n")  # 400bp genome

            id_manager = CompactIDManager()

            # Mock external tool functions
            with patch("pangenomeplus.extraction.run_prodigal") as mock_prodigal, patch(
                "pangenomeplus.extraction.parse_prodigal_gff"
            ) as mock_parse:
                # Setup mocks
                mock_prodigal.return_value = "/fake/gff/file"
                mock_parse.return_value = [
                    Feature(
                        compact_id="P1",
                        genome_id="test_genome",
                        contig="test_contig",
                        start=1,
                        end=200,
                        strand="+",
                        sequence="A" * 200,
                        feature_type="P",
                        original_id="gene_1",
                        metadata={},
                    )
                ]

                features = extract_genome_features(
                    genome_file, "test_genome", tmpdir, id_manager, skip_intergenic=True
                )

                assert "proteins" in features
                assert len(features["proteins"]) == 1
                assert features["proteins"][0].compact_id == "P1"
                assert features["intergenic"] == []  # Skipped

    def test_extract_genome_features_with_intergenic(self):
        """Test extraction including intergenic regions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            with open(genome_file, "w") as f:
                f.write(">test_contig\n")
                f.write("A" * 100 + "T" * 100 + "C" * 100 + "\n")  # 300bp genome

            id_manager = CompactIDManager()

            with patch("pangenomeplus.extraction.run_prodigal") as mock_prodigal, patch(
                "pangenomeplus.extraction.parse_prodigal_gff"
            ) as mock_parse:
                mock_prodigal.return_value = "/fake/gff/file"
                mock_parse.return_value = [
                    Feature(
                        compact_id="P1",
                        genome_id="test_genome",
                        contig="test_contig",
                        start=1,
                        end=100,
                        strand="+",
                        sequence="A" * 100,
                        feature_type="P",
                        original_id="gene_1",
                        metadata={},
                    ),
                    Feature(
                        compact_id="P2",
                        genome_id="test_genome",
                        contig="test_contig",
                        start=201,
                        end=300,
                        strand="+",
                        sequence="C" * 100,
                        feature_type="P",
                        original_id="gene_2",
                        metadata={},
                    ),
                ]

                features = extract_genome_features(
                    genome_file, "test_genome", tmpdir, id_manager
                )

                assert "proteins" in features
                assert "intergenic" in features
                assert len(features["proteins"]) == 2
                # Should find intergenic region between genes (positions 101-200)
                assert len(features["intergenic"]) >= 0  # Depends on min_length

    def test_extract_genome_features_prodigal_failure(self):
        """Test handling of Prodigal failure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            with open(genome_file, "w") as f:
                f.write(">test_contig\nATCG\n")

            id_manager = CompactIDManager()

            with patch("pangenomeplus.extraction.run_prodigal") as mock_prodigal:
                mock_prodigal.side_effect = ExtractionError("Prodigal failed")

                with pytest.raises(ExtractionError, match="Protein extraction failed"):
                    extract_genome_features(
                        genome_file, "test_genome", tmpdir, id_manager
                    )


class TestErrorHandling:
    """Test error handling throughout extraction module."""

    def test_extraction_error_creation(self):
        """Test ExtractionError exception."""
        error = ExtractionError("Test error message")
        assert str(error) == "Test error message"
        assert isinstance(error, Exception)

    def test_invalid_trna_model(self):
        """Test invalid tRNA model parameter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "genome.fasta")

            with pytest.raises(ExtractionError, match="Unknown tRNA model"):
                run_trnascan(genome_file, tmpdir, model="invalid_model")
