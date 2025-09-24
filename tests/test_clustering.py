"""Tests for clustering functionality."""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from pangenomeplus.clustering import (
    ClusteringError,
    ClusteringParams,
    FamilyStats,
    _get_reverse_complement,
    assign_family_ids,
    classify_families,
    cluster_features_by_type,
    create_feature_fasta,
    get_clustering_params,
    parse_cluster_results,
    run_mmseqs2_clustering,
)
from pangenomeplus.compact_ids import CompactIDManager
from pangenomeplus.core import Feature


class TestClusteringParams:
    """Test clustering parameter configuration."""

    def test_get_clustering_params_proteins(self):
        """Test parameter settings for proteins."""
        params = get_clustering_params("P")
        assert params.min_seq_id == 0.8
        assert params.coverage == 0.8
        assert params.coverage_mode == 1

    def test_get_clustering_params_intergenic(self):
        """Test relaxed parameters for intergenic regions."""
        params = get_clustering_params("I")
        assert params.min_seq_id == 0.7
        assert params.coverage == 0.5

    def test_get_clustering_params_trna(self):
        """Test strict parameters for tRNAs."""
        params = get_clustering_params("T")
        assert params.min_seq_id == 0.9
        assert params.coverage == 0.9

    def test_get_clustering_params_rrna(self):
        """Test very strict parameters for rRNAs."""
        params = get_clustering_params("R")
        assert params.min_seq_id == 0.95
        assert params.coverage == 0.95

    def test_get_clustering_params_crispr(self):
        """Test parameters for CRISPR spacers."""
        params = get_clustering_params("C")
        assert params.min_seq_id == 0.8
        assert params.coverage == 0.6

    def test_get_clustering_params_invalid(self):
        """Test error for invalid feature type."""
        with pytest.raises(ClusteringError, match="Invalid feature type: X"):
            get_clustering_params("X")


class TestReverseComplement:
    """Test reverse complement functionality."""

    def test_simple_sequence(self):
        """Test reverse complement of simple DNA sequence."""
        seq = "ATCG"
        result = _get_reverse_complement(seq)
        assert result == "CGAT"

    def test_complex_sequence(self):
        """Test reverse complement with ambiguous bases."""
        seq = "ATCGRYWSMKBDHVN"
        result = _get_reverse_complement(seq)
        assert result == "NBDHVMKSWRYCGAT"

    def test_lowercase(self):
        """Test reverse complement preserves case."""
        seq = "atcg"
        result = _get_reverse_complement(seq)
        assert result == "cgat"

    def test_mixed_case(self):
        """Test reverse complement with mixed case."""
        seq = "AtCg"
        result = _get_reverse_complement(seq)
        assert result == "cGaT"

    def test_empty_sequence(self):
        """Test empty sequence."""
        result = _get_reverse_complement("")
        assert result == ""


class TestFastaGeneration:
    """Test FASTA file generation."""

    def test_create_feature_fasta_unidirectional(self):
        """Test basic FASTA generation."""
        features = [
            Feature(
                compact_id="P1",
                genome_id="test",
                contig="contig1",
                start=1,
                end=100,
                strand="+",
                sequence="ATGAAACCC",
                feature_type="P",
                original_id="gene1",
                metadata={},
            ),
            Feature(
                compact_id="P2",
                genome_id="test",
                contig="contig1",
                start=200,
                end=300,
                strand="+",
                sequence="ATGCCCAAA",
                feature_type="P",
                original_id="gene2",
                metadata={},
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "test.fasta")
            result = create_feature_fasta(features, output_file)

            assert result == output_file
            assert os.path.exists(output_file)

            with open(output_file, "r") as f:
                content = f.read()

            expected_content = ">P1\nATGAAACCC\n>P2\nATGCCCAA\n"
            assert ">P1\nATGAAACCC" in content
            assert ">P2\nATGCCCAA" in content

    def test_create_feature_fasta_bidirectional(self):
        """Test FASTA generation with reverse complements."""
        features = [
            Feature(
                compact_id="I1",
                genome_id="test",
                contig="contig1",
                start=1,
                end=10,
                strand=".",
                sequence="ATCG",
                feature_type="I",
                original_id="inter1",
                metadata={},
            )
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "test.fasta")
            create_feature_fasta(features, output_file, bidirectional=True)

            with open(output_file, "r") as f:
                content = f.read()

            assert ">I1\nATCG\n" in content
            assert ">I1_RC\nCGAT\n" in content

    def test_create_feature_fasta_no_features(self):
        """Test error when no features provided."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "test.fasta")

            with pytest.raises(ClusteringError, match="No features provided"):
                create_feature_fasta([], output_file)

    def test_create_feature_fasta_empty_sequence(self):
        """Test handling of features with empty sequences."""
        features = [
            Feature(
                compact_id="P1",
                genome_id="test",
                contig="contig1",
                start=1,
                end=100,
                strand="+",
                sequence="ATGAAA",
                feature_type="P",
                original_id="gene1",
                metadata={},
            ),
            Feature(
                compact_id="P2",
                genome_id="test",
                contig="contig1",
                start=200,
                end=300,
                strand="+",
                sequence="",  # Empty sequence
                feature_type="P",
                original_id="gene2",
                metadata={},
            ),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "test.fasta")
            create_feature_fasta(features, output_file)

            with open(output_file, "r") as f:
                content = f.read()

            # Should only contain P1, not P2 (empty sequence)
            assert ">P1\nATGAAA" in content
            assert "P2" not in content


class TestMMseqs2Integration:
    """Test MMseqs2 integration (mocked)."""

    @patch("pangenomeplus.clustering.subprocess.run")
    def test_run_mmseqs2_clustering_success(self, mock_run):
        """Test successful MMseqs2 execution."""
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = ""

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input FASTA file
            fasta_file = os.path.join(tmpdir, "input.fasta")
            with open(fasta_file, "w") as f:
                f.write(">P1\nATGAAA\n>P2\nATGCCC\n")

            # Create mock TSV output
            tsv_file = os.path.join(tmpdir, "cluster_results.tsv")
            with open(tsv_file, "w") as f:
                f.write("P1\tP1\n")
                f.write("P1\tP2\n")

            params = ClusteringParams()
            result = run_mmseqs2_clustering(fasta_file, tmpdir, params)

            assert result == tsv_file
            assert mock_run.call_count == 3  # createdb, cluster, createtsv

    @patch("pangenomeplus.clustering.subprocess.run")
    def test_run_mmseqs2_clustering_failure(self, mock_run):
        """Test MMseqs2 execution failure."""
        mock_run.return_value.returncode = 1
        mock_run.return_value.stderr = "MMseqs2 error"

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_file = os.path.join(tmpdir, "input.fasta")
            with open(fasta_file, "w") as f:
                f.write(">P1\nATGAAA\n")

            params = ClusteringParams()

            with pytest.raises(ClusteringError, match="mmseqs createdb failed"):
                run_mmseqs2_clustering(fasta_file, tmpdir, params)

    def test_run_mmseqs2_clustering_missing_input(self):
        """Test error when input file missing."""
        params = ClusteringParams()

        with pytest.raises(ClusteringError, match="FASTA file not found"):
            run_mmseqs2_clustering("/nonexistent/file.fasta", "/tmp", params)


class TestClusterResultParsing:
    """Test cluster result parsing."""

    def test_parse_cluster_results_simple(self):
        """Test parsing simple cluster results."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = os.path.join(tmpdir, "clusters.tsv")
            with open(tsv_file, "w") as f:
                f.write("P1\tP1\n")
                f.write("P1\tP2\n")
                f.write("P3\tP3\n")

            clusters = parse_cluster_results(tsv_file)

            expected = {"P1": ["P1", "P2"], "P3": ["P3"]}
            assert clusters == expected

    def test_parse_cluster_results_bidirectional(self):
        """Test parsing bidirectional cluster results."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = os.path.join(tmpdir, "clusters.tsv")
            with open(tsv_file, "w") as f:
                f.write("I1\tI1\n")
                f.write("I1\tI2_RC\n")
                f.write("I3_RC\tI3_RC\n")

            clusters = parse_cluster_results(tsv_file, bidirectional=True)

            expected = {"I1": ["I1", "I2"], "I3": ["I3"]}
            assert clusters == expected

    def test_parse_cluster_results_missing_file(self):
        """Test error when cluster file missing."""
        with pytest.raises(ClusteringError, match="Cluster results file not found"):
            parse_cluster_results("/nonexistent/file.tsv")

    def test_parse_cluster_results_empty_lines(self):
        """Test handling of empty lines."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = os.path.join(tmpdir, "clusters.tsv")
            with open(tsv_file, "w") as f:
                f.write("P1\tP1\n")
                f.write("\n")  # Empty line
                f.write("P1\tP2\n")

            clusters = parse_cluster_results(tsv_file)
            assert clusters == {"P1": ["P1", "P2"]}


class TestFamilyAssignment:
    """Test family ID assignment."""

    def test_assign_family_ids_proteins(self):
        """Test family assignment for proteins."""
        clusters = {
            "P1": ["P1", "P2", "P3"],  # Multi-member family
            "P4": ["P4"],  # Singleton
        }

        compact_to_family, family_to_members = assign_family_ids(clusters, "P")

        # Check family assignments
        assert compact_to_family["P1"] == "FAM_P1"
        assert compact_to_family["P2"] == "FAM_P1"
        assert compact_to_family["P3"] == "FAM_P1"
        assert compact_to_family["P4"] == "SING_P4"

        # Check reverse mapping
        assert family_to_members["FAM_P1"] == ["P1", "P2", "P3"]
        assert family_to_members["SING_P4"] == ["P4"]

    def test_assign_family_ids_multiple_families(self):
        """Test assignment of multiple families."""
        clusters = {"P1": ["P1", "P2"], "P3": ["P3", "P4", "P5"], "P6": ["P6"]}

        compact_to_family, family_to_members = assign_family_ids(clusters, "P")

        assert compact_to_family["P1"] == "FAM_P1"
        assert compact_to_family["P3"] == "FAM_P2"
        assert compact_to_family["P6"] == "SING_P6"

    def test_assign_family_ids_invalid_type(self):
        """Test error for invalid feature type."""
        clusters = {"P1": ["P1"]}

        with pytest.raises(ClusteringError, match="Invalid feature type: X"):
            assign_family_ids(clusters, "X")


class TestFamilyClassification:
    """Test family classification."""

    def test_classify_families_core_accessory_cloud(self):
        """Test classification into core, accessory, and cloud families."""
        # Set up mock ID manager
        id_manager = CompactIDManager()

        # Create features in different genomes
        features = [
            Feature("P1", "genome1", "contig1", 1, 100, "+", "ATG", "P", "gene1", {}),
            Feature("P2", "genome2", "contig1", 1, 100, "+", "ATG", "P", "gene2", {}),
            Feature("P3", "genome3", "contig1", 1, 100, "+", "ATG", "P", "gene3", {}),
            Feature("P4", "genome4", "contig1", 1, 100, "+", "ATG", "P", "gene4", {}),
            Feature("P5", "genome1", "contig1", 200, 300, "+", "ATG", "P", "gene5", {}),
            Feature("P6", "genome1", "contig1", 400, 500, "+", "ATG", "P", "gene6", {}),
        ]

        for feature in features:
            id_manager.register_feature(feature)

        family_to_members = {
            "FAM_P1": ["P1", "P2", "P3", "P4"],  # Core (4/4 genomes)
            "FAM_P2": ["P5", "P6"],  # Cloud (1/4 genomes)
            "SING_P7": ["P7"],  # Singleton
        }

        family_stats = classify_families(
            family_to_members,
            id_manager,
            total_genomes=4,
            core_threshold=0.95,
            cloud_threshold=0.30,  # Set higher threshold so 25% is cloud
        )

        assert family_stats["FAM_P1"].classification == "core"
        assert family_stats["FAM_P1"].genome_count == 4

        assert family_stats["FAM_P2"].classification == "cloud"
        assert family_stats["FAM_P2"].genome_count == 1

        assert family_stats["SING_P7"].classification == "singleton"

    def test_classify_families_accessory(self):
        """Test accessory family classification."""
        id_manager = CompactIDManager()

        # Create features for accessory classification
        features = [
            Feature("P1", "genome1", "contig1", 1, 100, "+", "ATG", "P", "gene1", {}),
            Feature("P2", "genome2", "contig1", 1, 100, "+", "ATG", "P", "gene2", {}),
        ]

        for feature in features:
            id_manager.register_feature(feature)

        family_to_members = {
            "FAM_P1": ["P1", "P2"]  # Present in 2/4 genomes = 50% (accessory)
        }

        family_stats = classify_families(family_to_members, id_manager, total_genomes=4)

        assert family_stats["FAM_P1"].classification == "accessory"
        assert family_stats["FAM_P1"].genome_count == 2


class TestIntegrationClustering:
    """Integration tests for complete clustering workflow."""

    def test_cluster_features_by_type_mock(self):
        """Test complete clustering workflow with mocked MMseqs2."""
        # Create test features
        features = [
            Feature(
                "P1", "genome1", "contig1", 1, 100, "+", "ATGAAA", "P", "gene1", {}
            ),
            Feature(
                "P2", "genome2", "contig1", 1, 100, "+", "ATGCCC", "P", "gene2", {}
            ),
            Feature(
                "P3", "genome1", "contig1", 200, 300, "+", "ATGTTT", "P", "gene3", {}
            ),
        ]

        id_manager = CompactIDManager()
        for feature in features:
            id_manager.register_feature(feature)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Mock MMseqs2 execution
            with patch(
                "pangenomeplus.clustering.run_mmseqs2_clustering"
            ) as mock_mmseqs:
                mock_tsv = os.path.join(tmpdir, "mock_results.tsv")
                with open(mock_tsv, "w") as f:
                    f.write("P1\tP1\n")
                    f.write("P1\tP2\n")
                    f.write("P3\tP3\n")

                mock_mmseqs.return_value = mock_tsv

                compact_to_family, family_stats = cluster_features_by_type(
                    features=features,
                    feature_type="P",
                    output_dir=tmpdir,
                    id_manager=id_manager,
                    total_genomes=2,
                )

                # Check family assignments
                assert compact_to_family["P1"] == "FAM_P1"
                assert compact_to_family["P2"] == "FAM_P1"
                assert compact_to_family["P3"] == "SING_P3"

                # Check family statistics
                assert "FAM_P1" in family_stats
                assert "SING_P3" in family_stats

                fam_p1_stats = family_stats["FAM_P1"]
                assert fam_p1_stats.member_count == 2
                assert fam_p1_stats.genome_count == 2
                assert fam_p1_stats.classification == "core"

    def test_cluster_features_by_type_empty(self):
        """Test clustering with no features."""
        id_manager = CompactIDManager()

        with tempfile.TemporaryDirectory() as tmpdir:
            compact_to_family, family_stats = cluster_features_by_type(
                features=[],
                feature_type="P",
                output_dir=tmpdir,
                id_manager=id_manager,
                total_genomes=1,
            )

            assert compact_to_family == {}
            assert family_stats == {}
