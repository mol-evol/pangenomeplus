"""Unit tests for pipeline orchestration module."""

import json
import logging
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from pangenomeplus.compact_ids import CompactIDManager
from pangenomeplus.core import Feature
from pangenomeplus.extraction import ExtractionError
from pangenomeplus.pipeline import (
    PipelineConfig,
    PipelineError,
    ProcessingStats,
    create_checkpoint,
    create_pipeline_config,
    discover_genomes,
    load_checkpoint,
    process_single_genome,
    save_feature_mappings,
    setup_logging,
)


class TestPipelineConfig:
    """Test PipelineConfig dataclass and validation."""

    def test_config_creation_defaults(self):
        """Test config creation with default parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(genome_dir=tmpdir, output_dir=tmpdir)

            assert config.genome_dir == tmpdir
            assert config.output_dir == tmpdir
            assert config.skip_trna is False
            assert config.skip_rrna is False
            assert config.skip_crispr is False
            assert config.skip_intergenic is False
            assert config.protein_only is False
            assert config.non_coding_only is False
            assert config.resume is True
            assert config.keep_intermediates is True
            assert config.log_level == "INFO"

            # Check default parameters are initialized
            assert config.prodigal_params is not None
            assert config.trna_params is not None
            assert config.rrna_params is not None
            assert config.crispr_params is not None

    def test_config_custom_parameters(self):
        """Test config with custom parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            custom_prodigal = {"mode": "meta", "translation_table": 4}
            custom_trna = {"model": "archaea", "score_cutoff": 15.0}

            config = PipelineConfig(
                genome_dir=tmpdir,
                output_dir=tmpdir,
                skip_trna=True,
                protein_only=True,
                log_level="DEBUG",
                prodigal_params=custom_prodigal,
                trna_params=custom_trna,
            )

            assert config.skip_trna is True
            assert config.protein_only is True
            assert config.log_level == "DEBUG"
            assert config.prodigal_params == custom_prodigal
            assert config.trna_params == custom_trna

    def test_config_validation_success(self):
        """Test successful config validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(genome_dir=tmpdir, output_dir=tmpdir)
            config.validate()  # Should not raise

    def test_config_validation_bad_directory(self):
        """Test config validation with bad genome directory."""
        config = PipelineConfig(genome_dir="/nonexistent", output_dir="/tmp")

        with pytest.raises(ValueError, match="Genome directory does not exist"):
            config.validate()

    def test_config_validation_conflicting_flags(self):
        """Test config validation with conflicting feature flags."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(
                genome_dir=tmpdir,
                output_dir=tmpdir,
                protein_only=True,
                non_coding_only=True,
            )

            with pytest.raises(ValueError, match="Cannot specify both"):
                config.validate()

    def test_config_validation_no_features(self):
        """Test config validation when all features are skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(
                genome_dir=tmpdir,
                output_dir=tmpdir,
                skip_intergenic=True,
                skip_trna=True,
                skip_rrna=True,
                skip_crispr=True,
                protein_only=True,  # This combination leaves nothing
            )

            with pytest.raises(
                ValueError, match="Must analyze at least one feature type"
            ):
                config.validate()


class TestProcessingStats:
    """Test ProcessingStats dataclass."""

    def test_stats_initialization(self):
        """Test stats initialization."""
        stats = ProcessingStats()

        assert stats.total_genomes == 0
        assert stats.processed_genomes == 0
        assert stats.failed_genomes == 0
        assert stats.total_features == 0
        assert stats.feature_counts is not None
        assert len(stats.feature_counts) == 5

    def test_processing_time_calculation(self):
        """Test processing time calculations."""
        import time

        stats = ProcessingStats()
        stats.start_time = time.time()

        # Add small delay to ensure measurable time difference
        time.sleep(0.01)  # 10ms delay

        # Test ongoing processing time
        assert stats.processing_time > 0

        # Test completed processing time
        stats.end_time = stats.start_time + 10
        assert stats.processing_time == 10

    def test_genomes_per_second_calculation(self):
        """Test rate calculations."""
        import time

        stats = ProcessingStats()
        stats.start_time = time.time()
        stats.processed_genomes = 5
        stats.end_time = stats.start_time + 10

        assert stats.genomes_per_second == 0.5

        # Test with zero time
        stats.start_time = stats.end_time
        assert stats.genomes_per_second == 0.0

    def test_add_genome_features(self):
        """Test adding features to statistics."""
        stats = ProcessingStats()

        features = {
            "proteins": [MagicMock(), MagicMock(), MagicMock()],
            "intergenic": [MagicMock()],
            "tRNAs": [],
            "rRNAs": [MagicMock()],
            "CRISPR": [],
        }

        stats.add_genome_features(features)

        assert stats.feature_counts["proteins"] == 3
        assert stats.feature_counts["intergenic"] == 1
        assert stats.feature_counts["tRNAs"] == 0
        assert stats.feature_counts["rRNAs"] == 1
        assert stats.feature_counts["CRISPR"] == 0
        assert stats.total_features == 5


class TestUtilityFunctions:
    """Test utility functions."""

    def test_setup_logging_console_only(self):
        """Test logging setup with console handler only."""
        logger = setup_logging("WARNING")

        assert logger.level == logging.WARNING
        assert len(logger.handlers) == 1
        assert isinstance(logger.handlers[0], logging.StreamHandler)

    def test_setup_logging_with_file(self):
        """Test logging setup with file handler."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")

            logger = setup_logging("DEBUG", log_file)

            assert logger.level == logging.DEBUG
            assert len(logger.handlers) == 2  # Console + file
            assert os.path.exists(log_file)

            # Test actual logging
            logger.info("Test message")

            with open(log_file, "r") as f:
                content = f.read()
                assert "Test message" in content

    def test_discover_genomes_success(self):
        """Test successful genome discovery."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test genome files
            for i, ext in enumerate([".fasta", ".fna", ".fa"]):
                with open(os.path.join(tmpdir, f"genome_{i}{ext}"), "w") as f:
                    f.write(">test\nATCG\n")

            genomes = discover_genomes(tmpdir)

            assert len(genomes) == 3
            assert all(os.path.exists(g) for g in genomes)
            assert genomes == sorted(genomes)  # Should be sorted

    def test_discover_genomes_no_files(self):
        """Test genome discovery with no files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(PipelineError, match="No genome files found"):
                discover_genomes(tmpdir)

    def test_discover_genomes_mixed_files(self):
        """Test genome discovery ignoring non-genome files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create genome files
            with open(os.path.join(tmpdir, "genome1.fna"), "w") as f:
                f.write(">test\nATCG\n")

            # Create non-genome files
            with open(os.path.join(tmpdir, "readme.txt"), "w") as f:
                f.write("Not a genome")
            with open(os.path.join(tmpdir, "data.csv"), "w") as f:
                f.write("header,data\n")

            genomes = discover_genomes(tmpdir)

            assert len(genomes) == 1
            assert genomes[0].endswith("genome1.fna")

    def test_create_checkpoint(self):
        """Test checkpoint creation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint_file = os.path.join(tmpdir, "subdir", "checkpoint.json")
            data = {"test": "data", "number": 42}

            create_checkpoint(checkpoint_file, data)

            assert os.path.exists(checkpoint_file)

            with open(checkpoint_file, "r") as f:
                loaded = json.load(f)
                assert loaded == data

    def test_load_checkpoint_success(self):
        """Test successful checkpoint loading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint_file = os.path.join(tmpdir, "checkpoint.json")
            data = {"processed": ["genome1", "genome2"], "count": 2}

            # Create checkpoint
            with open(checkpoint_file, "w") as f:
                json.dump(data, f)

            loaded = load_checkpoint(checkpoint_file)

            assert loaded == data

    def test_load_checkpoint_missing_file(self):
        """Test loading non-existent checkpoint."""
        result = load_checkpoint("/nonexistent/checkpoint.json")
        assert result is None

    def test_load_checkpoint_corrupted_file(self):
        """Test loading corrupted checkpoint."""
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint_file = os.path.join(tmpdir, "checkpoint.json")

            # Create corrupted JSON
            with open(checkpoint_file, "w") as f:
                f.write("invalid json content")

            result = load_checkpoint(checkpoint_file)
            assert result is None

    def test_save_feature_mappings(self):
        """Test saving feature mappings."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Create mock ID manager with features
            id_manager = CompactIDManager()

            feature1 = Feature(
                compact_id="P1",
                genome_id="genome1",
                contig="contig1",
                start=1,
                end=100,
                strand="+",
                sequence="ATCG",
                feature_type="P",
                original_id="gene1",
                metadata={"product": "test protein"},
            )

            id_manager.register_feature(feature1)

            processed_genomes = ["genome1"]

            mappings_file = save_feature_mappings(
                output_dir, id_manager, processed_genomes
            )

            assert os.path.exists(mappings_file)

            with open(mappings_file, "r") as f:
                mappings = json.load(f)

                assert mappings["version"] == "1.0"
                assert mappings["processed_genomes"] == processed_genomes
                assert "P1" in mappings["compact_to_full"]
                assert mappings["compact_to_full"]["P1"]["genome_id"] == "genome1"

    def test_create_pipeline_config_convenience(self):
        """Test convenience function for creating pipeline config."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = create_pipeline_config(
                tmpdir, tmpdir, protein_only=True, log_level="DEBUG"
            )

            assert isinstance(config, PipelineConfig)
            assert config.genome_dir == tmpdir
            assert config.output_dir == tmpdir
            assert config.protein_only is True
            assert config.log_level == "DEBUG"


class TestSingleGenomeProcessing:
    """Test single genome processing function."""

    def test_process_single_genome_success(self):
        """Test successful single genome processing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create mock genome file
            genome_file = os.path.join(tmpdir, "test_genome.fna")
            with open(genome_file, "w") as f:
                f.write(">contig1\nATCGATCGATCG\n")

            # Create config
            config = PipelineConfig(genome_dir=tmpdir, output_dir=tmpdir)

            # Create ID manager
            id_manager = CompactIDManager()

            # Create logger
            logger = logging.getLogger("test")

            # Mock extraction function
            mock_features = {
                "proteins": [MagicMock()],
                "intergenic": [],
                "tRNAs": [],
                "rRNAs": [],
                "CRISPR": [],
            }

            with patch(
                "pangenomeplus.pipeline.extract_genome_features"
            ) as mock_extract:
                mock_extract.return_value = mock_features

                genome_id, features = process_single_genome(
                    genome_file, config, id_manager, logger
                )

                assert genome_id == "test_genome"
                assert features == mock_features

                # Verify extract_genome_features was called correctly
                mock_extract.assert_called_once()
                call_args = mock_extract.call_args

                # Check keyword arguments since the function uses keywords
                assert call_args.kwargs["genome_file"] == genome_file
                assert call_args.kwargs["genome_id"] == "test_genome"
                assert call_args.kwargs["id_manager"] == id_manager

    def test_process_single_genome_failure(self):
        """Test single genome processing failure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            genome_file = os.path.join(tmpdir, "test_genome.fna")
            with open(genome_file, "w") as f:
                f.write(">contig1\nATCG\n")

            config = PipelineConfig(genome_dir=tmpdir, output_dir=tmpdir)
            id_manager = CompactIDManager()
            logger = logging.getLogger("test")

            # Mock extraction to fail
            with patch(
                "pangenomeplus.pipeline.extract_genome_features"
            ) as mock_extract:
                mock_extract.side_effect = Exception("Extraction failed")

                with pytest.raises(
                    ExtractionError, match="Genome test_genome processing failed"
                ):
                    process_single_genome(genome_file, config, id_manager, logger)


class TestErrorHandling:
    """Test error handling throughout pipeline."""

    def test_pipeline_error_creation(self):
        """Test PipelineError exception."""
        error = PipelineError("Test error")
        assert str(error) == "Test error"
        assert isinstance(error, Exception)

    def test_extraction_error_propagation(self):
        """Test that ExtractionError is properly handled."""
        # This is tested in the process_single_genome_failure test
        # but we can add specific error propagation tests here if needed
        pass
