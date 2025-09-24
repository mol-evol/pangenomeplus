"""Integration tests for PanGenomePlus pipeline with real E. coli data."""

import json
import logging
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from pangenomeplus.compact_ids import CompactIDManager
from pangenomeplus.pipeline import (
    PipelineConfig,
    PipelineError,
    create_pipeline_config,
    discover_genomes,
    process_genomes,
    setup_logging,
)


class TestPipelineIntegration:
    """Integration tests for the complete pipeline."""

    @pytest.fixture
    def test_genome_dir(self):
        """Path to test E. coli genomes."""
        current_dir = Path(__file__).parent.parent
        genome_dir = current_dir / "test_data" / "genomes"

        if not genome_dir.exists():
            pytest.skip("Test genomes not available - run download script first")

        return str(genome_dir)

    @pytest.fixture
    def small_genome_set(self, test_genome_dir):
        """Create a small set of 2 genomes for quick testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy first 2 genomes to temp directory
            import shutil

            source_genomes = list(Path(test_genome_dir).glob("*.fna"))[:2]
            if len(source_genomes) < 2:
                pytest.skip("Need at least 2 test genomes")

            for genome in source_genomes:
                shutil.copy(genome, tmpdir)

            yield tmpdir

    def test_discover_genomes_success(self, test_genome_dir):
        """Test genome file discovery."""
        genomes = discover_genomes(test_genome_dir)

        assert len(genomes) > 0
        assert all(g.endswith(".fna") for g in genomes)
        assert all(os.path.exists(g) for g in genomes)

        # Should be sorted
        assert genomes == sorted(genomes)

    def test_discover_genomes_no_files(self):
        """Test genome discovery with no files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(PipelineError, match="No genome files found"):
                discover_genomes(tmpdir)

    def test_pipeline_config_validation(self):
        """Test pipeline configuration validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Valid config
            config = create_pipeline_config(tmpdir, tmpdir)
            config.validate()  # Should not raise

            # Invalid genome directory
            config_bad = create_pipeline_config("/nonexistent", tmpdir)
            with pytest.raises(ValueError, match="Genome directory does not exist"):
                config_bad.validate()

            # Conflicting flags
            config_conflict = create_pipeline_config(
                tmpdir, tmpdir, protein_only=True, non_coding_only=True
            )
            with pytest.raises(ValueError, match="Cannot specify both"):
                config_conflict.validate()

    def test_logging_setup(self):
        """Test logging configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")

            logger = setup_logging("INFO", log_file)

            assert logger.level == logging.INFO
            assert os.path.exists(log_file)

            # Test logging
            logger.info("Test message")

            with open(log_file, "r") as f:
                content = f.read()
                assert "Test message" in content

    @pytest.mark.slow
    def test_end_to_end_pipeline_small_set(self, small_genome_set):
        """Test complete pipeline with 2 genomes."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Configure for protein-only analysis (faster)
            config = create_pipeline_config(
                genome_dir=small_genome_set,
                output_dir=output_dir,
                protein_only=True,
                log_level="DEBUG",
            )

            # Run pipeline
            stats = process_genomes(config)

            # Validate results
            assert stats.total_genomes == 2
            assert stats.processed_genomes == 2
            assert stats.failed_genomes == 0
            assert stats.total_features > 0
            assert stats.feature_counts["proteins"] > 0

            # Check output files exist
            mappings_file = os.path.join(
                output_dir, "mappings", "compact_id_mappings.json"
            )
            assert os.path.exists(mappings_file)

            # Validate mappings file
            with open(mappings_file, "r") as f:
                mappings = json.load(f)
                assert "version" in mappings
                assert "compact_to_full" in mappings
                assert len(mappings["compact_to_full"]) > 0
                assert "processed_genomes" in mappings
                assert len(mappings["processed_genomes"]) == 2

    def test_checkpoint_system(self, small_genome_set):
        """Test checkpoint and resume functionality."""
        with tempfile.TemporaryDirectory() as output_dir:
            config = create_pipeline_config(
                genome_dir=small_genome_set,
                output_dir=output_dir,
                protein_only=True,
                resume=True,
            )

            # Mock extraction to fail on second genome
            call_count = 0

            def mock_extract_fail_second(*args, **kwargs):
                nonlocal call_count
                call_count += 1
                if call_count == 2:  # Fail on second genome
                    raise Exception("Simulated failure")
                # Return mock successful extraction for first genome
                from unittest.mock import MagicMock

                return {
                    "proteins": [MagicMock()],
                    "intergenic": [],
                    "tRNAs": [],
                    "rRNAs": [],
                    "CRISPR": [],
                }

            # First run - should fail partway through
            with patch(
                "pangenomeplus.pipeline.extract_genome_features"
            ) as mock_extract:
                mock_extract.side_effect = mock_extract_fail_second

                stats1 = process_genomes(config)

                # Should have processed 1 genome, failed 1
                assert stats1.processed_genomes == 1
                assert stats1.failed_genomes == 1

            # Reset mock for second run
            call_count = 0

            # Second run - should resume and complete
            stats2 = process_genomes(config)

            # Should complete all genomes
            assert stats2.processed_genomes == 2
            assert stats2.failed_genomes == 0

    def test_pipeline_with_feature_selection(self, small_genome_set):
        """Test pipeline with different feature selection options."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Test protein-only
            config_proteins = create_pipeline_config(
                genome_dir=small_genome_set,
                output_dir=os.path.join(output_dir, "proteins"),
                protein_only=True,
            )

            stats = process_genomes(config_proteins)

            assert stats.feature_counts["proteins"] > 0
            assert stats.feature_counts["intergenic"] == 0
            assert stats.feature_counts["tRNAs"] == 0

    def test_error_handling_bad_genome(self):
        """Test error handling with corrupted genome file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create bad genome file
            bad_genome = os.path.join(tmpdir, "bad_genome.fna")
            with open(bad_genome, "w") as f:
                f.write("This is not a valid FASTA file")

            with tempfile.TemporaryDirectory() as output_dir:
                config = create_pipeline_config(tmpdir, output_dir)

                # Should handle error gracefully
                stats = process_genomes(config)

                assert stats.total_genomes == 1
                assert stats.processed_genomes == 0
                assert stats.failed_genomes == 1

    def test_performance_timing(self, small_genome_set):
        """Test performance timing and statistics."""
        with tempfile.TemporaryDirectory() as output_dir:
            config = create_pipeline_config(
                genome_dir=small_genome_set,
                output_dir=output_dir,
                protein_only=True,  # Faster processing
            )

            stats = process_genomes(config)

            # Verify timing calculations
            assert stats.processing_time > 0
            assert stats.genomes_per_second >= 0
            assert stats.start_time > 0
            assert stats.end_time > stats.start_time

    @pytest.mark.slow
    def test_full_feature_extraction(self, test_genome_dir):
        """Test extraction of all feature types with one genome."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy just one genome for this test
            import shutil

            source_genomes = list(Path(test_genome_dir).glob("*.fna"))
            if not source_genomes:
                pytest.skip("No test genomes available")

            # Use K-12 MG1655 as it's well-characterized
            k12_genome = None
            for genome in source_genomes:
                if "K12" in genome.name or "MG1655" in genome.name:
                    k12_genome = genome
                    break

            if not k12_genome:
                k12_genome = source_genomes[0]  # Use first available

            shutil.copy(k12_genome, tmpdir)

            with tempfile.TemporaryDirectory() as output_dir:
                config = create_pipeline_config(
                    genome_dir=tmpdir,
                    output_dir=output_dir,
                    # Enable all features (this will be slow with real tools)
                    skip_trna=False,
                    skip_rrna=False,
                    skip_crispr=False,
                    skip_intergenic=False,
                )

                # Mock external tools for speed (integration test focuses on pipeline orchestration)
                with patch(
                    "pangenomeplus.extraction.run_prodigal"
                ) as mock_prodigal, patch(
                    "pangenomeplus.extraction.parse_prodigal_gff"
                ) as mock_parse:
                    # Mock successful tool execution
                    mock_prodigal.return_value = "/mock/output.gff"
                    mock_parse.return_value = []  # Empty results for speed

                    stats = process_genomes(config)

                    # Verify pipeline completed
                    assert stats.total_genomes == 1
                    assert stats.processed_genomes == 1
                    assert stats.failed_genomes == 0


class TestRealToolIntegration:
    """Tests that actually run external tools (marked as slow)."""

    @pytest.fixture
    def single_genome_dir(self):
        """Create directory with single small E. coli genome for tool testing."""
        current_dir = Path(__file__).parent.parent
        genome_dir = current_dir / "test_data" / "genomes"

        if not genome_dir.exists():
            pytest.skip("Test genomes not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy smallest genome
            import shutil

            genomes = list(genome_dir.glob("*.fna"))
            if not genomes:
                pytest.skip("No test genomes available")

            # Find smallest genome for fastest testing
            smallest = min(genomes, key=lambda x: x.stat().st_size)
            shutil.copy(smallest, tmpdir)

            yield tmpdir

    @pytest.mark.slow
    @pytest.mark.real_tools
    def test_real_prodigal_execution(self, single_genome_dir):
        """Test pipeline with actual Prodigal execution."""
        with tempfile.TemporaryDirectory() as output_dir:
            config = create_pipeline_config(
                genome_dir=single_genome_dir,
                output_dir=output_dir,
                protein_only=True,  # Only run Prodigal
                log_level="DEBUG",
            )

            stats = process_genomes(config)

            # Should successfully extract proteins
            assert stats.processed_genomes == 1
            assert stats.feature_counts["proteins"] > 1000  # E. coli has ~4000 genes

            # Verify output structure
            mappings_file = os.path.join(
                output_dir, "mappings", "compact_id_mappings.json"
            )
            assert os.path.exists(mappings_file)

            # Check intermediate files were created
            intermediate_dir = os.path.join(output_dir, "intermediate")
            assert os.path.exists(intermediate_dir)

    @pytest.mark.slow
    @pytest.mark.real_tools
    def test_tool_failure_handling(self, single_genome_dir):
        """Test handling of external tool failures."""
        with tempfile.TemporaryDirectory() as output_dir:
            config = create_pipeline_config(
                genome_dir=single_genome_dir, output_dir=output_dir, protein_only=True
            )

            # Mock Prodigal to fail
            with patch("pangenomeplus.extraction.subprocess.run") as mock_run:
                mock_run.side_effect = Exception("Tool failure")

                stats = process_genomes(config)

                # Should handle failure gracefully
                assert stats.processed_genomes == 0
                assert stats.failed_genomes == 1
