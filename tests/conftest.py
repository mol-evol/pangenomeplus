"""Pytest configuration and fixtures for PanGenomePlus tests.

All fixtures use real genomes and real tool outputs - no mocking.
"""

import os
import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def test_genomes_dir():
    """Path to directory containing real test genomes."""
    return Path(__file__).parent / "test_genomes"


@pytest.fixture
def single_ecoli_genome(test_genomes_dir):
    """Path to single E. coli genome for quick tests."""
    genome = test_genomes_dir / "single_ecoli" / "GCF_000005845.2.fasta"
    if not genome.exists():
        pytest.skip(f"Test genome not found: {genome}")
    return genome


@pytest.fixture
def three_ecoli_dir(test_genomes_dir):
    """Directory with three E. coli genomes for clustering tests."""
    genome_dir = test_genomes_dir / "three_genomes"
    if not genome_dir.exists():
        pytest.skip(f"Test genome directory not found: {genome_dir}")

    # Verify all three genomes are present
    expected = ["GCF_000005845.2.fasta", "GCF_000008865.2.fasta", "GCF_001280405.1.fasta"]
    for genome in expected:
        if not (genome_dir / genome).exists():
            pytest.skip(f"Missing test genome: {genome}")

    return genome_dir


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory for test results."""
    temp_dir = tempfile.mkdtemp(prefix="pangenome_test_")
    yield Path(temp_dir)
    # Cleanup after test
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def corrupted_genome_file(temp_output_dir):
    """Create a corrupted genome file for error testing."""
    corrupted = temp_output_dir / "corrupted.fasta"
    corrupted.write_text(">Incomplete_genome\nATCGATCG")  # No complete sequence
    return corrupted


@pytest.fixture
def empty_genome_file(temp_output_dir):
    """Create an empty file for error testing."""
    empty = temp_output_dir / "empty.fasta"
    empty.touch()
    return empty


@pytest.fixture
def invalid_fasta_file(temp_output_dir):
    """Create an invalid FASTA file for error testing."""
    invalid = temp_output_dir / "invalid.fasta"
    invalid.write_text("This is not a FASTA file\nJust random text")
    return invalid


@pytest.fixture
def pipeline_config():
    """Standard pipeline configuration for tests."""
    from pangenomeplus.pipeline import PipelineConfig

    return PipelineConfig(
        genome_dir="",  # Will be set per test
        output_dir="",  # Will be set per test
        skip_proteins=False,
        skip_intergenic=False,
        skip_trna=False,
        skip_rrna=False,
        skip_crispr=False,
        prodigal_mode="single",
        translation_table=11,
        trna_model="bacteria",
        rrna_kingdom="bac",
        threads=2,  # Limit threads for testing
        log_level="DEBUG",
        core_threshold=0.95,
        cloud_threshold=0.15,
        min_intergenic_length=50,
    )


@pytest.fixture(scope="session")
def check_tools_available():
    """Check that all required external tools are available."""
    required_tools = ["prodigal", "mmseqs", "barrnap", "tRNAscan-SE", "minced"]
    missing = []

    for tool in required_tools:
        if shutil.which(tool) is None:
            # Try alternative names
            if tool == "tRNAscan-SE" and shutil.which("trnascan-se"):
                continue
            missing.append(tool)

    if missing:
        pytest.skip(f"Required tools not found in PATH: {', '.join(missing)}")


# Biological expectation fixtures
@pytest.fixture
def ecoli_expectations():
    """Expected values for E. coli based on literature."""
    return {
        "total_genes": (4000, 4500),  # Range
        "proteins": (3900, 4400),
        "trnas": (80, 90),
        "rrnas": (21, 22),  # 7 operons × 3 types
        "core_percentage_3_genomes": (0.75, 0.85),  # For 3 E. coli strains
        "core_percentage_6_genomes": (0.70, 0.80),  # For 6 E. coli strains
    }


@pytest.fixture
def performance_expectations():
    """Performance expectations for the pipeline."""
    return {
        "extraction_per_genome_seconds": 60,  # Max 1 minute per genome
        "clustering_per_1000_features_seconds": 10,
        "memory_limit_mb": 1000,  # Max 1GB for typical genome
        "features_per_second": 100,  # Minimum processing rate
    }


# Markers for test categorization
def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "quick: Quick tests with single genome")
    config.addinivalue_line("markers", "slow: Slow integration tests")
    config.addinivalue_line("markers", "edge_cases: Edge case tests")
    config.addinivalue_line("markers", "requires_tools: Requires external tools")