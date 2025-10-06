"""Tests for utility functions using real genome data."""

import pytest
from pathlib import Path

from pangenomeplus.utils import (
    load_genome,
    validate_fasta,
    find_genome_files,
    format_percentage,
    safe_divide,
)


class TestLoadGenome:
    """Test genome loading with real E. coli genome."""

    def test_load_real_ecoli_genome(self, single_ecoli_genome):
        """Test loading a real E. coli genome."""
        record = load_genome(single_ecoli_genome)

        # Verify it's a real E. coli genome
        assert record is not None
        assert len(record.seq) > 4_000_000  # E. coli is ~4.6 Mb
        assert len(record.seq) < 6_000_000
        assert "coli" in record.description.lower()

    def test_load_genome_file_not_found(self):
        """Test FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError, match="Genome file not found"):
            load_genome("nonexistent.fasta")

    def test_load_empty_file(self, empty_genome_file):
        """Test loading an empty file."""
        with pytest.raises(ValueError, match="No valid FASTA records"):
            load_genome(empty_genome_file)

    def test_load_invalid_fasta(self, invalid_fasta_file):
        """Test loading an invalid FASTA file."""
        with pytest.raises(ValueError, match="No valid FASTA records"):
            load_genome(invalid_fasta_file)

    def test_load_genome_path_types(self, single_ecoli_genome):
        """Test loading with both string and Path objects."""
        # Test with Path object
        record1 = load_genome(single_ecoli_genome)
        # Test with string
        record2 = load_genome(str(single_ecoli_genome))

        assert len(record1.seq) == len(record2.seq)


class TestValidateFasta:
    """Test FASTA validation with real and invalid files."""

    def test_validate_real_genome(self, single_ecoli_genome):
        """Test validation of real E. coli genome."""
        assert validate_fasta(single_ecoli_genome) is True

    def test_validate_nonexistent_file(self):
        """Test validation of non-existent file."""
        assert validate_fasta("nonexistent.fasta") is False

    def test_validate_empty_file(self, empty_genome_file):
        """Test validation of empty file."""
        assert validate_fasta(empty_genome_file) is False

    def test_validate_invalid_fasta(self, invalid_fasta_file):
        """Test validation of invalid FASTA."""
        assert validate_fasta(invalid_fasta_file) is False

    def test_validate_corrupted_genome(self, corrupted_genome_file):
        """Test validation of corrupted genome."""
        # Corrupted but still valid FASTA format
        assert validate_fasta(corrupted_genome_file) is True


class TestFindGenomeFiles:
    """Test genome file discovery with real test genomes."""

    def test_find_genomes_in_directory(self, three_ecoli_dir):
        """Test finding genomes in directory with 3 E. coli."""
        genomes = find_genome_files(three_ecoli_dir)

        assert len(genomes) == 3
        assert all(g.suffix == ".fasta" for g in genomes)
        assert all(g.exists() for g in genomes)

    def test_find_genomes_with_extensions(self, test_genomes_dir):
        """Test finding genomes with specific extensions."""
        # Only .fasta files
        fasta_files = find_genome_files(test_genomes_dir / "three_genomes", [".fasta"])
        assert len(fasta_files) == 3

        # Non-existent extension
        no_files = find_genome_files(test_genomes_dir / "three_genomes", [".gbk"])
        assert len(no_files) == 0

    def test_find_genomes_nonexistent_dir(self):
        """Test error when directory doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Genome directory not found"):
            find_genome_files("nonexistent_directory")

    def test_find_genomes_empty_dir(self, temp_output_dir):
        """Test finding genomes in empty directory."""
        genomes = find_genome_files(temp_output_dir)
        assert len(genomes) == 0


class TestFormattingUtilities:
    """Test formatting and math utility functions."""

    def test_format_percentage(self):
        """Test percentage formatting."""
        assert format_percentage(0.95) == "95.0%"
        assert format_percentage(0.95, 2) == "95.00%"
        assert format_percentage(0.12345, 1) == "12.3%"
        assert format_percentage(1.0) == "100.0%"
        assert format_percentage(0.0) == "0.0%"

    def test_safe_divide(self):
        """Test safe division."""
        assert safe_divide(10, 2) == 5.0
        assert safe_divide(10, 0) == 0.0  # Default
        assert safe_divide(10, 0, -1) == -1  # Custom default
        assert safe_divide(0, 10) == 0.0
        assert safe_divide(1, 3) == pytest.approx(0.333333, rel=1e-5)


class TestRealGenomeProperties:
    """Test utilities with properties of real genomes."""

    def test_ecoli_genome_properties(self, single_ecoli_genome):
        """Validate E. coli genome has expected properties."""
        record = load_genome(single_ecoli_genome)

        # Check GC content (E. coli ~50.8%)
        gc_count = record.seq.count('G') + record.seq.count('C')
        total = len(record.seq)
        gc_content = safe_divide(gc_count, total)

        assert 0.48 < gc_content < 0.52, f"Unexpected GC content: {gc_content}"

        # Check it's a complete genome (single contig)
        assert record.id.startswith("NC_")  # RefSeq complete genome

    def test_multiple_genomes_in_dir(self, three_ecoli_dir):
        """Test loading multiple real genomes from directory."""
        genomes = find_genome_files(three_ecoli_dir)
        records = []

        for genome_file in genomes:
            record = load_genome(genome_file)
            records.append(record)

        # All should be E. coli genomes of similar size
        sizes = [len(r.seq) for r in records]
        assert all(4_000_000 < s < 6_000_000 for s in sizes)

        # Check variation between strains (should be <10%)
        avg_size = sum(sizes) / len(sizes)
        variations = [abs(s - avg_size) / avg_size for s in sizes]
        assert all(v < 0.1 for v in variations)