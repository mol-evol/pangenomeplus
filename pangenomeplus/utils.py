"""Utility functions for PanGenomePlus.

This module provides common utility functions used across the pipeline,
following DRY principles by centralizing frequently used operations.
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_genome(filepath: Union[str, Path]) -> SeqRecord:
    """Single source of truth for genome loading.

    Args:
        filepath: Path to genome FASTA file

    Returns:
        SeqRecord object containing the genome sequence

    Raises:
        FileNotFoundError: If genome file doesn't exist
        ValueError: If file is not valid FASTA format
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Genome file not found: {filepath}")

    try:
        # Try to read as single record first (most common case)
        record = SeqIO.read(filepath, "fasta")
        return record
    except ValueError:
        # If multiple records, read all and return as list
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            raise ValueError(f"No valid FASTA records found in {filepath}")
        if len(records) == 1:
            return records[0]
        # For multiple contigs, this should be handled differently
        # For now, log warning and return first record
        logging.warning(f"Multiple contigs found in {filepath}, using first record only")
        return records[0]


def validate_fasta(filepath: Union[str, Path]) -> bool:
    """Validate that a file is valid FASTA format.

    Args:
        filepath: Path to file to validate

    Returns:
        True if valid FASTA, False otherwise
    """
    filepath = Path(filepath)

    if not filepath.exists():
        return False

    try:
        records = list(SeqIO.parse(filepath, "fasta"))
        return len(records) > 0
    except (IOError, OSError, ValueError, UnicodeDecodeError) as e:
        logging.warning(f"Failed to validate FASTA file {filepath}: {e}")
        return False


def find_genome_files(genome_dir: Union[str, Path], extensions: Optional[List[str]] = None) -> List[Path]:
    """Find all genome files in a directory.

    Args:
        genome_dir: Directory to search
        extensions: List of valid extensions (default: [".fasta", ".fna", ".fa"])

    Returns:
        List of Path objects for found genome files
    """
    from .constants import GENOME_EXTENSIONS

    genome_dir = Path(genome_dir)
    if not genome_dir.exists():
        raise FileNotFoundError(f"Genome directory not found: {genome_dir}")

    if extensions is None:
        extensions = GENOME_EXTENSIONS

    genome_files = []
    for ext in extensions:
        pattern = f"*{ext}"
        genome_files.extend(genome_dir.glob(pattern))

    return sorted(genome_files)


def format_percentage(value: float, decimal_places: int = 1) -> str:
    """Format a float as a percentage string.

    Args:
        value: Float value between 0 and 1
        decimal_places: Number of decimal places

    Returns:
        Formatted percentage string
    """
    return f"{value * 100:.{decimal_places}f}%"


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """Safely divide two numbers, returning default if denominator is zero.

    Args:
        numerator: The numerator
        denominator: The denominator
        default: Value to return if denominator is zero

    Returns:
        Result of division or default value
    """
    if denominator == 0:
        return default
    return numerator / denominator