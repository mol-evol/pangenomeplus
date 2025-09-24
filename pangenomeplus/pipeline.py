"""Pipeline orchestration and workflow control for PanGenomePlus.

This module provides the main orchestration functions for processing genomes
through the complete PanGenomePlus pipeline, following KISS principles with
function-first architecture and comprehensive checkpoint/recovery system.

Pipeline Flow:
Raw Genomes → Annotation Tools → Feature Extraction → Compact ID Assignment → Output

Key Features:
- Checkpoint system for recovery from failures
- Progress monitoring and structured logging
- Graceful error handling with partial results
- Memory-efficient genome-by-genome processing
"""

import json
import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .compact_ids import CompactIDManager
from .core import Feature
from .extraction import ExtractionError, extract_genome_features


@dataclass
class PipelineConfig:
    """Configuration for pipeline execution."""

    # Input/Output paths
    genome_dir: str
    output_dir: str

    # Feature selection
    skip_trna: bool = False
    skip_rrna: bool = False
    skip_crispr: bool = False
    skip_intergenic: bool = False
    protein_only: bool = False
    non_coding_only: bool = False

    # Tool parameters
    prodigal_params: Dict[str, Any] = None
    trna_params: Dict[str, Any] = None
    rrna_params: Dict[str, Any] = None
    crispr_params: Dict[str, Any] = None

    # Performance settings
    resume: bool = True
    keep_intermediates: bool = True
    log_level: str = "INFO"

    def __post_init__(self):
        """Initialize default parameters if not provided."""
        if self.prodigal_params is None:
            self.prodigal_params = {"mode": "single", "translation_table": 11}
        if self.trna_params is None:
            self.trna_params = {"model": "bacteria", "score_cutoff": 20.0}
        if self.rrna_params is None:
            self.rrna_params = {"kingdom": "bac", "evalue": 1e-6}
        if self.crispr_params is None:
            self.crispr_params = {
                "min_repeats": 3,
                "min_spacer_length": 26,
                "max_spacer_length": 50,
            }

    def validate(self) -> None:
        """Validate configuration parameters."""
        if not os.path.isdir(self.genome_dir):
            raise ValueError(f"Genome directory does not exist: {self.genome_dir}")

        if self.protein_only and self.non_coding_only:
            raise ValueError("Cannot specify both --protein-only and --non-coding-only")

        if self.skip_intergenic and all(
            [self.skip_trna, self.skip_rrna, self.skip_crispr, self.protein_only]
        ):
            raise ValueError("Must analyze at least one feature type")


@dataclass
class ProcessingStats:
    """Statistics for pipeline processing."""

    total_genomes: int = 0
    processed_genomes: int = 0
    failed_genomes: int = 0
    total_features: int = 0
    feature_counts: Dict[str, int] = None
    start_time: float = 0
    end_time: float = 0

    def __post_init__(self):
        """Initialize feature counts."""
        if self.feature_counts is None:
            self.feature_counts = {
                "proteins": 0,
                "intergenic": 0,
                "tRNAs": 0,
                "rRNAs": 0,
                "CRISPR": 0,
            }

    @property
    def processing_time(self) -> float:
        """Total processing time in seconds."""
        if self.end_time > 0:
            return self.end_time - self.start_time
        return time.time() - self.start_time

    @property
    def genomes_per_second(self) -> float:
        """Processing rate in genomes per second."""
        if self.processing_time > 0:
            return self.processed_genomes / self.processing_time
        return 0.0

    def add_genome_features(self, features: Dict[str, List[Feature]]) -> None:
        """Add features from a processed genome to statistics."""
        for feature_type, feature_list in features.items():
            if feature_type in self.feature_counts:
                count = len(feature_list)
                self.feature_counts[feature_type] += count
                self.total_features += count


class PipelineError(Exception):
    """Exception raised during pipeline execution."""

    pass


def setup_logging(
    log_level: str = "INFO", log_file: Optional[str] = None
) -> logging.Logger:
    """Set up structured logging for pipeline execution.

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional log file path

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger("pangenomeplus")
    logger.setLevel(getattr(logging, log_level.upper()))

    # Remove existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def discover_genomes(genome_dir: str) -> List[str]:
    """Discover genome FASTA files in the specified directory.

    Args:
        genome_dir: Directory containing genome files

    Returns:
        List of genome file paths

    Raises:
        PipelineError: If no genome files found
    """
    extensions = [".fasta", ".fna", ".fa"]
    genome_files = []

    for ext in extensions:
        pattern = f"*{ext}"
        files = list(Path(genome_dir).glob(pattern))
        genome_files.extend([str(f) for f in files])

    if not genome_files:
        raise PipelineError(f"No genome files found in {genome_dir}")

    # Sort for consistent processing order
    return sorted(genome_files)


def create_checkpoint(checkpoint_file: str, data: Dict[str, Any]) -> None:
    """Create a checkpoint file with processing state.

    Args:
        checkpoint_file: Path to checkpoint file
        data: Data to save in checkpoint
    """
    os.makedirs(os.path.dirname(checkpoint_file), exist_ok=True)
    with open(checkpoint_file, "w") as f:
        json.dump(data, f, indent=2)


def load_checkpoint(checkpoint_file: str) -> Optional[Dict[str, Any]]:
    """Load checkpoint data if it exists.

    Args:
        checkpoint_file: Path to checkpoint file

    Returns:
        Checkpoint data or None if file doesn't exist
    """
    if not os.path.exists(checkpoint_file):
        return None

    try:
        with open(checkpoint_file, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return None


def save_feature_mappings(
    output_dir: str, id_manager: CompactIDManager, processed_genomes: List[str]
) -> str:
    """Save compact ID mappings to output directory.

    Args:
        output_dir: Output directory
        id_manager: CompactIDManager with mappings
        processed_genomes: List of processed genome IDs

    Returns:
        Path to saved mappings file
    """
    mappings_dir = os.path.join(output_dir, "mappings")
    os.makedirs(mappings_dir, exist_ok=True)

    mappings_file = os.path.join(mappings_dir, "compact_id_mappings.json")

    # Create serializable mapping data
    mapping_data = {
        "version": "1.0",
        "timestamp": time.time(),
        "processed_genomes": processed_genomes,
        "compact_to_full": {},
        "location_to_compact": {},
        "genome_features": dict(id_manager.genome_features),
        "counters": dict(id_manager._counters),
    }

    # Convert Feature objects to dicts for serialization
    for compact_id, feature in id_manager.compact_to_full.items():
        mapping_data["compact_to_full"][compact_id] = {
            "compact_id": feature.compact_id,
            "genome_id": feature.genome_id,
            "contig": feature.contig,
            "start": feature.start,
            "end": feature.end,
            "strand": feature.strand,
            "sequence": feature.sequence,
            "feature_type": feature.feature_type,
            "original_id": feature.original_id,
            "metadata": feature.metadata,
        }

    # Convert location tuples to strings for JSON serialization
    for (genome_id, start, end), compact_id in id_manager.location_to_compact.items():
        key = f"{genome_id}:{start}:{end}"
        mapping_data["location_to_compact"][key] = compact_id

    with open(mappings_file, "w") as f:
        json.dump(mapping_data, f, indent=2)

    return mappings_file


def process_single_genome(
    genome_file: str,
    config: PipelineConfig,
    id_manager: CompactIDManager,
    logger: logging.Logger,
) -> Tuple[str, Dict[str, List[Feature]]]:
    """Process a single genome through the extraction pipeline.

    Args:
        genome_file: Path to genome FASTA file
        config: Pipeline configuration
        id_manager: CompactIDManager for ID assignment
        logger: Logger instance

    Returns:
        Tuple of (genome_id, extracted_features)

    Raises:
        ExtractionError: If genome processing fails
    """
    genome_id = Path(genome_file).stem
    logger.info(f"Processing genome: {genome_id}")

    # Create genome-specific output directory
    genome_output_dir = os.path.join(config.output_dir, "intermediate", genome_id)

    try:
        # Extract features using existing extraction functions
        # Apply protein_only logic by setting skip flags appropriately
        skip_intergenic = config.skip_intergenic or config.protein_only
        skip_trna = config.skip_trna or config.protein_only
        skip_rrna = config.skip_rrna or config.protein_only
        skip_crispr = config.skip_crispr or config.protein_only

        features = extract_genome_features(
            genome_file=genome_file,
            genome_id=genome_id,
            output_dir=genome_output_dir,
            id_manager=id_manager,
            skip_trna=skip_trna,
            skip_rrna=skip_rrna,
            skip_crispr=skip_crispr,
            skip_intergenic=skip_intergenic,
            prodigal=config.prodigal_params,
            trna=config.trna_params,
            rrna=config.rrna_params,
            crispr=config.crispr_params,
        )

        # Log feature counts
        for feature_type, feature_list in features.items():
            count = len(feature_list)
            if count > 0:
                logger.info(f"  {feature_type}: {count} features")

        return genome_id, features

    except Exception as e:
        logger.error(f"Failed to process {genome_id}: {e}")
        raise ExtractionError(f"Genome {genome_id} processing failed: {e}")


def process_genomes(config: PipelineConfig) -> ProcessingStats:
    """Main orchestration function for processing multiple genomes.

    Implements the complete PanGenomePlus pipeline with checkpoint/recovery:
    1. Discover genome files
    2. Process each genome through extraction pipeline
    3. Assign compact IDs universally
    4. Save intermediate results with checkpoints
    5. Generate final mappings and statistics

    Args:
        config: Pipeline configuration

    Returns:
        Processing statistics

    Raises:
        PipelineError: If pipeline setup or critical errors occur
    """
    # Validate configuration
    config.validate()

    # Set up output directories
    os.makedirs(config.output_dir, exist_ok=True)

    # Set up logging
    log_file = os.path.join(config.output_dir, "logs", "pipeline.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = setup_logging(config.log_level, log_file)

    logger.info("Starting PanGenomePlus pipeline")
    logger.info(f"Output directory: {config.output_dir}")

    # Initialize statistics
    stats = ProcessingStats()
    stats.start_time = time.time()

    # Discover genome files
    try:
        genome_files = discover_genomes(config.genome_dir)
        stats.total_genomes = len(genome_files)
        logger.info(f"Discovered {stats.total_genomes} genome files")
    except PipelineError as e:
        logger.error(f"Genome discovery failed: {e}")
        raise

    # Initialize compact ID manager
    id_manager = CompactIDManager()

    # Checkpoint file for recovery
    checkpoint_file = os.path.join(config.output_dir, "checkpoint.json")

    # Load previous progress if resuming
    processed_genomes = []
    if config.resume:
        checkpoint_data = load_checkpoint(checkpoint_file)
        if checkpoint_data:
            processed_genomes = checkpoint_data.get("processed_genomes", [])
            stats.processed_genomes = checkpoint_data.get("processed_count", 0)
            stats.failed_genomes = checkpoint_data.get("failed_count", 0)
            logger.info(
                f"Resuming from checkpoint: {stats.processed_genomes} genomes completed"
            )

    # Process each genome
    failed_genomes = []

    for i, genome_file in enumerate(genome_files):
        genome_id = Path(genome_file).stem

        # Skip if already processed (resume functionality)
        if genome_id in processed_genomes:
            logger.debug(f"Skipping already processed genome: {genome_id}")
            continue

        try:
            # Process single genome
            processed_id, features = process_single_genome(
                genome_file, config, id_manager, logger
            )

            # Update statistics
            stats.add_genome_features(features)
            stats.processed_genomes += 1
            processed_genomes.append(processed_id)

            # Progress reporting
            percent_complete = ((i + 1) / stats.total_genomes) * 100
            rate = stats.genomes_per_second
            logger.info(
                f"Progress: {percent_complete:.1f}% complete "
                f"({stats.processed_genomes}/{stats.total_genomes}) "
                f"Rate: {rate:.2f} genomes/sec"
            )

        except ExtractionError as e:
            logger.warning(f"Genome processing failed, continuing: {e}")
            failed_genomes.append(genome_id)
            stats.failed_genomes += 1
            continue

        # Create checkpoint every 10 genomes
        if (i + 1) % 10 == 0 or (i + 1) == len(genome_files):
            checkpoint_data = {
                "processed_genomes": processed_genomes,
                "processed_count": stats.processed_genomes,
                "failed_count": stats.failed_genomes,
                "failed_genomes": failed_genomes,
                "feature_counts": stats.feature_counts,
                "timestamp": time.time(),
            }
            create_checkpoint(checkpoint_file, checkpoint_data)
            logger.debug(f"Checkpoint created at genome {i + 1}")

    # Finalize processing
    stats.end_time = time.time()

    # Save final mappings
    if stats.processed_genomes > 0:
        mappings_file = save_feature_mappings(
            config.output_dir, id_manager, processed_genomes
        )
        logger.info(f"Compact ID mappings saved to: {mappings_file}")

    # Final statistics
    logger.info("Pipeline completed successfully")
    logger.info(f"Processed genomes: {stats.processed_genomes}")
    logger.info(f"Failed genomes: {stats.failed_genomes}")
    logger.info(f"Total features extracted: {stats.total_features}")
    logger.info(f"Processing time: {stats.processing_time:.1f} seconds")
    logger.info(f"Rate: {stats.genomes_per_second:.2f} genomes/second")

    for feature_type, count in stats.feature_counts.items():
        if count > 0:
            logger.info(f"  {feature_type}: {count:,} features")

    if failed_genomes:
        logger.warning(f"Failed genomes: {', '.join(failed_genomes)}")

    # Clean up checkpoint if all successful
    if stats.failed_genomes == 0 and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
        logger.debug("Checkpoint file removed after successful completion")

    return stats


def create_pipeline_config(
    genome_dir: str, output_dir: str, **kwargs
) -> PipelineConfig:
    """Convenience function to create pipeline configuration.

    Args:
        genome_dir: Directory containing genome FASTA files
        output_dir: Output directory for results
        **kwargs: Additional configuration parameters

    Returns:
        Configured PipelineConfig instance
    """
    return PipelineConfig(genome_dir=genome_dir, output_dir=output_dir, **kwargs)
