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
import random
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .clustering import ClusteringError, FamilyStats, cluster_features_by_type
from .compact_ids import CompactIDManager
from .core import Feature
from .extraction import ExtractionError, extract_genome_features


def extract_genome_name(genome_file: str) -> str:
    """Extract a display-friendly genome name from file path."""
    return Path(genome_file).stem


def format_time(seconds: float) -> str:
    """Format seconds into MM:SS string."""
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    return f"{minutes:02d}:{secs:02d}"


def generate_playful_message() -> str:
    """Generate a dynamic playful status message for genomic analysis.

    Returns:
        A randomly generated phrase in format: "verb modifier noun"
        Example: "wrangling prokaryotic widgets"
    """
    # Actions the software performs
    verbs = [
        "wrangling",
        "herding",
        "clustering",
        "sorting",
        "fishing",
        "cataloging",
        "hunting",
        "mining",
        "juggling",
        "assembling",
        "aligning",
        "mapping",
        "partitioning",
        "annotating",
        "detecting",
        "untangling",
        "collecting",
        "organizing",
        "indexing",
        "comparing",
        "disassembling",
        "decombobulating",
        "discombobulating",
        "bamboozling",
        "flummoxing",
        "befuddling",
        "hornswoggling",
        "finagling",
        "jiggery-pokering",
        "kerfuffling",
        "rambunctifying",
        "shenaniganning",
        "hullaballooing",
        "perambulating",
        "confabulating",
        "recalcitrating",
        "obfuscating",
        "triangulating",
        "extrapolating",
        "prognosticating",
        "pontificating",
        "defenestrating",
        "discombobifying",
        "flibbertigibbeting",
        "skedaddling",
        "gallivanting",
        "traipsing",
        "meandering",
        "floundering",
        "fumbling",
        "tinkering",
        "fiddling",
        "doodling",
        "noodling",
    ]

    # Prokaryotic/genomic descriptors
    modifiers = [
        "bug",
        "microbe",
        "prokaryotic",
        "homologous",
        "core-genome",
        "accessory",
        "plasmid",
        "chromosomal",
        "prophage",
        "CRISPR",
        "operon",
        "genetic",
        "metagenomic",
        "syntenic",
        "pangenomic",
        "xenologous",
        "mobile",
        "conjugative",
        "transposon",
        "genomic",
        "buggy",
        "germy",
        "critter",
        "microbe-y",
        "beastie",
        "essential",
        "bonus",
        "extra",
        "spare",
        "optional",
        "freeloading",
        "hitchhiking",
        "wandering",
        "roaming",
        "nomadic",
        "zombie",
        "lurking",
        "sneaky",
        "hidden",
        "dormant",
        "cousin",
        "sibling",
        "family",
        "kinfolk",
        "related",
        "neighborly",
        "adjacent",
        "side-by-side",
        "clustered",
        "huddled",
        "choppy",
        "snippy",
        "scissor-happy",
        "editor",
        "pruning",
        "copy-paste",
        "duplicated",
        "redundant",
        "repetitive",
        "cloned",
        "jumping",
        "hopping",
        "vagrant",
        "drifting",
        "core",
        "backbone",
        "skeleton",
        "fundamental",
        "bread-and-butter",
    ]

    # Whimsical objects
    nouns = [
        "confetti",
        "sprinkles",
        "jellybeans",
        "breadcrumbs",
        "marbles",
        "beads",
        "nuggets",
        "clusters",
        "widgets",
        "doodads",
        "gizmos",
        "baubles",
        "trinkets",
        "bits",
        "chunks",
        "morsels",
        "tidbits",
        "specks",
        "crumbs",
        "fragments",
        "modules",
        "cassettes",
        "islands",
        "patches",
        "segments",
        "cartridges",
        "packets",
        "bundles",
        "arrays",
        "scaffolds",
        "contigs",
        "motifs",
        "domains",
        "elements",
        "loci",
        "regions",
        "zones",
        "stretches",
        "blocks",
        "strings",
        "strands",
        "units",
        "subunits",
        "components",
        "assemblies",
        "constructs",
        "repertoires",
        "collections",
        "libraries",
        "catalogs",
        "inventories",
    ]

    return f"{random.choice(verbs)} {random.choice(modifiers)} {random.choice(nouns)}"


class SingleLineProgress:
    """Single-line progress indicator with spinner and metrics."""

    def __init__(self, total: int, playful: bool = True) -> None:
        """Initialize progress tracker.

        Args:
            total: Total number of items to process
            playful: Whether to show playful status messages
        """
        self.total = total
        self.current = 0
        self.playful = playful
        self.spinners = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        self.spinner_idx = 0
        self.start_time = time.time()
        self.terminal_width = 120  # Conservative terminal width

    def update(self, genome_name: str, status_msg: str = "") -> None:
        """Update progress display.

        Args:
            genome_name: Name of current genome being processed
            status_msg: Playful status message (optional)
        """
        # Calculate metrics
        elapsed = time.time() - self.start_time
        rate = self.current / elapsed if elapsed > 0 else 0
        remaining_items = self.total - self.current
        remaining_time = remaining_items / rate if rate > 0 else 0
        pct = int(100 * self.current / self.total) if self.total > 0 else 0

        # Get spinner character
        spinner = self.spinners[self.spinner_idx % len(self.spinners)]
        self.spinner_idx += 1

        # Build message components
        if self.playful and status_msg:
            status_part = f"{status_msg} | "
        else:
            status_part = "Processing | "

        # Shorten genome name if needed
        max_genome_len = 30
        if len(genome_name) > max_genome_len:
            genome_name = genome_name[: max_genome_len - 3] + "..."

        progress_part = f"[{self.current}/{self.total}] {pct}%"
        rate_part = f"{rate:.1f} g/s"
        time_part = format_time(remaining_time)

        # Build complete message
        msg = (
            f"{spinner} {status_part}{genome_name} {progress_part} • "
            f"{rate_part} • {time_part} remaining"
        )

        # Clear line and write message
        sys.stdout.write("\r" + " " * self.terminal_width + "\r")
        sys.stdout.write(msg[: self.terminal_width])
        sys.stdout.flush()

    def advance(self) -> None:
        """Increment progress counter."""
        self.current += 1

    def finish(self) -> None:
        """Complete progress and move to next line."""
        sys.stdout.write("\n")
        sys.stdout.flush()


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
    prodigal_params: Dict[str, Any] = field(default_factory=dict)
    trna_params: Dict[str, Any] = field(default_factory=dict)
    rrna_params: Dict[str, Any] = field(default_factory=dict)
    crispr_params: Dict[str, Any] = field(default_factory=dict)

    # Performance settings
    resume: bool = True
    keep_intermediates: bool = True
    log_level: str = "INFO"

    # Annotation options
    use_existing_annotations: bool = False

    # Downstream analysis options
    enable_downstream_analysis: bool = False
    rarefaction_iterations: int = 100
    rarefaction_step_size: int = 1

    # UI options
    playful_mode: bool = True

    def __post_init__(self) -> None:
        """Initialize default parameters if not provided."""
        if not self.prodigal_params:
            self.prodigal_params = {"mode": "single", "translation_table": 11}
        if not self.trna_params:
            self.trna_params = {"model": "bacteria", "score_cutoff": 20.0}
        if not self.rrna_params:
            self.rrna_params = {"kingdom": "bac", "evalue": 1e-6}
        if not self.crispr_params:
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
    feature_counts: Dict[str, int] = field(default_factory=dict)
    start_time: float = 0
    end_time: float = 0

    # Clustering statistics
    total_families: int = 0
    family_counts: Dict[str, int] = field(default_factory=dict)
    singleton_counts: Dict[str, int] = field(default_factory=dict)
    core_families: int = 0
    accessory_families: int = 0
    cloud_families: int = 0

    def __post_init__(self) -> None:
        """Initialize feature counts."""
        if not self.feature_counts:
            self.feature_counts = {
                "proteins": 0,
                "intergenic": 0,
                "tRNAs": 0,
                "rRNAs": 0,
                "CRISPR": 0,
            }
        if not self.family_counts:
            self.family_counts = {
                "P": 0,  # Protein families
                "I": 0,  # Intergenic families
                "T": 0,  # tRNA families
                "R": 0,  # rRNA families
                "C": 0,  # CRISPR families
            }
        if not self.singleton_counts:
            self.singleton_counts = {
                "P": 0,  # Protein singletons
                "I": 0,  # Intergenic singletons
                "T": 0,  # tRNA singletons
                "R": 0,  # rRNA singletons
                "C": 0,  # CRISPR singletons
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

    def add_family_stats(self, family_stats: Dict[str, FamilyStats]) -> None:
        """Add family statistics from clustering results."""
        for family_id, stats in family_stats.items():
            feature_type = stats.feature_type

            if feature_type in self.family_counts:
                if stats.classification == "singleton":
                    self.singleton_counts[feature_type] += 1
                else:
                    self.family_counts[feature_type] += 1
                    self.total_families += 1

                    # Update classification counts
                    if stats.classification == "core":
                        self.core_families += 1
                    elif stats.classification == "accessory":
                        self.accessory_families += 1
                    elif stats.classification == "cloud":
                        self.cloud_families += 1


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
            return json.load(f)  # type: ignore[no-any-return]
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
        mapping_data["compact_to_full"][compact_id] = {  # type: ignore[index]
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
    for (
        genome_id,
        contig,
        start,
        end,
    ), compact_id in id_manager.location_to_compact.items():
        key = f"{genome_id}: {contig}: {start}: {end}"
        mapping_data["location_to_compact"][key] = compact_id  # type: ignore[index]

    with open(mappings_file, "w") as f:
        json.dump(mapping_data, f, indent=2)

    return mappings_file


def run_clustering_stage(
    id_manager: CompactIDManager,
    config: PipelineConfig,
    total_genomes: int,
    logger: logging.Logger,
) -> Dict[str, Dict[str, FamilyStats]]:
    """Run clustering stage on all extracted features.

    Args:
        id_manager: CompactIDManager with all features
        config: Pipeline configuration
        total_genomes: Total number of genomes for classification
        logger: Logger instance

    Returns:
        Dictionary mapping feature type to family statistics

    Raises:
        ClusteringError: If clustering fails
    """
    logger.info("Starting clustering stage")

    clustering_output_dir = os.path.join(config.output_dir, "clustering")
    os.makedirs(clustering_output_dir, exist_ok=True)

    # Feature type mapping for clustering
    feature_type_map = {
        "proteins": "P",
        "intergenic": "I",
        "tRNAs": "T",
        "rRNAs": "R",
        "CRISPR": "C",
    }

    all_family_stats = {}

    for feature_name, feature_type in feature_type_map.items():
        # Skip if feature type disabled by configuration
        if _should_skip_feature_type(feature_type, config):
            logger.debug(f"Skipping {feature_name} clustering (disabled by config)")
            continue

        # Get all features of this type
        features = []
        for compact_id, feature in id_manager.compact_to_full.items():
            if feature.feature_type == feature_type:
                features.append(feature)

        if not features:
            logger.info(f"No {feature_name} features found for clustering")
            continue

        logger.info(f"Clustering {len(features)} {feature_name} features")

        try:
            # Run clustering for this feature type
            compact_to_family, family_stats = cluster_features_by_type(
                features=features,
                feature_type=feature_type,
                output_dir=clustering_output_dir,
                id_manager=id_manager,
                total_genomes=total_genomes,
            )

            all_family_stats[feature_type] = family_stats

            # Update id_manager with family assignments (if we extend it later)
            multi_member = [
                f for f in family_stats.values() if f.classification != "singleton"
            ]
            logger.info(
                f"  {feature_name}: {len(family_stats)} families "
                f"({len(multi_member)} multi-member families)"
            )

        except ClusteringError as e:
            logger.error(f"Clustering failed for {feature_name}: {e}")
            # Continue with other feature types
            continue

    logger.info("Clustering stage completed")
    return all_family_stats


def _should_skip_feature_type(feature_type: str, config: PipelineConfig) -> bool:
    """Check if feature type should be skipped based on configuration."""
    skip_map = {
        "P": config.non_coding_only,  # Skip proteins if non-coding only
        "I": config.skip_intergenic or config.protein_only,
        "T": config.skip_trna or config.protein_only,
        "R": config.skip_rrna or config.protein_only,
        "C": config.skip_crispr or config.protein_only,
    }
    return skip_map.get(feature_type, False)


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
            use_existing_annotations=config.use_existing_annotations,
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
    5. Run MMseqs2 clustering on all features
    6. Assign family IDs and classify as core/accessory/cloud
    7. Generate final mappings and statistics

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

    # Display header if playful mode is enabled
    if config.playful_mode:
        print(f"\nPanGenomePlus | Analyzing {stats.total_genomes} genomes")
        print(f"Output directory: {config.output_dir}\n")

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

    if config.playful_mode:
        # Use single-line progress indicator for playful mode
        progress = SingleLineProgress(
            total=len(genome_files), playful=config.playful_mode
        )

        for i, genome_file in enumerate(genome_files):
            genome_id = Path(genome_file).stem

            # Skip if already processed (resume functionality)
            if genome_id in processed_genomes:
                logger.debug(f"Skipping already processed genome: {genome_id}")
                progress.advance()
                continue

            genome_name = extract_genome_name(genome_file)
            status = generate_playful_message()

            # Update progress display
            progress.update(genome_name, status)

            try:
                # Process single genome
                processed_id, features = process_single_genome(
                    genome_file, config, id_manager, logger
                )

                # Update statistics
                stats.add_genome_features(features)
                stats.processed_genomes += 1
                processed_genomes.append(processed_id)

                # Log completion with feature counts to file
                protein_count = len(features.get("proteins", []))
                intergenic_count = len(features.get("intergenic", []))
                trna_count = len(features.get("tRNAs", []))
                rrna_count = len(features.get("rRNAs", []))
                crispr_count = len(features.get("CRISPR", []))

                logger.info(
                    f"✓ {genome_name}: {protein_count:,} proteins | "
                    f"{intergenic_count:,} intergenic | {trna_count} tRNAs | "
                    f"{rrna_count} rRNAs | {crispr_count} CRISPR"
                )

            except ExtractionError as e:
                logger.warning(f"Genome processing failed, continuing: {e}")
                logger.warning(f"✗ {genome_name}: Processing failed")
                failed_genomes.append(genome_id)
                stats.failed_genomes += 1

            progress.advance()

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

        # Finish progress and move to new line
        progress.finish()

    else:
        # Traditional mode for production environments
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

    # Run clustering stage if any genomes were processed
    if stats.processed_genomes > 0 and stats.total_features > 0:
        logger.info("Starting clustering stage")
        try:
            all_family_stats = run_clustering_stage(
                id_manager=id_manager,
                config=config,
                total_genomes=stats.processed_genomes,
                logger=logger,
            )

            # Update statistics with clustering results
            for feature_type, family_stats in all_family_stats.items():
                stats.add_family_stats(family_stats)

            logger.info("Clustering stage completed successfully")

        except ClusteringError as e:
            logger.error(f"Clustering stage failed: {e}")
            logger.warning("Continuing without clustering results")

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

    # Display completion message if playful mode is enabled
    if config.playful_mode:
        print("\n✓ Analysis complete")

    logger.info("Feature extraction summary:")
    for feature_type, count in stats.feature_counts.items():
        if count > 0:
            logger.info(f"  {feature_type}: {count:,} features")

    # Clustering statistics
    if stats.total_families > 0 or any(stats.family_counts.values()):
        logger.info("Clustering summary:")
        logger.info(f"  Total families: {stats.total_families:,}")
        logger.info(f"  Core families: {stats.core_families:,}")
        logger.info(f"  Accessory families: {stats.accessory_families:,}")
        logger.info(f"  Cloud families: {stats.cloud_families:,}")

        logger.info("Families by feature type:")
        feature_names = {
            "P": "Proteins",
            "I": "Intergenic",
            "T": "tRNAs",
            "R": "rRNAs",
            "C": "CRISPR",
        }
        for feature_type, count in stats.family_counts.items():
            if count > 0:
                singletons = stats.singleton_counts.get(feature_type, 0)
                logger.info(
                    f"  {feature_names.get(feature_type, feature_type)}: "
                    f"{count:,} families, {singletons:,} singletons"
                )

    if failed_genomes:
        logger.warning(f"Failed genomes: {', '.join(failed_genomes)}")

    # Clean up checkpoint if all successful
    if stats.failed_genomes == 0 and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
        logger.debug("Checkpoint file removed after successful completion")

    return stats


def create_pipeline_config(
    genome_dir: str, output_dir: str, **kwargs: Any
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
