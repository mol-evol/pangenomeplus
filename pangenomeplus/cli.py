#!/usr/bin/env python3
"""Command-line interface for PanGenomePlus.

This module provides the complete command-line interface for PanGenomePlus
with 24 parameters for comprehensive pangenome analysis configuration.

Usage:
    pangenomeplus --genome-dir input_genomes/ --output-dir results/
    pangenomeplus --help
"""

import argparse
import logging
import os
import sys
from pathlib import Path

from .constants import FeatureType
from .pipeline import PipelineConfig, process_genomes


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser with all 24 parameters."""
    parser = argparse.ArgumentParser(
        prog="pangenomeplus",
        description="Comprehensive pangenome analysis pipeline for "
        "all genomic features",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic analysis with 8 genomes
    pangenomeplus --genome-dir genomes/ --output-dir results/

    # Protein-only analysis
    pangenomeplus --genome-dir genomes/ --output-dir results/ --protein-only

    # Skip specific feature types
    pangenomeplus --genome-dir genomes/ --output-dir results/ --skip-trna --skip-crispr

    # Custom clustering parameters
    pangenomeplus --genome-dir genomes/ --output-dir results/ \\
        --clustering-identity 0.9 --clustering-coverage 0.9

For more information, visit: https://github.com/your-org/pangenomeplus
        """,
    )

    # Required arguments
    required = parser.add_argument_group("Required Arguments")
    required.add_argument(
        "--genome-dir",
        type=str,
        required=True,
        help="Directory containing genome FASTA files",
    )
    required.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for analysis results",
    )

    # Prodigal parameters (2 parameters)
    prodigal = parser.add_argument_group("Prodigal Gene Prediction Parameters")
    prodigal.add_argument(
        "--prodigal-mode",
        choices=["single", "meta"],
        default="single",
        help="Prodigal prediction mode (default: single)",
    )
    prodigal.add_argument(
        "--translation-table",
        type=int,
        choices=list(range(1, 26)),
        default=11,
        help="Genetic code translation table (default: 11 for bacteria)",
    )

    # MMseqs2 clustering parameters (6 parameters)
    clustering = parser.add_argument_group("MMseqs2 Clustering Parameters")
    clustering.add_argument(
        "--clustering-identity",
        type=float,
        default=0.8,
        help="Sequence identity threshold for clustering (default: 0.8)",
    )
    clustering.add_argument(
        "--clustering-coverage",
        type=float,
        default=0.8,
        help="Coverage threshold for clustering (default: 0.8)",
    )
    clustering.add_argument(
        "--clustering-sensitivity",
        type=float,
        default=4.0,
        help="MMseqs2 sensitivity level (default: 4.0)",
    )
    clustering.add_argument(
        "--coverage-mode",
        choices=["query", "target", "bidirectional"],
        default="bidirectional",
        help="Coverage calculation method (default: bidirectional)",
    )
    clustering.add_argument(
        "--cluster-mode",
        choices=["set-cover", "connected-component", "greedy"],
        default="set-cover",
        help="MMseqs2 clustering algorithm (default: set-cover)",
    )
    clustering.add_argument(
        "--max-seqs",
        type=int,
        default=10000,
        help="Maximum sequences per query for memory management (default: 10000)",
    )

    # tRNA detection parameters (3 parameters)
    trna = parser.add_argument_group("tRNAscan-SE Parameters")
    trna.add_argument(
        "--trna-model",
        choices=["bacteria", "archaea", "eukaryota", "mitochondrial"],
        default="bacteria",
        help="Organism-specific tRNA model (default: bacteria)",
    )
    trna.add_argument(
        "--trna-score-cutoff",
        type=float,
        default=20.0,
        help="Score threshold for tRNA detection (default: 20.0)",
    )
    trna.add_argument(
        "--trna-search-mode",
        choices=["relaxed", "normal", "strict"],
        default="normal",
        help="tRNA search stringency (default: normal)",
    )

    # rRNA detection parameters (2 parameters)
    rrna = parser.add_argument_group("Barrnap rRNA Detection Parameters")
    rrna.add_argument(
        "--rrna-kingdom",
        choices=["bac", "arc", "euk", "mito"],
        default="bac",
        help="Kingdom-specific rRNA model (default: bac)",
    )
    rrna.add_argument(
        "--rrna-evalue",
        type=float,
        default=1e-6,
        help="E-value threshold for rRNA detection (default: 1e-6)",
    )

    # CRISPR detection parameters (3 parameters)
    crispr = parser.add_argument_group("MINCED CRISPR Detection Parameters")
    crispr.add_argument(
        "--crispr-min-repeats",
        type=int,
        default=3,
        help="Minimum repeat count for CRISPR array (default: 3)",
    )
    crispr.add_argument(
        "--crispr-repeat-length-range",
        type=str,
        default="23,47",
        help="Repeat length range as MIN,MAX (default: 23,47)",
    )
    crispr.add_argument(
        "--crispr-spacer-length-range",
        type=str,
        default="26,50",
        help="Spacer length range as MIN,MAX (default: 26,50)",
    )

    # Performance control parameters (2 parameters)
    performance = parser.add_argument_group("Performance Control Parameters")
    performance.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Number of parallel threads (default: auto-detect)",
    )
    performance.add_argument(
        "--memory-limit",
        type=str,
        default="auto",
        help="Memory limit for MMseqs2 operations (default: auto)",
    )

    # Annotation input parameters
    annotation = parser.add_argument_group("Annotation Input Parameters")
    annotation.add_argument(
        "--use-existing-annotations",
        action="store_true",
        help="Use existing GFF3 files if found (default: False)",
    )

    # Feature selection parameters (6 parameters)
    features = parser.add_argument_group("Feature Selection Parameters")
    features.add_argument(
        "--skip-trna",
        action="store_true",
        help="Skip tRNA analysis (excludes tRNAscan-SE)",
    )
    features.add_argument(
        "--skip-rrna", action="store_true", help="Skip rRNA analysis (excludes Barrnap)"
    )
    features.add_argument(
        "--skip-crispr",
        action="store_true",
        help="Skip CRISPR analysis (excludes MINCED)",
    )
    features.add_argument(
        "--skip-intergenic",
        action="store_true",
        help="Skip intergenic region extraction",
    )
    features.add_argument(
        "--protein-only",
        action="store_true",
        help="Analyze only protein-coding genes (excludes all non-coding features)",
    )
    features.add_argument(
        "--non-coding-only",
        action="store_true",
        help="Analyze only non-coding features (excludes proteins)",
    )

    # Downstream analysis parameters (3 parameters)
    downstream = parser.add_argument_group("Downstream Analysis Parameters")
    downstream.add_argument(
        "--enable-downstream-analysis",
        action="store_true",
        help=(
            "Enable downstream pangenome analyses (rarefaction curves, "
            "openness analysis, summary reports)"
        ),
    )
    downstream.add_argument(
        "--rarefaction-iterations",
        type=int,
        default=100,
        help="Number of iterations for rarefaction curve calculation (default: 100)",
    )
    downstream.add_argument(
        "--rarefaction-step-size",
        type=int,
        default=1,
        help="Step size for genome count in rarefaction analysis (default: 1)",
    )

    # Pipeline control parameters
    control = parser.add_argument_group("Pipeline Control Parameters")
    control.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="Resume from checkpoint if available (default: True)",
    )
    control.add_argument(
        "--no-resume", action="store_true", help="Force restart from beginning"
    )
    control.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=True,
        help="Keep intermediate files (default: True)",
    )
    control.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity level (default: INFO)",
    )
    control.add_argument(
        "--playful",
        action="store_true",
        default=True,
        help="Enable playful status messages and rich output (default: True)",
    )
    control.add_argument(
        "--no-playful",
        action="store_true",
        help="Disable playful mode for production environments",
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate command-line arguments and parameter combinations."""
    errors = []

    # Validate required directories
    if not os.path.isdir(args.genome_dir):
        errors.append(f"Genome directory does not exist: {args.genome_dir}")

    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Validate mutually exclusive parameters
    if args.protein_only and args.non_coding_only:
        errors.append("Cannot specify both --protein-only and --non-coding-only")

    if args.playful and args.no_playful:
        errors.append("Cannot specify both --playful and --no-playful")

    # Validate parameter ranges
    if not (0.0 <= args.clustering_identity <= 1.0):
        errors.append("Clustering identity must be between 0.0 and 1.0")

    if not (0.0 <= args.clustering_coverage <= 1.0):
        errors.append("Clustering coverage must be between 0.0 and 1.0")

    if not (1.0 <= args.clustering_sensitivity <= 8.0):
        errors.append("Clustering sensitivity must be between 1.0 and 8.0")

    if args.threads < 0:
        errors.append("Thread count must be non-negative")

    if args.trna_score_cutoff < 0:
        errors.append("tRNA score cutoff must be non-negative")

    if args.rrna_evalue <= 0:
        errors.append("rRNA E-value must be positive")

    if args.crispr_min_repeats < 2:
        errors.append("CRISPR minimum repeats must be at least 2")

    # Validate downstream analysis parameters
    if args.rarefaction_iterations < 1:
        errors.append("Rarefaction iterations must be positive (≥1)")

    if args.rarefaction_step_size < 1:
        errors.append("Rarefaction step size must be positive (≥1)")

    # Validate range parameters
    try:
        repeat_parts = args.crispr_repeat_length_range.split(",")
        if len(repeat_parts) != 2:
            errors.append(
                f"CRISPR repeat length range must have exactly 2 values, "
                f"got {len(repeat_parts)}: '{args.crispr_repeat_length_range}'"
            )
        else:
            repeat_min, repeat_max = map(int, repeat_parts)
            if repeat_min >= repeat_max or repeat_min < 1:
                errors.append("CRISPR repeat length range: MIN must be < MAX and > 0")
    except ValueError:
        errors.append(
            "CRISPR repeat length range must be in format MIN,MAX (e.g., 23,47)"
        )

    try:
        spacer_parts = args.crispr_spacer_length_range.split(",")
        if len(spacer_parts) != 2:
            errors.append(
                f"CRISPR spacer length range must have exactly 2 values, "
                f"got {len(spacer_parts)}: '{args.crispr_spacer_length_range}'"
            )
        else:
            spacer_min, spacer_max = map(int, spacer_parts)
            if spacer_min >= spacer_max or spacer_min < 1:
                errors.append("CRISPR spacer length range: MIN must be < MAX and > 0")
    except ValueError:
        errors.append(
            "CRISPR spacer length range must be in format MIN,MAX (e.g., 26,50)"
        )

    # Check for genomes in input directory
    genome_files = (
        list(Path(args.genome_dir).glob("*.fasta"))
        + list(Path(args.genome_dir).glob("*.fna"))
        + list(Path(args.genome_dir).glob("*.fa"))
    )

    if not genome_files:
        errors.append(
            f"No genome files found in {args.genome_dir} "
            "(expected .fasta, .fna, or .fa)"
        )

    if errors:
        print("Error: Invalid arguments:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)


def args_to_config(args: argparse.Namespace) -> PipelineConfig:
    """Convert command-line arguments to PipelineConfig object."""
    # Parse CRISPR range parameters
    repeat_parts = args.crispr_repeat_length_range.split(",")
    if len(repeat_parts) != 2:
        raise ValueError(
            f"CRISPR repeat length range must have exactly 2 values, "
            f"got {len(repeat_parts)}: '{args.crispr_repeat_length_range}'"
        )
    repeat_min, repeat_max = map(int, repeat_parts)

    spacer_parts = args.crispr_spacer_length_range.split(",")
    if len(spacer_parts) != 2:
        raise ValueError(
            f"CRISPR spacer length range must have exactly 2 values, "
            f"got {len(spacer_parts)}: '{args.crispr_spacer_length_range}'"
        )
    spacer_min, spacer_max = map(int, spacer_parts)

    # Handle resume logic
    resume = args.resume and not args.no_resume

    # Handle playful mode logic
    playful_mode = args.playful and not args.no_playful

    # Build parameter dictionaries
    prodigal_params = {
        "mode": args.prodigal_mode,
        "translation_table": args.translation_table,
    }

    trna_params = {
        "model": args.trna_model,
        "score_cutoff": args.trna_score_cutoff,
        "search_mode": args.trna_search_mode,
    }

    rrna_params = {"kingdom": args.rrna_kingdom, "evalue": args.rrna_evalue}

    crispr_params = {
        "min_repeats": args.crispr_min_repeats,
        "min_repeat_length": repeat_min,
        "max_repeat_length": repeat_max,
        "min_spacer_length": spacer_min,
        "max_spacer_length": spacer_max,
    }

    # Create configuration
    config = PipelineConfig(
        genome_dir=args.genome_dir,
        output_dir=args.output_dir,
        skip_trna=args.skip_trna,
        skip_rrna=args.skip_rrna,
        skip_crispr=args.skip_crispr,
        skip_intergenic=args.skip_intergenic,
        protein_only=args.protein_only,
        non_coding_only=args.non_coding_only,
        prodigal_params=prodigal_params,
        trna_params=trna_params,
        rrna_params=rrna_params,
        crispr_params=crispr_params,
        resume=resume,
        keep_intermediates=args.keep_intermediates,
        log_level=args.log_level,
        enable_downstream_analysis=args.enable_downstream_analysis,
        rarefaction_iterations=args.rarefaction_iterations,
        rarefaction_step_size=args.rarefaction_step_size,
        use_existing_annotations=args.use_existing_annotations,
        playful_mode=playful_mode,
    )

    return config


def setup_logging(level: str) -> logging.Logger:
    """Set up logging configuration."""
    numeric_level = getattr(logging, level.upper())
    logging.basicConfig(
        level=numeric_level,
        format="%(message)s",  # Simple format - just the message
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(__name__)


def main() -> int:
    """Main entry point for PanGenomePlus CLI."""
    try:
        # Parse arguments
        parser = create_parser()
        args = parser.parse_args()

        # Set up logging
        logger = setup_logging(args.log_level)

        # Validate arguments
        validate_arguments(args)

        # Convert to configuration
        config = args_to_config(args)
        config.validate()

        logger.info("Starting PanGenomePlus analysis")
        logger.info(f"Genome directory: {config.genome_dir}")
        logger.info(f"Output directory: {config.output_dir}")

        # Run pipeline
        stats = process_genomes(config)

        logger.info("Pipeline completed successfully!")
        logger.info(f"Processed {stats.processed_genomes} genomes")
        logger.info(f"Generated {stats.total_families} gene families")
        logger.info(
            f"Total - Core: {stats.core_families}, "
            f"Accessory: {stats.accessory_families}, "
            f"Cloud: {stats.cloud_families}"
        )

        # Per-feature-type breakdown
        logger.info("By feature type:")
        for feature_type in ["P", "I", "T", "R", "C"]:
            if stats.family_counts.get(feature_type, 0) > 0:
                core = stats.core_families_by_type.get(feature_type, 0)
                accessory = stats.accessory_families_by_type.get(feature_type, 0)
                cloud = stats.cloud_families_by_type.get(feature_type, 0)
                logger.info(
                    f"  {FeatureType.NAMES[feature_type]}: "
                    f"Core={core}, Accessory={accessory}, Cloud={cloud}"
                )

        return 0

    except KeyboardInterrupt:
        print("\\nAnalysis interrupted by user", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
