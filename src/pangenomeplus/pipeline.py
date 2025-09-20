"""Main pipeline orchestration module."""

import logging
from pathlib import Path
from typing import Dict, Any, Union, Optional, Literal, List
import psutil
import time

from pangenomeplus.core.types import ValidationResult
from pangenomeplus.core.exceptions import PipelineError, ValidationError
from pangenomeplus.utils.config import validate_configuration_schema
from pangenomeplus.modules.annotation import annotate_genomes
from pangenomeplus.modules.extraction import extract_features_from_annotations
from pangenomeplus.modules.clustering import run_clustering
from pangenomeplus.modules.families import assign_gene_families
from pangenomeplus.modules.output import generate_outputs


def run_pangenome_pipeline(
    genomes_dir: Union[str, Path],
    output_dir: Union[str, Path],
    config: Dict[str, Any],
    resume: bool = False,
    stop_after: Optional[Literal["annotation", "clustering", "families", "analysis"]] = None,
    validate_inputs: bool = True,
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = "INFO"
) -> Dict[str, Any]:
    """
    Simple pangenome analysis pipeline.

    Args:
        genomes_dir: Directory containing input genome FASTA files
        output_dir: Directory for pipeline outputs
        config: Complete pipeline configuration dictionary
        resume: Whether to resume from existing checkpoints (simplified)
        stop_after: Optional stage to stop after
        validate_inputs: Whether to validate inputs before processing
        log_level: Logging verbosity level

    Returns:
        Dictionary containing pipeline results and metadata
    """
    # Convert paths
    genomes_dir = Path(genomes_dir)
    output_dir = Path(output_dir)

    # Create output directory structure first
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)

    # Setup logging
    setup_logging(log_level, output_dir / "pipeline.log")
    logger = logging.getLogger(__name__)

    logger.info("Starting pangenome pipeline")
    logger.info(f"Input genomes: {genomes_dir}")
    logger.info(f"Output directory: {output_dir}")

    start_time = time.time()
    results = {
        "start_time": start_time,
        "pipeline_version": "1.0.0",
        "config": config
    }

    try:
        # Validation stage
        if validate_inputs:
            logger.info("Validating inputs...")
            validation_result = validate_pipeline_inputs(genomes_dir, config)
            if not validation_result.is_valid:
                raise ValidationError(f"Input validation failed: {validation_result.errors}")
            logger.info("Input validation passed")

        # No strategy selection needed - MMseqs2 handles scaling automatically
        logger.info("Dataset will be processed using MMseqs2 easy-cluster")

        # Annotation stage
        logger.info("Starting genome annotation...")
        annotation_dir = output_dir / "annotation"
        annotation_dir.mkdir(exist_ok=True)

        genome_files = list_genome_files(genomes_dir)
        logger.info(f"Found {len(genome_files)} genome files")

        annotation_results = annotate_genomes(
            genome_paths=genome_files,
            output_dir=annotation_dir,
            tools=config.get("annotation", {}).get("tools", ["prodigal", "trnascan", "barrnap"]),
            threads=config.get("resources", {}).get("threads", 4),
            config=config.get("annotation", {})
        )

        results["n_genomes"] = len(genome_files)
        results["annotation_results"] = annotation_results

        if stop_after == "annotation":
            logger.info("Stopping after annotation stage as requested")
            return results

        # Feature extraction stage
        logger.info("Extracting features from annotations...")
        features_dir = output_dir / "features"
        features_dir.mkdir(exist_ok=True)

        # Get intergenic extraction settings from config
        intergenic_config = config.get("annotation", {}).get("intergenic", {})
        extract_intergenic = intergenic_config.get("extract", True)
        min_intergenic_length = intergenic_config.get("min_length", 50)

        features, sequence_indices = extract_features_from_annotations(
            annotation_results=annotation_results,
            output_dir=features_dir,
            extract_intergenic=extract_intergenic,
            min_intergenic_length=min_intergenic_length
        )

        results["n_features"] = len(features)
        logger.info(f"Extracted {len(features)} features")

        if stop_after == "clustering":
            logger.info("Stopping after feature extraction as requested")
            return results

        # Clustering stage
        logger.info("Starting clustering...")
        clustering_dir = output_dir / "clustering"
        clustering_dir.mkdir(exist_ok=True)

        # Group features by type
        features_by_type = {}
        for feature in features:
            if feature.feature_type not in features_by_type:
                features_by_type[feature.feature_type] = []
            features_by_type[feature.feature_type].append(feature)

        logger.info(
            f"Features grouped by type: {[(ft.value, len(feats)) for ft, feats in features_by_type.items()]}"
        )

        # Execute simplified clustering using MMseqs2 easy-cluster
        cluster_results = run_clustering(
            features_by_type=features_by_type,
            output_dir=clustering_dir,
            config=config
        )

        results["cluster_results"] = cluster_results

        if stop_after == "families":
            logger.info("Stopping after clustering as requested")
            return results

        # Gene family assignment
        logger.info("Assigning gene families...")
        families_dir = output_dir / "families"
        families_dir.mkdir(exist_ok=True)

        # Create ID mappings
        id_mappings = {feature.pangenomeplus_id: feature.original_id for feature in features}

        gene_families, family_members = assign_gene_families(
            cluster_results=cluster_results,
            id_mappings=id_mappings,
            output_dir=families_dir
        )

        results["n_families"] = len(set(gene_families.values()))
        logger.info(f"Created {results['n_families']} gene families")

        # Calculate runtime for reporting (before output generation)
        end_time = time.time()
        results["end_time"] = end_time
        results["runtime_seconds"] = end_time - start_time
        results["runtime_formatted"] = f"{(end_time - start_time) / 60:.1f} minutes"

        # Output generation
        logger.info("Generating outputs...")
        output_formats_dir = output_dir / "outputs"
        output_formats_dir.mkdir(exist_ok=True)

        output_files = generate_outputs(
            gene_families=gene_families,
            family_members=family_members,
            features=features,
            id_mappings={"pangenomeplus_to_original": id_mappings},
            output_dir=output_formats_dir,
            formats=config.get("output", {}).get("formats", ["transformer", "presence_absence"]),
            config=config,
            pipeline_results=results,
            annotation_results=annotation_results,
            cluster_results=cluster_results
        )

        results["output_files"] = {k: str(v) for k, v in output_files.items()}

        logger.info(f"Pipeline completed successfully in {results['runtime_formatted']}")
        logger.info(f"Processed {results['n_genomes']} genomes")
        logger.info(f"Generated {results['n_families']} gene families")

        return results

    except Exception as e:
        end_time = time.time()
        results["end_time"] = end_time
        results["runtime_seconds"] = end_time - start_time
        results["error"] = str(e)

        logger.error(f"Pipeline failed after {(end_time - start_time) / 60:.1f} minutes: {e}")
        raise


def analyze_dataset_and_choose_strategy(
    genomes_dir: Union[str, Path],
    sample_size: int = 10
) -> Dict[str, Any]:
    """
    Analyze dataset characteristics for informational purposes.

    MMseqs2 handles scaling automatically, so no strategy selection needed.

    Args:
        genomes_dir: Directory containing genome files
        sample_size: Number of genomes to sample for analysis

    Returns:
        Dictionary with dataset information
    """
    genomes_dir = Path(genomes_dir)
    genome_files = list_genome_files(genomes_dir)
    n_genomes = len(genome_files)

    memory_available = psutil.virtual_memory().total // (1024**3)  # Convert to GB

    return {
        "n_genomes": n_genomes,
        "memory_available_gb": memory_available,
        "description": (
            f"Dataset with {n_genomes} genomes - MMseqs2 will handle scaling automatically"
        )
    }


def validate_pipeline_inputs(
    genomes_dir: Union[str, Path],
    config: Dict[str, Any],
    check_tools: bool = True
) -> ValidationResult:
    """
    Comprehensive validation of pipeline inputs.

    Args:
        genomes_dir: Directory containing genome files
        config: Pipeline configuration
        check_tools: Whether to verify external tool availability

    Returns:
        ValidationResult with validation status and details
    """
    errors = []
    warnings = []
    details = {}

    # Validate genome directory
    genomes_dir = Path(genomes_dir)
    if not genomes_dir.exists():
        errors.append(f"Genomes directory does not exist: {genomes_dir}")
    elif not genomes_dir.is_dir():
        errors.append(f"Genomes path is not a directory: {genomes_dir}")
    else:
        genome_files = list_genome_files(genomes_dir)
        if not genome_files:
            errors.append(f"No genome files found in {genomes_dir}")
        else:
            details["n_genomes"] = len(genome_files)

    # Validate configuration
    config_validation = validate_configuration_schema(config)
    if not config_validation.is_valid:
        errors.extend(config_validation.errors)

    # Check system resources
    memory_gb = psutil.virtual_memory().total // (1024**3)
    if memory_gb < 8:
        warnings.append(f"Low system memory: {memory_gb}GB (recommended: 32GB+)")
    details["system_memory_gb"] = memory_gb

    # Tool validation (simplified for now)
    if check_tools:
        # This would check for external tools like mmseqs2, prodigal, etc.
        # For now, just log that we would check
        details["tool_check"] = "skipped in development mode"

    is_valid = len(errors) == 0

    return ValidationResult(
        is_valid=is_valid,
        errors=errors,
        warnings=warnings,
        details=details
    )


def list_genome_files(genomes_dir: Path) -> List[Path]:
    """
    List all genome FASTA files in the directory.

    Args:
        genomes_dir: Directory to search

    Returns:
        List of genome file paths
    """
    extensions = [".fasta", ".fa", ".fna", ".fas"]
    genome_files = []

    for ext in extensions:
        genome_files.extend(genomes_dir.glob(f"*{ext}"))
        genome_files.extend(genomes_dir.glob(f"*{ext}.gz"))

    return sorted(genome_files)


def setup_logging(level: str, log_file: Optional[Path] = None) -> None:
    """
    Setup logging configuration.

    Args:
        level: Logging level
        log_file: Optional log file path
    """
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file) if log_file else logging.NullHandler()
        ]
    )