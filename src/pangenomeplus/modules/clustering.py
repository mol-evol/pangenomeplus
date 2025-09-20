"""Ultra-simplified clustering using MMseqs2 easy-cluster."""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Any
from pangenomeplus.core.types import FeatureType, GenomeFeature
from pangenomeplus.core.exceptions import ToolExecutionError

logger = logging.getLogger(__name__)


def run_clustering(
    features_by_type: Dict[FeatureType, List[GenomeFeature]],
    output_dir: Path,
    config: Dict[str, Any]
) -> Dict[FeatureType, Dict[str, List[str]]]:
    """
    Run MMseqs2 clustering using easy-cluster for maximum simplicity.

    Args:
        features_by_type: Features grouped by type
        output_dir: Output directory
        config: Configuration containing clustering parameters

    Returns:
        Clustering results by feature type
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {}

    # Extract clustering parameters from config
    clustering_config = config.get("clustering", {})
    threads = config.get("resources", {}).get("threads", 4)

    for feature_type, features in features_by_type.items():
        if not features:
            logger.info(f"No features of type {feature_type.value}, skipping clustering")
            results[feature_type] = {}
            continue

        logger.info(f"Clustering {len(features)} features of type {feature_type.value}")

        # Create feature type specific output directory
        type_output_dir = output_dir / feature_type.value
        type_output_dir.mkdir(exist_ok=True)

        # Get clustering parameters for this feature type
        type_config = get_clustering_config_for_type(feature_type, clustering_config)

        # Find existing FASTA file from annotation step
        fasta_file = find_annotation_fasta(features, feature_type, config)

        if not fasta_file:
            logger.error(f"No FASTA file found for {feature_type.value} - skipping clustering")
            results[feature_type] = {}
            continue

        # Run MMseqs2 easy-cluster
        cluster_results = run_mmseqs2_easy_cluster(
            fasta_file=fasta_file,
            feature_type=feature_type,
            output_dir=type_output_dir,
            identity=type_config.get("identity", 0.8),
            coverage=type_config.get("coverage", 0.8),
            sensitivity=type_config.get("sensitivity", 7.5),
            threads=threads,
            use_gpu=type_config.get("use_gpu", False)
        )

        results[feature_type] = cluster_results
        logger.info(f"Created {len(cluster_results)} clusters for {feature_type.value}")

    return results


def find_annotation_fasta(
    features: List[GenomeFeature],
    feature_type: FeatureType,
    config: Dict[str, Any]
) -> Path:
    """
    Find existing FASTA file from annotation step instead of recreating.

    For proteins: Use Prodigal .faa files
    For other features: Use appropriate annotation outputs
    """
    if feature_type == FeatureType.PROTEIN_CODING_GENE:
        # Look for existing Prodigal FASTA files
        # We can concatenate all genome .faa files into one
        # Get the pipeline output directory from the first feature's path context
        if features:
            # Extract output directory from feature context or use fallback
            # Will be resolved relative to current pipeline output
            annotation_dir = Path("annotation")
            # Try to find annotation directory from current working structure
            possible_dirs = [
                Path.cwd() / "annotation",
                Path.cwd().parent / "annotation",
                Path("../annotation"),
                Path("../../annotation")
            ]

            for possible_dir in possible_dirs:
                if possible_dir.exists():
                    annotation_dir = possible_dir
                    break

        if annotation_dir.exists():
            # Find first .faa file as template
            faa_files = list(annotation_dir.glob("*/*.prodigal.faa"))
            if faa_files:
                logger.info(f"Found {len(faa_files)} Prodigal FASTA files to combine")
                combined_fasta = annotation_dir.parent / "clustering" / "combined_proteins.fasta"
                combined_fasta.parent.mkdir(exist_ok=True)

                # Combine all FASTA files
                with open(combined_fasta, 'w') as outf:
                    for faa_file in faa_files:
                        with open(faa_file, 'r') as inf:
                            outf.write(inf.read())

                logger.info(f"Created combined FASTA: {combined_fasta}")
                return combined_fasta

    logger.warning(f"No existing FASTA found for {feature_type.value}")
    return None


def run_mmseqs2_easy_cluster(
    fasta_file: Path,
    feature_type: FeatureType,
    output_dir: Path,
    identity: float,
    coverage: float,
    sensitivity: float,
    threads: int,
    use_gpu: bool
) -> Dict[str, List[str]]:
    """
    Run MMseqs2 easy-cluster command directly.

    This replaces 700+ lines of manual workflow with a single command.
    """
    # Check if MMseqs2 is available
    if not shutil.which("mmseqs"):
        logger.error("MMseqs2 not found - cannot perform clustering")
        raise ToolExecutionError("mmseqs", "mmseqs not found in PATH", 127, "")

    # Create temporary directory
    tmp_dir = output_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Output prefix for easy-cluster
    output_prefix = output_dir / "clusters"

    # Build easy-cluster command
    cmd = [
        "mmseqs", "easy-cluster",
        str(fasta_file),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(identity),
        "-c", str(coverage),
        "-s", str(sensitivity),
        "--threads", str(threads)
    ]

    # Add GPU support if requested
    if use_gpu:
        cmd.extend(["--gpu", "1"])
        logger.info(f"GPU acceleration enabled for {feature_type.value}")

    logger.info(f"Running MMseqs2 easy-cluster for {feature_type.value}")
    logger.debug(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"MMseqs2 easy-cluster completed for {feature_type.value}")

        # Parse the cluster TSV output
        cluster_tsv = output_prefix.with_suffix("_cluster.tsv")
        if cluster_tsv.exists():
            return parse_easy_cluster_results(cluster_tsv)
        else:
            logger.error(f"Expected cluster output not found: {cluster_tsv}")
            return {}

    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 easy-cluster failed: {e.stderr}")
        raise ToolExecutionError("mmseqs", " ".join(cmd), e.returncode, e.stderr)


def parse_easy_cluster_results(cluster_tsv: Path) -> Dict[str, List[str]]:
    """
    Parse MMseqs2 easy-cluster TSV output.

    Format: representative_id<tab>member_id
    """
    clusters = {}

    try:
        with open(cluster_tsv, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    representative = parts[0]
                    member = parts[1]

                    if representative not in clusters:
                        clusters[representative] = []
                    clusters[representative].append(member)

    except FileNotFoundError:
        logger.error(f"Cluster TSV file not found: {cluster_tsv}")

    logger.info(f"Parsed {len(clusters)} clusters from {cluster_tsv}")
    return clusters


def get_clustering_config_for_type(
    feature_type: FeatureType,
    clustering_config: Dict[str, Any]
) -> Dict[str, Any]:
    """Get clustering configuration for specific feature type."""
    if feature_type == FeatureType.PROTEIN_CODING_GENE:
        return clustering_config.get("protein", {})
    elif feature_type == FeatureType.TRNA:
        return clustering_config.get("trna", {})
    elif feature_type == FeatureType.RRNA:
        return clustering_config.get("rrna", {})
    elif feature_type == FeatureType.INTERGENIC_REGION:
        return clustering_config.get("intergenic", {})
    else:
        return clustering_config.get("protein", {})  # Default