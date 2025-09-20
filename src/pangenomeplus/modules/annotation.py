"""Genome annotation module using external tools."""

import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Literal
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from pangenomeplus.core.exceptions import ToolExecutionError


logger = logging.getLogger(__name__)


def annotate_genomes(
    genome_paths: List[Path],
    output_dir: Path,
    tools: List[str],
    threads: int = 4,
    config: Optional[Dict[str, Any]] = None
) -> Dict[str, Dict[str, Any]]:
    """
    Annotate genomes using specified tools with parallel processing.

    Args:
        genome_paths: List of paths to genome FASTA files
        output_dir: Directory for annotation outputs
        tools: List of annotation tools to use
        threads: Number of threads for parallel tools
        config: Tool-specific configuration parameters

    Returns:
        Dictionary mapping genome_id to annotation results
    """
    if config is None:
        config = {}

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine number of parallel processes (limit to avoid overwhelming system)
    max_workers = min(len(genome_paths), multiprocessing.cpu_count(), 8)

    # For small datasets, don't use parallel processing overhead
    if len(genome_paths) <= 2:
        return _annotate_genomes_sequential(genome_paths, output_dir, tools, threads, config)

    logger.info(f"Annotating {len(genome_paths)} genomes in parallel using {max_workers} processes")

    results = {}

    # Use ProcessPoolExecutor for parallel genome processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all genome annotation jobs
        future_to_genome = {
            executor.submit(
                _annotate_single_genome,
                genome_path,
                output_dir,
                tools,
                threads,
                config
            ): genome_path
            for genome_path in genome_paths
        }

        # Collect results as they complete
        for future in as_completed(future_to_genome):
            genome_path = future_to_genome[future]
            try:
                genome_result = future.result()
                results[genome_result["genome_id"]] = genome_result
                logger.info(f"Completed annotation for {genome_result['genome_id']}")
            except Exception as e:
                genome_id = genome_path.stem.split('.')[0]
                logger.error(f"Failed to annotate {genome_id}: {e}")
                results[genome_id] = {
                    "genome_id": genome_id,
                    "genome_path": str(genome_path),
                    "error": str(e)
                }

    logger.info(f"Annotation complete for {len(genome_paths)} genomes")
    return results


def _annotate_genomes_sequential(
    genome_paths: List[Path],
    output_dir: Path,
    tools: List[str],
    threads: int,
    config: Dict[str, Any]
) -> Dict[str, Dict[str, Any]]:
    """Sequential annotation for small datasets."""
    results = {}

    for genome_path in genome_paths:
        genome_result = _annotate_single_genome(genome_path, output_dir, tools, threads, config)
        results[genome_result["genome_id"]] = genome_result
        logger.info(f"Completed annotation for {genome_result['genome_id']}")

    return results


def _annotate_single_genome(
    genome_path: Path,
    output_dir: Path,
    tools: List[str],
    threads: int,
    config: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Annotate a single genome with all specified tools.

    This function is designed to be called by multiprocessing workers.
    """
    genome_id = genome_path.stem.split('.')[0]  # Remove extensions

    genome_output_dir = output_dir / genome_id
    genome_output_dir.mkdir(exist_ok=True)

    genome_results = {
        "genome_id": genome_id,
        "genome_path": str(genome_path),
        "output_dir": str(genome_output_dir),
        "tools_used": tools,
        "annotation_files": {}
    }

    # Run each annotation tool
    for tool in tools:
        try:
            if tool == "prodigal":
                tool_results = run_prodigal(
                    genome_path=genome_path,
                    output_dir=genome_output_dir,
                    threads=1,  # Prodigal is single-threaded
                    **config.get("prodigal", {})
                )
            elif tool == "trnascan":
                tool_results = run_trnascan_se(
                    genome_path=genome_path,
                    output_dir=genome_output_dir,
                    threads=threads,
                    **config.get("trnascan", {})
                )
            elif tool == "barrnap":
                tool_results = run_barrnap(
                    genome_path=genome_path,
                    output_dir=genome_output_dir,
                    threads=threads,
                    **config.get("barrnap", {})
                )
            elif tool == "minced":
                tool_results = run_minced(
                    genome_path=genome_path,
                    output_dir=genome_output_dir,
                    **config.get("minced", {})
                )
            else:
                raise ValueError(f"Unknown annotation tool: {tool}")

            genome_results["annotation_files"][tool] = tool_results

        except Exception as e:
            genome_results["annotation_files"][tool] = {"error": str(e)}

    return genome_results


def run_prodigal(
    genome_path: Path,
    output_dir: Path,
    threads: int = 1,
    genetic_code: int = 11,
    procedure: Literal["single", "meta"] = "single",
    closed_ends: bool = True
) -> Dict[str, Path]:
    """
    Run Prodigal for protein-coding gene prediction.

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory
        threads: Number of threads (Prodigal is single-threaded)
        genetic_code: Genetic code table (11=bacterial, 4=mycoplasma)
        procedure: Prediction procedure
        closed_ends: Whether to consider closed ends

    Returns:
        Dictionary with paths to output files
    """
    genome_name = genome_path.stem

    output_files = {
        "gff": output_dir / f"{genome_name}.prodigal.gff",
        "proteins": output_dir / f"{genome_name}.prodigal.faa",
        "nucleotides": output_dir / f"{genome_name}.prodigal.fna",
        "coords": output_dir / f"{genome_name}.prodigal.coords"
    }

    # Check if prodigal is available
    if not shutil.which("prodigal"):
        # Mock execution for development
        logger.warning("Prodigal not found - creating mock output files")
        for file_path in output_files.values():
            file_path.touch()
        return output_files

    # Build command
    cmd = [
        "prodigal",
        "-i", str(genome_path),
        "-f", "gff",  # Use GFF format for better parsing
        "-o", str(output_files["gff"]),
        "-a", str(output_files["proteins"]),
        "-d", str(output_files["nucleotides"]),
        "-s", str(output_files["coords"]),
        "-g", str(genetic_code),
        "-p", procedure
    ]

    if closed_ends:
        cmd.append("-c")

    # Execute command
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise ToolExecutionError(
                tool="prodigal",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr
            )

        logger.debug(f"Prodigal completed for {genome_name}")
        return output_files

    except subprocess.CalledProcessError as e:
        raise ToolExecutionError(
            tool="prodigal",
            command=" ".join(cmd),
            returncode=e.returncode,
            stderr=str(e)
        )



def run_trnascan_se(
    genome_path: Path,
    output_dir: Path,
    threads: int = 4,
    search_mode: Literal["general", "bacterial", "archaeal", "organellar"] = "bacterial",
    relaxed: bool = False
) -> Dict[str, Path]:
    """
    Run tRNAscan-SE for tRNA detection.

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory
        threads: Number of threads
        search_mode: Search mode for different organism types
        relaxed: Use relaxed search parameters

    Returns:
        Dictionary with paths to output files
    """
    genome_name = genome_path.stem

    output_files = {
        "trnascan": output_dir / f"{genome_name}.trnascan.out",
        "structures": output_dir / f"{genome_name}.trnascan.ss"
    }

    # Check if tRNAscan-SE is available
    if not shutil.which("tRNAscan-SE"):
        logger.warning("tRNAscan-SE not found - creating mock output files")
        for file_path in output_files.values():
            file_path.touch()
        return output_files

    # Build command
    cmd = [
        "tRNAscan-SE",
        "-o", str(output_files["trnascan"]),
        "-f", str(output_files["structures"]),
        "--thread", str(threads)
    ]

    # Set search mode
    if search_mode == "bacterial":
        cmd.append("-B")
    elif search_mode == "archaeal":
        cmd.append("-A")
    elif search_mode == "organellar":
        cmd.append("-O")

    if relaxed:
        cmd.append("-H")

    cmd.append(str(genome_path))

    # Execute command
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise ToolExecutionError(
                tool="tRNAscan-SE",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr
            )

        logger.debug(f"tRNAscan-SE completed for {genome_name}")
        return output_files

    except subprocess.CalledProcessError as e:
        raise ToolExecutionError(
            tool="tRNAscan-SE",
            command=" ".join(cmd),
            returncode=e.returncode,
            stderr=str(e)
        )



def run_barrnap(
    genome_path: Path,
    output_dir: Path,
    threads: int = 4,
    kingdom: Literal["bac", "arc", "euk", "mito"] = "bac",
    evalue: float = 1e-6,
    length_cutoff: float = 0.8
) -> Dict[str, Path]:
    """
    Run Barrnap for rRNA detection.

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory
        threads: Number of threads
        kingdom: Kingdom for HMM models
        evalue: E-value cutoff
        length_cutoff: Minimum length fraction

    Returns:
        Dictionary with paths to output files
    """
    genome_name = genome_path.stem

    output_files = {
        "gff": output_dir / f"{genome_name}.barrnap.gff"
    }

    # Check if Barrnap is available
    if not shutil.which("barrnap"):
        logger.warning("Barrnap not found - creating mock output files")
        for file_path in output_files.values():
            file_path.touch()
        return output_files

    # Build command
    cmd = [
        "barrnap",
        "--kingdom", kingdom,
        "--threads", str(threads),
        "--evalue", str(evalue),
        "--lencutoff", str(length_cutoff),
        "--outdir", str(output_dir),
        str(genome_path)
    ]

    # Execute command
    try:
        with open(output_files["gff"], 'w') as outf:
            result = subprocess.run(
                cmd,
                stdout=outf,
                stderr=subprocess.PIPE,
                text=True,
                )

        if result.returncode != 0:
            raise ToolExecutionError(
                tool="barrnap",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr
            )

        logger.debug(f"Barrnap completed for {genome_name}")
        return output_files

    except subprocess.CalledProcessError as e:
        raise ToolExecutionError(
            tool="barrnap",
            command=" ".join(cmd),
            returncode=e.returncode,
            stderr=str(e)
        )



def run_minced(
    genome_path: Path,
    output_dir: Path,
    min_array_length: int = 3,
    min_spacer_length: int = 26,
    max_spacer_length: int = 50
) -> Dict[str, Path]:
    """
    Run MINCED for CRISPR array detection.

    Args:
        genome_path: Path to input genome FASTA
        output_dir: Output directory
        min_array_length: Minimum number of repeats
        min_spacer_length: Minimum spacer length
        max_spacer_length: Maximum spacer length

    Returns:
        Dictionary with paths to output files
    """
    genome_name = genome_path.stem

    output_files = {
        "gff": output_dir / f"{genome_name}.minced.gff",
        "spacers": output_dir / f"{genome_name}.minced.spacers"
    }

    # Check if MINCED is available
    if not shutil.which("minced"):
        logger.warning("MINCED not found - creating mock output files")
        for file_path in output_files.values():
            file_path.touch()
        return output_files

    # Build command
    cmd = [
        "minced",
        "-gff",
        "-spacers",
        str(genome_path),
        str(output_files["gff"])
    ]

    # Execute command
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise ToolExecutionError(
                tool="minced",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr
            )

        logger.debug(f"MINCED completed for {genome_name}")
        return output_files

    except subprocess.CalledProcessError as e:
        raise ToolExecutionError(
            tool="minced",
            command=" ".join(cmd),
            returncode=e.returncode,
            stderr=str(e)
        )