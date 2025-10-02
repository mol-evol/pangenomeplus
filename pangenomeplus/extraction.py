"""External tool integration and feature extraction functions.

This module provides simple function-based wrappers around external bioinformatics
tools following KISS principles. Each function handles tool execution and output
parsing with minimal abstractions.

External tools integrated:
- Prodigal: Protein-coding gene prediction
- tRNAscan-SE: tRNA detection
- Barrnap: rRNA detection
- MINCED: CRISPR detection

All functions follow the pattern:
1. Execute external tool with appropriate parameters
2. Parse native tool output format
3. Return Feature objects with compact IDs assigned
"""

import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List

from Bio import SeqIO

from .compact_ids import CompactIDManager
from .constants import MIN_INTERGENIC_LENGTH
from .core import Feature


class ExtractionError(Exception):
    """Exception raised during feature extraction."""

    pass


def run_prodigal(
    genome_file: str, output_dir: str, mode: str = "single", translation_table: int = 11
) -> str:
    """Run Prodigal gene prediction on a genome file.

    Args:
        genome_file: Path to input FASTA genome file
        output_dir: Directory to write output files
        mode: Prodigal mode ('single' or 'meta')
        translation_table: Genetic code table number

    Returns:
        Path to generated GFF3 file

    Raises:
        ExtractionError: If Prodigal execution fails
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(genome_file).stem
    gff_file = os.path.join(output_dir, f"{base_name}.gff")
    faa_file = os.path.join(output_dir, f"{base_name}.faa")
    fna_file = os.path.join(output_dir, f"{base_name}.fna")

    cmd = [
        "prodigal",
        "-i",
        genome_file,
        "-o",
        gff_file,
        "-a",
        faa_file,
        "-d",
        fna_file,
        "-f",
        "gff",
        "-p",
        mode,
        "-g",
        str(translation_table),
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise ExtractionError(f"Prodigal failed: {e.stderr}")
    except FileNotFoundError:
        raise ExtractionError("Prodigal not found in PATH")

    if not os.path.exists(gff_file):
        raise ExtractionError(
            f"Prodigal did not generate expected GFF file: {gff_file}"
        )

    return gff_file


def run_trnascan(
    genome_file: str,
    output_dir: str,
    model: str = "bacteria",
    score_cutoff: float = 20.0,
    search_mode: str = "normal",
) -> str:
    """Run tRNAscan-SE for tRNA detection.

    Args:
        genome_file: Path to input FASTA genome file
        output_dir: Directory to write output files
        model: Organism model ('bacteria', 'archaea', 'eukaryota')
        score_cutoff: Score threshold for tRNA detection

    Returns:
        Path to generated tab-delimited output file

    Raises:
        ExtractionError: If tRNAscan-SE execution fails
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(genome_file).stem
    output_file = os.path.join(output_dir, f"{base_name}_trna.txt")

    # Map model names to tRNAscan-SE parameters
    model_flags = {
        "bacteria": ["-B"],
        "archaea": ["-A"],
        "eukaryota": ["-E"],
        "mitochondrial": ["-M", "mamm"],
    }

    if model not in model_flags:
        raise ExtractionError(f"Unknown tRNA model: {model}")

    cmd = [
        "trnascan-se",
        *model_flags[model],
        "-X",
        str(score_cutoff),
        "-o",
        output_file,
        genome_file,
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise ExtractionError(f"tRNAscan-SE failed: {e.stderr}")
    except FileNotFoundError:
        raise ExtractionError("tRNAscan-SE not found in PATH")

    if not os.path.exists(output_file):
        raise ExtractionError(
            f"tRNAscan-SE did not generate expected output: {output_file}"
        )

    return output_file


def run_barrnap(
    genome_file: str, output_dir: str, kingdom: str = "bac", evalue: float = 1e-6
) -> str:
    """Run Barrnap for rRNA detection.

    Args:
        genome_file: Path to input FASTA genome file
        output_dir: Directory to write output files
        kingdom: Kingdom model ('bac', 'arc', 'euk', 'mito')
        evalue: E-value threshold

    Returns:
        Path to generated GFF3 file

    Raises:
        ExtractionError: If Barrnap execution fails
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(genome_file).stem
    output_file = os.path.join(output_dir, f"{base_name}_rrna.gff")

    cmd = [
        "barrnap",
        "--kingdom",
        kingdom,
        "--evalue",
        str(evalue),
        "--outseq",
        "/dev/null",  # We don't need sequence output
        genome_file,
    ]

    try:
        with open(output_file, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise ExtractionError(f"Barrnap failed: {e.stderr}")
    except FileNotFoundError:
        raise ExtractionError("Barrnap not found in PATH")

    if not os.path.exists(output_file):
        raise ExtractionError(
            f"Barrnap did not generate expected output: {output_file}"
        )

    return output_file


def run_minced(
    genome_file: str,
    output_dir: str,
    min_repeats: int = 3,
    min_spacer_length: int = 26,
    max_spacer_length: int = 50,
) -> str:
    """Run MINCED for CRISPR detection.

    Args:
        genome_file: Path to input FASTA genome file
        output_dir: Directory to write output files
        min_repeats: Minimum number of repeats for CRISPR array
        min_spacer_length: Minimum spacer length
        max_spacer_length: Maximum spacer length

    Returns:
        Path to generated text output file

    Raises:
        ExtractionError: If MINCED execution fails
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(genome_file).stem
    output_file = os.path.join(output_dir, f"{base_name}_crispr.txt")

    cmd = [
        "minced",
        "-minNR",
        str(min_repeats),
        "-minSL",
        str(min_spacer_length),
        "-maxSL",
        str(max_spacer_length),
        genome_file,
        output_file,
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise ExtractionError(f"MINCED failed: {e.stderr}")
    except FileNotFoundError:
        raise ExtractionError("MINCED not found in PATH")

    # MINCED always creates output file, even if no CRISPRs found
    if not os.path.exists(output_file):
        raise ExtractionError(f"MINCED did not generate expected output: {output_file}")

    return output_file


def parse_minced_output(
    minced_file: str, genome_file: str, genome_id: str, id_manager: CompactIDManager
) -> List[Feature]:
    """Parse MINCED text output and create Feature objects for CRISPR spacers.

    Args:
        minced_file: Path to MINCED text output file
        genome_file: Path to original genome FASTA file
        genome_id: Identifier for this genome
        id_manager: CompactIDManager for ID assignment

    Returns:
        List of Feature objects for CRISPR spacer sequences

    Raises:
        ExtractionError: If parsing fails
    """
    features = []

    # Validate genome file exists and is readable
    try:
        genome_records = list(SeqIO.parse(genome_file, "fasta"))
        if not genome_records:
            raise ValueError("No sequences found in genome file")
    except Exception as e:
        raise ExtractionError(f"Failed to load genome: {e}")

    try:
        with open(minced_file, "r") as f:
            lines = f.readlines()

        current_contig = None
        current_crispr_start = None
        spacer_count = 0

        for line in lines:
            line = line.strip()

            # Parse sequence header to get contig name
            if line.startswith("Sequence '") and "' (" in line:
                current_contig = line.split("'")[1]
                continue

            # Parse CRISPR array header
            if line.startswith("CRISPR ") and "Range:" in line:
                # Extract range: "CRISPR 1   Range: 50 - 120"
                range_part = line.split("Range:")[1].strip()
                start, end = map(int, range_part.split(" - "))
                current_crispr_start = start
                continue

            # Parse spacer lines (contain both repeat and spacer positions)
            if current_contig and current_crispr_start and line and line[0].isdigit():
                parts = line.split("\t")
                if len(parts) >= 5:  # Position, tabs, repeat, tabs, spacer
                    try:
                        position = int(parts[0])
                        repeat_seq = parts[2].strip() if len(parts) > 2 else ""
                        spacer_seq = parts[4].strip() if len(parts) > 4 else ""

                        # Only process lines that have spacer sequences
                        if (
                            spacer_seq
                            and spacer_seq != ""
                            and spacer_seq not in ["------"]
                        ):
                            # Calculate spacer coordinates
                            # Spacer starts after the repeat at this position
                            repeat_len = len(repeat_seq)
                            spacer_start = position + repeat_len
                            spacer_end = spacer_start + len(spacer_seq) - 1

                            # Generate compact ID
                            compact_id = id_manager.generate_compact_id("C")

                            crispr_feature = Feature(
                                compact_id=compact_id,
                                genome_id=genome_id,
                                contig=current_contig,
                                start=spacer_start,
                                end=spacer_end,
                                strand=".",  # CRISPR spacers don't have inherent strand
                                sequence=spacer_seq,
                                feature_type="C",
                                original_id=(f"crispr_spacer_{spacer_count + 1}"),
                                metadata={
                                    "product": "CRISPR_spacer",
                                    "array_start": current_crispr_start,
                                    "repeat_sequence": repeat_seq,
                                    "spacer_position": position,
                                },
                            )

                            features.append(crispr_feature)
                            id_manager.register_feature(crispr_feature)
                            spacer_count += 1

                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue

    except Exception as e:
        raise ExtractionError(f"Failed to parse MINCED output {minced_file}: {e}")

    return features


def parse_prodigal_gff(
    gff_file: str, genome_file: str, genome_id: str, id_manager: CompactIDManager
) -> List[Feature]:
    """Parse Prodigal GFF3 output and create Feature objects.

    Uses linear time O(n) algorithm - GFF3 features are already coordinate-sorted.

    Args:
        gff_file: Path to Prodigal GFF3 output
        genome_file: Path to original genome FASTA file
        genome_id: Identifier for this genome
        id_manager: CompactIDManager for ID assignment

    Returns:
        List of Feature objects for protein-coding genes

    Raises:
        ExtractionError: If parsing fails
    """
    features = []

    # Load genome sequences - handle multi-contig genomes
    try:
        genome_records = list(SeqIO.parse(genome_file, "fasta"))
        if not genome_records:
            raise ValueError("No sequences found")

        # Create contig lookup for multi-contig genomes
        contig_sequences = {record.id: str(record.seq) for record in genome_records}
    except Exception as e:
        raise ExtractionError(f"Failed to load genome sequence: {e}")

    try:
        with open(gff_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) != 9:
                    continue

                # Parse GFF3 fields
                contig = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]

                # Only process CDS features (protein-coding genes)
                if feature_type != "CDS":
                    continue

                # Extract sequence from appropriate contig
                try:
                    if contig not in contig_sequences:
                        raise ExtractionError(f"Contig {contig} not found in genome")

                    contig_seq = contig_sequences[contig]
                    sequence = contig_seq[start - 1 : end]
                    if strand == "-":
                        # Reverse complement for minus strand
                        complement_map = str.maketrans("ATCG", "TAGC")
                        sequence = sequence.translate(complement_map)[::-1]
                except IndexError:
                    contig_len = len(contig_sequences.get(contig, ""))
                    raise ExtractionError(
                        f"Invalid coordinates: {start}-{end} for contig {contig} "
                        f"length {contig_len}"
                    )

                # Generate compact ID
                compact_id = id_manager.generate_compact_id("P")

                # Parse attributes for metadata
                metadata = {}
                original_id = ""
                for attr in attributes.split(";"):
                    if "=" in attr:
                        key, value = attr.split("=", 1)
                        metadata[key] = value
                        if key == "ID":
                            original_id = value

                feature = Feature(
                    compact_id=compact_id,
                    genome_id=genome_id,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    sequence=sequence,
                    feature_type="P",
                    original_id=original_id,
                    metadata=metadata,
                )

                features.append(feature)

                # Register in ID manager
                id_manager.register_feature(feature)

    except Exception as e:
        raise ExtractionError(f"Failed to parse Prodigal GFF: {e}")

    return features


def calculate_intergenic_regions(
    features: List[Feature],
    genome_file: str,
    genome_id: str,
    id_manager: CompactIDManager,
    min_length: int = MIN_INTERGENIC_LENGTH,
) -> List[Feature]:
    """Calculate intergenic regions between annotated features.

    Uses linear time O(n) algorithm - assumes features are coordinate-sorted.

    Args:
        features: List of annotated features (sorted by start coordinate)
        genome_file: Path to genome FASTA file
        genome_id: Identifier for this genome
        id_manager: CompactIDManager for ID assignment
        min_length: Minimum intergenic region length

    Returns:
        List of Feature objects for intergenic regions

    Raises:
        ExtractionError: If calculation fails
    """
    intergenic_features = []

    # Load genome sequences - handle multi-contig genomes
    try:
        genome_records = list(SeqIO.parse(genome_file, "fasta"))
        if not genome_records:
            raise ValueError("No sequences found")

        # For intergenic calculation, use first/largest contig
        genome_record = max(genome_records, key=lambda x: len(x.seq))
        genome_seq = str(genome_record.seq)
        contig_name = genome_record.id
    except Exception as e:
        raise ExtractionError(f"Failed to load genome sequence: {e}")

    # Features already coordinate-sorted from GFF3 parsing - use directly for O(n) algorithm
    # No need for O(n log n) sorting operation per CLAUDE.md performance guidelines

    # Linear pass through features to find gaps - O(n)
    prev_end = 1  # Start from position 1 (1-based coordinates)

    for feature in features:
        # Check for intergenic region before this feature
        gap_start = prev_end
        gap_end = feature.start - 1
        gap_length = gap_end - gap_start + 1

        if gap_length >= min_length:
            # Extract intergenic sequence
            try:
                sequence = genome_seq[
                    gap_start - 1 : gap_end
                ]  # Convert to 0-based for slicing
            except IndexError:
                continue  # Skip invalid coordinates

            # Generate compact ID
            compact_id = id_manager.generate_compact_id("I")

            intergenic_feature = Feature(
                compact_id=compact_id,
                genome_id=genome_id,
                contig=contig_name,
                start=gap_start,
                end=gap_end,
                strand=".",
                sequence=sequence,
                feature_type="I",
                original_id=f"intergenic_{gap_start}_{gap_end}",
                metadata={"product": "intergenic_region"},
            )

            intergenic_features.append(intergenic_feature)
            id_manager.register_feature(intergenic_feature)

        # Update end position for next iteration
        prev_end = max(prev_end, feature.end + 1)

    # Check for final intergenic region after last feature
    if prev_end <= len(genome_seq):
        gap_start = prev_end
        gap_end = len(genome_seq)
        gap_length = gap_end - gap_start + 1

        if gap_length >= min_length:
            try:
                sequence = genome_seq[gap_start - 1 : gap_end]
                compact_id = id_manager.generate_compact_id("I")

                final_feature = Feature(
                    compact_id=compact_id,
                    genome_id=genome_id,
                    contig=contig_name,
                    start=gap_start,
                    end=gap_end,
                    strand=".",
                    sequence=sequence,
                    feature_type="I",
                    original_id=f"intergenic_{gap_start}_{gap_end}",
                    metadata={"product": "intergenic_region"},
                )

                intergenic_features.append(final_feature)
                id_manager.register_feature(final_feature)
            except IndexError:
                pass  # Skip if coordinates are invalid

    return intergenic_features


def parse_trnascan_output(
    trnascan_file: str, genome_file: str, genome_id: str, id_manager: CompactIDManager
) -> List[Feature]:
    """Parse tRNAscan-SE output and create Feature objects.

    Args:
        trnascan_file: Path to tRNAscan-SE output file
        genome_file: Path to genome FASTA file
        genome_id: Genome identifier
        id_manager: CompactIDManager for ID assignment

    Returns:
        List of tRNA Feature objects
    """
    tRNA_features = []

    # Load genome sequences (handle multi-contig genomes)
    genome_seqs = {record.id: record for record in SeqIO.parse(genome_file, "fasta")}

    with open(trnascan_file, "r") as f:
        for line in f:
            line = line.strip()
            if (
                not line
                or line.startswith("Sequence")
                or line.startswith("Name")
                or line.startswith("---")
            ):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            # Strip whitespace from all fields (tRNAscan-SE adds trailing spaces)
            contig = fields[0].strip()
            start_raw = int(fields[2].strip())
            end_raw = int(fields[3].strip())
            trna_type = fields[4].strip()
            anticodon = fields[5].strip()
            score = float(fields[8].strip())

            # Determine strand and normalize coordinates
            # tRNAscan-SE outputs start > end for reverse strand
            if start_raw > end_raw:
                start = end_raw
                end = start_raw
                strand = "-"
            else:
                start = start_raw
                end = end_raw
                strand = "+"

            # Extract sequence from correct contig
            if contig not in genome_seqs:
                continue  # Skip if contig not found
            sequence = str(genome_seqs[contig].seq[start - 1 : end])

            # Create compact ID
            compact_id = id_manager.generate_compact_id("T")

            # Create feature
            feature = Feature(
                compact_id=compact_id,
                genome_id=genome_id,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                sequence=sequence,
                feature_type="T",
                original_id=f"tRNA-{trna_type}-{anticodon}",
                metadata={
                    "product": f"tRNA-{trna_type}",
                    "anticodon": anticodon,
                    "score": score,
                    "trna_type": trna_type,
                },
            )

            tRNA_features.append(feature)
            id_manager.register_feature(feature)

    return tRNA_features


def parse_barrnap_gff(
    barrnap_file: str, genome_file: str, genome_id: str, id_manager: CompactIDManager
) -> List[Feature]:
    """Parse Barrnap GFF output and create Feature objects.

    Args:
        barrnap_file: Path to Barrnap GFF output file
        genome_file: Path to genome FASTA file
        genome_id: Genome identifier
        id_manager: CompactIDManager for ID assignment

    Returns:
        List of rRNA Feature objects
    """
    rRNA_features = []

    # Load genome sequences (handle multi-contig genomes)
    genome_seqs = {record.id: record for record in SeqIO.parse(genome_file, "fasta")}

    with open(barrnap_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            contig = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            attributes = fields[8]

            # Parse attributes for Name and product
            name = ""
            product = ""
            for attr in attributes.split(";"):
                if attr.startswith("Name="):
                    name = attr.split("=")[1]
                elif attr.startswith("product="):
                    product = attr.split("=")[1]

            # Extract sequence from correct contig
            if contig not in genome_seqs:
                continue  # Skip if contig not found

            if strand == "-":
                sequence = str(
                    genome_seqs[contig].seq[start - 1 : end].reverse_complement()
                )
            else:
                sequence = str(genome_seqs[contig].seq[start - 1 : end])

            # Create compact ID
            compact_id = id_manager.generate_compact_id("R")

            # Create feature
            feature = Feature(
                compact_id=compact_id,
                genome_id=genome_id,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                sequence=sequence,
                feature_type="R",
                original_id=name,
                metadata={"product": product, "score": score, "rRNA_type": name},
            )

            rRNA_features.append(feature)
            id_manager.register_feature(feature)

    return rRNA_features


def extract_genome_features(
    genome_file: str,
    genome_id: str,
    output_dir: str,
    id_manager: CompactIDManager,
    skip_trna: bool = False,
    skip_rrna: bool = False,
    skip_crispr: bool = False,
    skip_intergenic: bool = False,
    use_existing_annotations: bool = False,
    **tool_params: Any,
) -> Dict[str, List[Feature]]:
    """Extract all features from a genome using external tools.

    Orchestrates external tool execution and feature extraction following
    the complete pipeline flow specified in CLAUDE.md.

    Args:
        genome_file: Path to genome FASTA file
        genome_id: Identifier for this genome
        output_dir: Directory for intermediate files
        id_manager: CompactIDManager for ID assignment
        skip_trna: Skip tRNA detection
        skip_rrna: Skip rRNA detection
        skip_crispr: Skip CRISPR detection
        skip_intergenic: Skip intergenic region calculation
        **tool_params: Parameters for external tools

    Returns:
        Dictionary mapping feature types to lists of Feature objects

    Raises:
        ExtractionError: If extraction fails
    """
    all_features: Dict[str, List[Feature]] = {
        "proteins": [],
        "intergenic": [],
        "tRNAs": [],
        "rRNAs": [],
        "CRISPR": [],
    }

    genome_output_dir = os.path.join(output_dir, genome_id)
    os.makedirs(genome_output_dir, exist_ok=True)

    # 1. Run Prodigal for protein-coding genes (always required for intergenic calculation)
    try:
        # Check for existing GFF3 file if requested
        if use_existing_annotations:
            genome_base = Path(genome_file).stem
            potential_gff = Path(genome_file).parent / f"{genome_base}.gff"
            potential_gff3 = Path(genome_file).parent / f"{genome_base}.gff3"

            existing_gff = None
            if potential_gff.exists():
                existing_gff = str(potential_gff)
            elif potential_gff3.exists():
                existing_gff = str(potential_gff3)

            if existing_gff:
                try:
                    protein_features = parse_prodigal_gff(
                        existing_gff, genome_file, genome_id, id_manager
                    )
                    all_features["proteins"] = protein_features
                except (IOError, OSError, ExtractionError) as e:
                    # Fall back to Prodigal on file or parsing errors
                    logger = logging.getLogger(__name__)
                    logger.warning(f"Failed to use existing GFF {existing_gff}: {e}")
                    logger.info("Falling back to Prodigal for gene prediction")
                    gff_file = run_prodigal(
                        genome_file,
                        genome_output_dir,
                        **tool_params.get("prodigal", {}),
                    )
                    protein_features = parse_prodigal_gff(
                        gff_file, genome_file, genome_id, id_manager
                    )
                    all_features["proteins"] = protein_features
            else:
                # No existing GFF found, run Prodigal
                gff_file = run_prodigal(
                    genome_file, genome_output_dir, **tool_params.get("prodigal", {})
                )
                protein_features = parse_prodigal_gff(
                    gff_file, genome_file, genome_id, id_manager
                )
                all_features["proteins"] = protein_features
        else:
            # Standard Prodigal run
            gff_file = run_prodigal(
                genome_file, genome_output_dir, **tool_params.get("prodigal", {})
            )
            protein_features = parse_prodigal_gff(
                gff_file, genome_file, genome_id, id_manager
            )
            all_features["proteins"] = protein_features
    except Exception as e:
        raise ExtractionError(f"Protein extraction failed for {genome_id}: {e}")

    # 2. Calculate intergenic regions (uses protein features as boundaries)
    if not skip_intergenic:
        try:
            intergenic_features = calculate_intergenic_regions(
                protein_features, genome_file, genome_id, id_manager
            )
            all_features["intergenic"] = intergenic_features
        except Exception as e:
            raise ExtractionError(f"Intergenic extraction failed for {genome_id}: {e}")

    # 3. CRISPR extraction
    if not skip_crispr:
        try:
            # Extract CRISPR parameters from tool_params
            crispr_params = {
                "min_repeats": tool_params.get("crispr_min_repeats", 3),
                "min_spacer_length": tool_params.get("crispr_min_spacer_length", 26),
                "max_spacer_length": tool_params.get("crispr_max_spacer_length", 50),
            }

            # Run MINCED to detect CRISPR arrays
            minced_output = run_minced(
                genome_file=genome_file, output_dir=output_dir, **crispr_params
            )

            # Parse MINCED output and create features
            crispr_features = parse_minced_output(
                minced_file=minced_output,
                genome_file=genome_file,
                genome_id=genome_id,
                id_manager=id_manager,
            )
            all_features["CRISPR"] = crispr_features

        except Exception as e:
            raise ExtractionError(f"CRISPR extraction failed for {genome_id}: {e}")

    # 4. tRNA detection with tRNAscan-SE
    if not skip_trna:
        try:
            trna_params = tool_params.get("trna", {})

            # Run tRNAscan-SE to detect tRNAs
            trnascan_output = run_trnascan(
                genome_file=genome_file, output_dir=genome_output_dir, **trna_params
            )

            # Parse tRNAscan-SE output and create features
            trna_features = parse_trnascan_output(
                trnascan_file=trnascan_output,
                genome_file=genome_file,
                genome_id=genome_id,
                id_manager=id_manager,
            )
            all_features["tRNAs"] = trna_features

        except Exception as e:
            raise ExtractionError(f"tRNA extraction failed for {genome_id}: {e}")

    # 5. rRNA detection with Barrnap
    if not skip_rrna:
        try:
            rrna_params = tool_params.get("rrna", {})

            # Run Barrnap to detect rRNAs
            barrnap_output = run_barrnap(
                genome_file=genome_file, output_dir=genome_output_dir, **rrna_params
            )

            # Parse Barrnap output and create features
            rrna_features = parse_barrnap_gff(
                barrnap_file=barrnap_output,
                genome_file=genome_file,
                genome_id=genome_id,
                id_manager=id_manager,
            )
            all_features["rRNAs"] = rrna_features

        except Exception as e:
            raise ExtractionError(f"rRNA extraction failed for {genome_id}: {e}")

    return all_features
