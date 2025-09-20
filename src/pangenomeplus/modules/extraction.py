"""Feature extraction and standardization module."""

import logging
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pangenomeplus.core.types import GenomeFeature, FeatureType
from pangenomeplus.core.id_generation import generate_pangenomeplus_id


logger = logging.getLogger(__name__)


def _extract_features_from_single_genome(
    genome_id: str,
    genome_data: Dict[str, Any],
    extract_intergenic: bool = True,
    min_intergenic_length: int = 50
) -> List[GenomeFeature]:
    """
    Extract features from a single genome's annotation data.
    This function is designed to be called by multiprocessing workers.

    Args:
        genome_id: Identifier for the genome
        genome_data: Annotation data for the genome
        extract_intergenic: Whether to extract intergenic regions
        min_intergenic_length: Minimum length for intergenic regions (bp)

    Returns:
        List of GenomeFeature objects for this genome
    """
    genome_features = []
    annotation_files = genome_data.get("annotation_files", {})

    # Extract protein-coding genes from Prodigal
    if "prodigal" in annotation_files:
        prodigal_files = annotation_files["prodigal"]
        if "error" not in prodigal_files:
            try:
                protein_features = parse_prodigal_results(
                    gff_file=Path(prodigal_files.get("gff", "")),
                    protein_file=Path(prodigal_files.get("proteins", "")),
                    nucleotide_file=Path(prodigal_files.get("nucleotides", "")),
                    genome_id=genome_id
                )
                genome_features.extend(protein_features)
                logger.debug(f"Extracted {len(protein_features)} protein features from {genome_id}")
            except Exception as e:
                logger.error(f"Error parsing Prodigal results for {genome_id}: {e}")

    # Extract tRNA features from tRNAscan-SE
    if "trnascan" in annotation_files:
        trnascan_files = annotation_files["trnascan"]
        if "error" not in trnascan_files:
            try:
                trna_features = parse_trnascan_results(
                    trnascan_file=Path(trnascan_files.get("trnascan", "")),
                    genome_id=genome_id
                )
                genome_features.extend(trna_features)
                logger.debug(f"Extracted {len(trna_features)} tRNA features from {genome_id}")
            except Exception as e:
                logger.error(f"Error parsing tRNAscan results for {genome_id}: {e}")

    # Extract rRNA features from Barrnap
    if "barrnap" in annotation_files:
        barrnap_files = annotation_files["barrnap"]
        if "error" not in barrnap_files:
            try:
                rrna_features = parse_barrnap_results(
                    gff_file=Path(barrnap_files.get("gff", "")),
                    genome_id=genome_id
                )
                genome_features.extend(rrna_features)
                logger.debug(f"Extracted {len(rrna_features)} rRNA features from {genome_id}")
            except Exception as e:
                logger.error(f"Error parsing Barrnap results for {genome_id}: {e}")

    # Extract intergenic regions if requested
    if extract_intergenic:
        genome_path = genome_data.get("genome_path")
        if genome_path:
            try:
                intergenic_features = extract_intergenic_regions(
                    annotated_features=genome_features,
                    genome_file=Path(genome_path),
                    genome_id=genome_id,
                    min_length=min_intergenic_length
                )
                genome_features.extend(intergenic_features)
                logger.debug(
                    f"Extracted {len(intergenic_features)} intergenic features from {genome_id}"
                )
            except Exception as e:
                logger.error(f"Error extracting intergenic regions for {genome_id}: {e}")
        else:
            logger.warning(f"No genome path found for {genome_id}, skipping intergenic extraction")

    logger.info(f"Extracted {len(genome_features)} total features from {genome_id}")
    return genome_features


def extract_features_from_annotations(
    annotation_results: Dict[str, Dict[str, Any]],
    output_dir: Path,
    create_indices: bool = True,
    extract_intergenic: bool = True,
    min_intergenic_length: int = 50
) -> Tuple[List[GenomeFeature], Dict[str, Any]]:
    """
    Extract standardized features from annotation results.

    Args:
        annotation_results: Results from annotation stage
        output_dir: Directory for extracted features
        create_indices: Whether to create sequence indices
        extract_intergenic: Whether to extract intergenic regions
        min_intergenic_length: Minimum length for intergenic regions (bp)

    Returns:
        Tuple of (feature_list, sequence_indices)

    Raises:
        ValueError: Invalid annotation data
        IOError: File I/O error
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    all_features = []
    sequence_indices = {}

    # Determine number of parallel processes (limit to avoid overwhelming system)
    max_workers = min(len(annotation_results), multiprocessing.cpu_count(), 8)

    # For small datasets, don't use parallel processing overhead
    if len(annotation_results) <= 2:
        logger.info(f"Processing {len(annotation_results)} genomes sequentially")
        for genome_id, genome_data in annotation_results.items():
            logger.info(f"Extracting features from {genome_id}")
            genome_features = _extract_features_from_single_genome(
                genome_id=genome_id,
                genome_data=genome_data,
                extract_intergenic=extract_intergenic,
                min_intergenic_length=min_intergenic_length
            )
            all_features.extend(genome_features)
    else:
        logger.info(
            f"Processing {len(annotation_results)} genomes in parallel using {max_workers} processes"
        )

        # Use ProcessPoolExecutor for parallel genome processing
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all genome extraction jobs
            future_to_genome = {
                executor.submit(
                    _extract_features_from_single_genome,
                    genome_id,
                    genome_data,
                    extract_intergenic,
                    min_intergenic_length
                ): genome_id
                for genome_id, genome_data in annotation_results.items()
            }

            # Collect results as they complete
            for future in as_completed(future_to_genome):
                genome_id = future_to_genome[future]
                try:
                    genome_features = future.result()
                    all_features.extend(genome_features)
                    logger.info(f"Completed feature extraction for {genome_id}")
                except Exception as e:
                    logger.error(f"Failed to extract features from {genome_id}: {e}")

        logger.info(f"Feature extraction complete for {len(annotation_results)} genomes")

    # Create sequence references if requested
    if create_indices:
        sequence_indices = create_sequence_references(
            features=all_features,
            annotation_results=annotation_results,
            output_dir=output_dir
        )

    # Save feature summary
    save_feature_summary(all_features, output_dir / "feature_summary.tsv")

    logger.info(
        f"Extracted {len(all_features)} total features from {len(annotation_results)} genomes"
    )
    return all_features, sequence_indices


def parse_prodigal_results(
    gff_file: Path,
    protein_file: Path,
    nucleotide_file: Path,
    genome_id: str
) -> List[GenomeFeature]:
    """
    Parse Prodigal GFF output into standardized features.

    Args:
        gff_file: Prodigal GFF output
        protein_file: Protein sequences FASTA
        nucleotide_file: Nucleotide sequences FASTA
        genome_id: Identifier for the genome

    Returns:
        List of GenomeFeature objects
    """
    features = []

    if not gff_file.exists():
        logger.warning(f"Prodigal GFF file not found: {gff_file}")
        return features

    try:
        # Parse GFF file for gene coordinates and annotations
        # GFF format: seqname source feature start end score strand frame attributes
        gene_counter = 0

        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue

                # Parse GFF line
                try:
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue

                    seqname, source, feature_type, start_str, end_str, score, strand, frame, attributes = (
                        parts
                    )

                    # Only process CDS features
                    if feature_type != 'CDS':
                        continue

                    gene_counter += 1
                    start = int(start_str)
                    end = int(end_str)
                    contig_id = seqname

                    # Parse attributes field (ID=1_1;partial=00;start_type=ATG;...)
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value

                    # Get gene ID from attributes
                    gene_id = attr_dict.get('ID', f'gene_{gene_counter}')

                    # Extract product information if available
                    product = "hypothetical protein"  # Default for Prodigal

                    # Generate PanGenomePlus ID
                    pgp_id = generate_pangenomeplus_id(
                        genome_id=genome_id,
                        feature_type=FeatureType.PROTEIN_CODING_GENE,
                        original_id=gene_id,
                        contig_id=contig_id,
                        start=start,
                        end=end
                    )

                    genome_feature = GenomeFeature(
                        pangenomeplus_id=pgp_id,
                        original_id=gene_id,
                        genome_id=genome_id,
                        feature_type=FeatureType.PROTEIN_CODING_GENE,
                        contig_id=contig_id,
                        start=start,
                        end=end,
                        strand=strand,
                        product=product,
                        source_tool="prodigal",
                        attributes=attr_dict
                    )

                    features.append(genome_feature)

                except (ValueError, IndexError) as parse_error:
                    logger.warning(f"Could not parse GFF line {line_num} in {gff_file}: {parse_error}")
                    logger.debug(f"Problematic line: {line}")
                    continue

    except Exception as e:
        logger.error(f"Error parsing Prodigal GFF file {gff_file}: {e}")

    logger.info(f"Parsed {len(features)} protein features from {gff_file}")
    return features


def parse_trnascan_results(
    trnascan_file: Path,
    genome_id: str
) -> List[GenomeFeature]:
    """
    Parse tRNAscan-SE output into standardized features.

    Args:
        trnascan_file: tRNAscan-SE output file
        genome_id: Identifier for the genome

    Returns:
        List of GenomeFeature objects
    """
    features = []

    if not trnascan_file.exists():
        logger.warning(f"tRNAscan file not found: {trnascan_file}")
        return features

    try:
        with open(trnascan_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('Sequence'):
                    continue

                parts = line.split()
                if len(parts) < 9:
                    continue

                # Parse tRNAscan output format
                contig_id = parts[0]
                tRNA_num = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                tRNA_type = parts[4]
                anticodon = parts[5]
                score = float(parts[8])

                # Determine strand
                strand = '+' if start < end else '-'
                if strand == '-':
                    start, end = end, start

                # Generate IDs
                original_id = f"tRNA_{tRNA_num}_{tRNA_type}"
                product = f"tRNA-{tRNA_type}({anticodon})"

                pgp_id = generate_pangenomeplus_id(
                    genome_id=genome_id,
                    feature_type=FeatureType.TRNA,
                    original_id=original_id,
                    contig_id=contig_id,
                    start=start,
                    end=end
                )

                genome_feature = GenomeFeature(
                    pangenomeplus_id=pgp_id,
                    original_id=original_id,
                    genome_id=genome_id,
                    feature_type=FeatureType.TRNA,
                    contig_id=contig_id,
                    start=start,
                    end=end,
                    strand=strand,
                    product=product,
                    source_tool="trnascan",
                    attributes={
                        "tRNA_type": tRNA_type,
                        "anticodon": anticodon,
                        "score": score
                    }
                )

                features.append(genome_feature)

    except Exception as e:
        logger.error(f"Error parsing tRNAscan file {trnascan_file}: {e}")

    return features


def parse_barrnap_results(
    gff_file: Path,
    genome_id: str
) -> List[GenomeFeature]:
    """
    Parse Barrnap GFF output into standardized features.

    Args:
        gff_file: Barrnap GFF output file
        genome_id: Identifier for the genome

    Returns:
        List of GenomeFeature objects
    """
    features = []

    if not gff_file.exists():
        logger.warning(f"Barrnap GFF file not found: {gff_file}")
        return features

    try:
        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                # Parse GFF format
                contig_id = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                score = parts[5]
                strand = parts[6]
                attributes_str = parts[8]

                # Parse attributes
                attributes = {}
                for attr in attributes_str.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value

                # Extract rRNA information
                gene_name = attributes.get('Name', f'rRNA_{line_num}')
                product = attributes.get('product', gene_name)

                pgp_id = generate_pangenomeplus_id(
                    genome_id=genome_id,
                    feature_type=FeatureType.RRNA,
                    original_id=gene_name,
                    contig_id=contig_id,
                    start=start,
                    end=end
                )

                genome_feature = GenomeFeature(
                    pangenomeplus_id=pgp_id,
                    original_id=gene_name,
                    genome_id=genome_id,
                    feature_type=FeatureType.RRNA,
                    contig_id=contig_id,
                    start=start,
                    end=end,
                    strand=strand,
                    product=product,
                    source_tool="barrnap",
                    attributes=attributes
                )

                features.append(genome_feature)

    except Exception as e:
        logger.error(f"Error parsing Barrnap GFF file {gff_file}: {e}")

    return features


def extract_intergenic_regions(
    annotated_features: List[GenomeFeature],
    genome_file: Path,
    genome_id: str,
    min_length: int = 50
) -> List[GenomeFeature]:
    """
    Extract intergenic regions (gaps between annotated features) from a genome.

    Args:
        annotated_features: List of annotated features for this genome
        genome_file: Path to the original genome FASTA file
        genome_id: Identifier for the genome
        min_length: Minimum length for intergenic regions (default: 50bp)

    Returns:
        List of GenomeFeature objects representing intergenic regions
    """
    if not genome_file.exists():
        logger.warning(f"Genome file not found: {genome_file}")
        return []

    intergenic_features = []

    try:
        # Load genome sequences
        genome_sequences = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_sequences[record.id] = record

        # Group features by contig and sort by position
        features_by_contig = {}
        for feature in annotated_features:
            contig_id = feature.contig_id
            if contig_id not in features_by_contig:
                features_by_contig[contig_id] = []
            features_by_contig[contig_id].append(feature)

        # Sort features by start position within each contig
        for contig_id in features_by_contig:
            features_by_contig[contig_id].sort(key=lambda f: f.start)

        # Process each contig to find intergenic regions
        for contig_id, contig_features in features_by_contig.items():
            if contig_id not in genome_sequences:
                logger.warning(f"Contig {contig_id} not found in genome file {genome_file}")
                continue

            contig_sequence = genome_sequences[contig_id]
            contig_length = len(contig_sequence.seq)

            # Find gaps between features
            gaps = []

            # Gap before first feature
            if contig_features and contig_features[0].start > 1:
                gaps.append((1, contig_features[0].start - 1))

            # Gaps between features
            for i in range(len(contig_features) - 1):
                current_feature = contig_features[i]
                next_feature = contig_features[i + 1]

                gap_start = current_feature.end + 1
                gap_end = next_feature.start - 1

                if gap_end >= gap_start:  # Valid gap
                    gaps.append((gap_start, gap_end))

            # Gap after last feature
            if contig_features and contig_features[-1].end < contig_length:
                gaps.append((contig_features[-1].end + 1, contig_length))

            # Create intergenic features for gaps >= min_length
            for gap_num, (gap_start, gap_end) in enumerate(gaps, 1):
                gap_length = gap_end - gap_start + 1

                if gap_length >= min_length:
                    # Extract sequence for this intergenic region
                    gap_sequence = str(contig_sequence.seq[gap_start-1:gap_end])  # Convert to 0-based

                    # Generate intergenic region ID
                    original_id = f"intergenic_{gap_num}"
                    pgp_id = generate_pangenomeplus_id(
                        genome_id=genome_id,
                        feature_type=FeatureType.INTERGENIC_REGION,
                        original_id=original_id,
                        contig_id=contig_id,
                        start=gap_start,
                        end=gap_end
                    )

                    # Create GenomeFeature for intergenic region
                    intergenic_feature = GenomeFeature(
                        pangenomeplus_id=pgp_id,
                        original_id=original_id,
                        genome_id=genome_id,
                        feature_type=FeatureType.INTERGENIC_REGION,
                        contig_id=contig_id,
                        start=gap_start,
                        end=gap_end,
                        strand='.',  # Intergenic regions don't have strand
                        product=f"intergenic region ({gap_length} bp)",
                        source_tool="intergenic_extraction",
                        attributes={
                            "length": gap_length,
                            "sequence": gap_sequence
                        }
                    )

                    intergenic_features.append(intergenic_feature)

        logger.info(f"Extracted {len(intergenic_features)} intergenic regions from {genome_id}")

    except Exception as e:
        logger.error(f"Error extracting intergenic regions from {genome_file}: {e}")

    return intergenic_features


def create_sequence_references(
    features: List[GenomeFeature],
    annotation_results: Dict[str, Dict[str, Any]],
    output_dir: Path
) -> Dict[str, Any]:
    """
    Create reference-only sequence storage system.

    Args:
        features: List of genome features
        annotation_results: Raw annotation data
        output_dir: Directory for sequence files

    Returns:
        Dictionary with sequence index information
    """
    sequence_dir = output_dir / "sequences"
    sequence_dir.mkdir(exist_ok=True)

    indices = {
        "protein_sequences": {},
        "nucleotide_sequences": {},
        "feature_locations": {}
    }

    # Group features by genome and type
    by_genome = {}
    for feature in features:
        if feature.genome_id not in by_genome:
            by_genome[feature.genome_id] = {}
        if feature.feature_type not in by_genome[feature.genome_id]:
            by_genome[feature.genome_id][feature.feature_type] = []
        by_genome[feature.genome_id][feature.feature_type].append(feature)

    # Create sequence files for each feature type
    for genome_id, feature_types in by_genome.items():
        genome_annotation = annotation_results.get(genome_id, {})

        for feature_type, genome_features in feature_types.items():
            if feature_type == FeatureType.PROTEIN_CODING_GENE:
                # Reference protein sequences from Prodigal output
                prodigal_files = genome_annotation.get("annotation_files", {}).get("prodigal", {})
                protein_file = prodigal_files.get("proteins")
                if protein_file and Path(protein_file).exists():
                    indices["protein_sequences"][genome_id] = protein_file

            # Store feature locations
            for feature in genome_features:
                indices["feature_locations"][feature.pangenomeplus_id] = {
                    "genome_id": feature.genome_id,
                    "contig_id": feature.contig_id,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand,
                    "feature_type": feature.feature_type.value
                }

    return indices


def save_feature_summary(features: List[GenomeFeature], output_file: Path) -> None:
    """
    Save feature summary to TSV file.

    Args:
        features: List of genome features
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        # Header
        f.write("pangenomeplus_id\toriginal_id\tgenome_id\tfeature_type\tcontig_id\t"
                "start\tend\tstrand\tproduct\tsource_tool\n")

        # Features
        for feature in features:
            f.write(f"{feature.pangenomeplus_id}\t{feature.original_id}\t{feature.genome_id}\t"
                   f"{feature.feature_type.value}\t{feature.contig_id}\t{feature.start}\t"
                   f"{feature.end}\t{feature.strand}\t{feature.product}\t{feature.source_tool}\n")

    logger.info(f"Saved feature summary to {output_file}")