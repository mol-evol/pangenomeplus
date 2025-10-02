"""Global constants for PanGenomePlus.

This module centralizes all magic numbers, thresholds, and repeated
dictionaries to ensure DRY compliance and single source of truth.
"""


class FeatureType:
    """Feature type identifiers and mappings."""

    # Single-character feature type codes
    PROTEIN = "P"
    INTERGENIC = "I"
    TRNA = "T"
    RRNA = "R"
    CRISPR = "C"

    # Set of all valid feature types
    ALL_TYPES = {PROTEIN, INTERGENIC, TRNA, RRNA, CRISPR}

    # Human-readable feature type names
    # Used across: pipeline.py, outputs.py, pangenome_analysis.py, cli.py
    NAMES = {
        PROTEIN: "Proteins",
        INTERGENIC: "Intergenic",
        TRNA: "tRNAs",
        RRNA: "rRNAs",
        CRISPR: "CRISPR",
    }


# Pangenome classification thresholds
# Used in: clustering.py, outputs.py
CORE_GENOME_THRESHOLD = 0.95  # Present in ≥95% of genomes
CLOUD_GENOME_THRESHOLD = 0.15  # Present in <15% of genomes
# Accessory families are implicitly between these thresholds (15-95%)

# Pipeline configuration defaults
MIN_INTERGENIC_LENGTH = 50  # Minimum intergenic region length in base pairs
CHECKPOINT_INTERVAL = 10  # Save checkpoint every N genomes during processing
