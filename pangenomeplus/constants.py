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

# Pangenome openness classification thresholds
# Used in: pangenome_analysis.py for classifying pangenome as open/closed/intermediate
# Based on growth rate (families gained per genome added)
PANGENOME_OPEN_THRESHOLD = 50  # Growth rate > 50 families/genome = open pangenome
PANGENOME_CLOSED_THRESHOLD = 10  # Growth rate < 10 families/genome = closed pangenome
# Intermediate pangenomes fall between these values (10-50 families/genome)

# Pipeline configuration defaults
MIN_INTERGENIC_LENGTH = 50  # Minimum intergenic region length in base pairs
CHECKPOINT_INTERVAL = 10  # Save checkpoint every N genomes during processing


class ClusteringDefaults:
    """Default clustering parameters for each feature type.

    Used in: clustering.py lines 73-77, cli.py for defaults
    """
    PROTEIN = {"identity": 0.8, "coverage": 0.8}
    INTERGENIC = {"identity": 0.7, "coverage": 0.5}  # Relaxed for high variability
    TRNA = {"identity": 0.9, "coverage": 0.9}  # Strict for conserved features
    RRNA = {"identity": 0.95, "coverage": 0.95}  # Very strict for highly conserved
    CRISPR = {"identity": 0.8, "coverage": 0.6}  # Moderate for spacer diversity

    # MMseqs2 general parameters
    SENSITIVITY = 4.0  # Default sensitivity level
    MAX_SEQS = 10000  # Maximum sequences per query


class VizParams:
    """Visualization parameters for plots and figures.

    Used in: pangenome_analysis.py for all visualizations
    """
    # Figure sizes
    DEFAULT_FIGSIZE = (10, 6)
    HEATMAP_FIGSIZE = (14, 8)
    PIE_FIGSIZE = (8, 6)
    DPI = 300

    # Data display limits
    TOP_VARIABLE_FAMILIES = 50  # Number of families to show in heatmap
    FAMILY_NAME_TRUNCATE = 15  # Truncate family names to this length

    # Font sizes
    FONT_SIZE_LABEL = 12  # Axis labels
    FONT_SIZE_TITLE = 14  # Plot titles
    FONT_SIZE_LEGEND = 10  # Legend text

    # Visual styling
    ALPHA_TRANSPARENCY = 0.2  # Alpha for uncertainty bands
    GRID_ALPHA = 0.3  # Alpha for grid lines
    PIE_START_ANGLE = 90  # Starting angle for pie chart
    PIE_EXPLODE_CORE = 0.05  # Explosion distance for core slice
    LABEL_ROTATION = 90  # Rotation angle for x-axis labels

    # Heatmap sizing
    HEATMAP_MIN_HEIGHT = 8  # Minimum heatmap height
    HEATMAP_HEIGHT_PER_GENOME = 0.4  # Additional height per genome


# File handling
# Used in: pipeline.py line 521, cli.py
GENOME_EXTENSIONS = [".fasta", ".fna", ".fa"]

# Output file names
# Used in: outputs.py, pipeline.py
OUTPUT_FILES = {
    "presence_absence": "presence_absence_matrix.csv",
    "family_summary": "family_summary.tsv",
    "transformer": "pangenome_transformer.txt",
    "mappings": "compact_id_mappings.json",
    "markdown_report": "pangenome_report.md"
}


class ToolDefaults:
    """Default parameters for external bioinformatics tools.

    Used in: extraction.py, cli.py for argument defaults
    """
    # Prodigal defaults
    PRODIGAL_MODE = "single"  # Default mode for gene prediction
    TRANSLATION_TABLE = 11  # Bacterial/Archaeal genetic code

    # tRNAscan-SE defaults
    TRNASCAN_SCORE_CUTOFF = 20.0  # Score threshold for tRNA detection
    TRNASCAN_MODEL = "bacteria"  # Default organism model

    # Barrnap defaults
    BARRNAP_EVALUE = 1e-6  # E-value threshold for rRNA detection
    BARRNAP_KINGDOM = "bac"  # Default kingdom (bacteria)

    # MINCED defaults
    MINCED_MIN_REPEATS = 3  # Minimum repeats for CRISPR array
    MINCED_REPEAT_LENGTH_MIN = 23
    MINCED_REPEAT_LENGTH_MAX = 47
    MINCED_SPACER_LENGTH_MIN = 26
    MINCED_SPACER_LENGTH_MAX = 50
