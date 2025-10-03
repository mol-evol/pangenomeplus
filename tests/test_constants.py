"""Tests for constants module validation."""

import pytest

from pangenomeplus.constants import (
    CHECKPOINT_INTERVAL,
    CLOUD_GENOME_THRESHOLD,
    CORE_GENOME_THRESHOLD,
    GENOME_EXTENSIONS,
    MIN_INTERGENIC_LENGTH,
    OUTPUT_FILES,
    ClusteringDefaults,
    FeatureType,
    ToolDefaults,
    VizParams,
)


class TestThresholdConstants:
    """Test pangenome classification thresholds."""

    def test_threshold_ranges(self):
        """Test thresholds are in valid ranges."""
        assert 0.0 < CORE_GENOME_THRESHOLD <= 1.0
        assert 0.0 <= CLOUD_GENOME_THRESHOLD < 1.0
        assert CLOUD_GENOME_THRESHOLD < CORE_GENOME_THRESHOLD

    def test_threshold_biological_validity(self):
        """Test thresholds match biological expectations."""
        # Core should be high (most genomes)
        assert CORE_GENOME_THRESHOLD >= 0.9

        # Cloud should be low (few genomes)
        assert CLOUD_GENOME_THRESHOLD <= 0.2

        # Accessory is implicitly between them
        accessory_min = CLOUD_GENOME_THRESHOLD
        accessory_max = CORE_GENOME_THRESHOLD
        assert accessory_max - accessory_min > 0.5  # Reasonable range

    def test_intergenic_length(self):
        """Test minimum intergenic length is reasonable."""
        assert MIN_INTERGENIC_LENGTH >= 30  # Biological minimum
        assert MIN_INTERGENIC_LENGTH <= 200  # Not too restrictive


class TestFeatureTypeConstants:
    """Test feature type definitions."""

    def test_feature_type_codes(self):
        """Test single-character feature codes."""
        assert FeatureType.PROTEIN == "P"
        assert FeatureType.INTERGENIC == "I"
        assert FeatureType.TRNA == "T"
        assert FeatureType.RRNA == "R"
        assert FeatureType.CRISPR == "C"

    def test_feature_type_uniqueness(self):
        """Test all feature type codes are unique."""
        codes = [
            FeatureType.PROTEIN,
            FeatureType.INTERGENIC,
            FeatureType.TRNA,
            FeatureType.RRNA,
            FeatureType.CRISPR,
        ]
        assert len(codes) == len(set(codes))

    def test_feature_type_names(self):
        """Test human-readable names exist for all types."""
        assert len(FeatureType.NAMES) == 5
        assert FeatureType.NAMES[FeatureType.PROTEIN] == "Proteins"
        assert FeatureType.NAMES[FeatureType.INTERGENIC] == "Intergenic"
        assert FeatureType.NAMES[FeatureType.TRNA] == "tRNAs"
        assert FeatureType.NAMES[FeatureType.RRNA] == "rRNAs"
        assert FeatureType.NAMES[FeatureType.CRISPR] == "CRISPR"

    def test_all_types_set(self):
        """Test ALL_TYPES contains all feature types."""
        assert len(FeatureType.ALL_TYPES) == 5
        assert FeatureType.PROTEIN in FeatureType.ALL_TYPES
        assert FeatureType.INTERGENIC in FeatureType.ALL_TYPES


class TestClusteringDefaults:
    """Test clustering parameter defaults."""

    def test_clustering_identity_ranges(self):
        """Test identity thresholds are valid."""
        for params in [
            ClusteringDefaults.PROTEIN,
            ClusteringDefaults.INTERGENIC,
            ClusteringDefaults.TRNA,
            ClusteringDefaults.RRNA,
            ClusteringDefaults.CRISPR,
        ]:
            assert 0.0 < params["identity"] <= 1.0
            assert 0.0 < params["coverage"] <= 1.0

    def test_clustering_biological_logic(self):
        """Test clustering parameters match biological expectations."""
        # rRNAs most conserved (highest threshold)
        assert ClusteringDefaults.RRNA["identity"] >= 0.95

        # tRNAs highly conserved
        assert ClusteringDefaults.TRNA["identity"] >= 0.9

        # Proteins moderate conservation
        assert 0.7 <= ClusteringDefaults.PROTEIN["identity"] <= 0.9

        # Intergenic most variable (lowest threshold)
        assert ClusteringDefaults.INTERGENIC["identity"] <= 0.8
        assert ClusteringDefaults.INTERGENIC["coverage"] <= 0.6

    def test_mmseqs_parameters(self):
        """Test MMseqs2 general parameters."""
        assert 1.0 <= ClusteringDefaults.SENSITIVITY <= 8.0  # MMseqs2 range
        assert ClusteringDefaults.MAX_SEQS >= 1000  # Reasonable minimum


class TestToolDefaults:
    """Test external tool default parameters."""

    def test_prodigal_defaults(self):
        """Test Prodigal parameters."""
        assert ToolDefaults.PRODIGAL_MODE in ["single", "meta"]
        assert ToolDefaults.TRANSLATION_TABLE in range(1, 26)  # Valid genetic codes

    def test_trnascan_defaults(self):
        """Test tRNAscan-SE parameters."""
        assert ToolDefaults.TRNASCAN_SCORE_CUTOFF > 0
        assert ToolDefaults.TRNASCAN_MODEL in ["bacteria", "archaea", "eukaryota"]

    def test_barrnap_defaults(self):
        """Test Barrnap parameters."""
        assert ToolDefaults.BARRNAP_EVALUE > 0
        assert ToolDefaults.BARRNAP_EVALUE < 1e-3  # Should be stringent
        assert ToolDefaults.BARRNAP_KINGDOM in ["bac", "arc", "euk", "mito"]

    def test_minced_defaults(self):
        """Test MINCED CRISPR detection parameters."""
        assert ToolDefaults.MINCED_MIN_REPEATS >= 2
        assert ToolDefaults.MINCED_REPEAT_LENGTH_MIN < ToolDefaults.MINCED_REPEAT_LENGTH_MAX
        assert ToolDefaults.MINCED_SPACER_LENGTH_MIN < ToolDefaults.MINCED_SPACER_LENGTH_MAX

        # Biological constraints for CRISPR
        assert 20 <= ToolDefaults.MINCED_REPEAT_LENGTH_MIN <= 30
        assert 40 <= ToolDefaults.MINCED_REPEAT_LENGTH_MAX <= 50
        assert 20 <= ToolDefaults.MINCED_SPACER_LENGTH_MIN <= 30
        assert 40 <= ToolDefaults.MINCED_SPACER_LENGTH_MAX <= 60


class TestVisualizationParams:
    """Test visualization parameters."""

    def test_figure_sizes(self):
        """Test figure size tuples."""
        assert len(VizParams.DEFAULT_FIGSIZE) == 2
        assert len(VizParams.HEATMAP_FIGSIZE) == 2
        assert len(VizParams.PIE_FIGSIZE) == 2

        # Reasonable sizes
        for size in [VizParams.DEFAULT_FIGSIZE, VizParams.HEATMAP_FIGSIZE, VizParams.PIE_FIGSIZE]:
            assert 5 <= size[0] <= 20  # Width
            assert 4 <= size[1] <= 15  # Height

    def test_dpi_setting(self):
        """Test DPI for publication quality."""
        assert VizParams.DPI >= 150  # Minimum for decent quality
        assert VizParams.DPI <= 600  # Maximum reasonable value

    def test_top_variable_families(self):
        """Test number of families to show."""
        assert 10 <= VizParams.TOP_VARIABLE_FAMILIES <= 100


class TestOutputFiles:
    """Test output file name constants."""

    def test_output_file_names(self):
        """Test all output files have proper extensions."""
        assert OUTPUT_FILES["presence_absence"].endswith(".csv")
        assert OUTPUT_FILES["family_summary"].endswith(".tsv")
        assert OUTPUT_FILES["transformer"].endswith(".txt")
        assert OUTPUT_FILES["mappings"].endswith(".json")
        assert OUTPUT_FILES["markdown_report"].endswith(".md")

    def test_output_file_uniqueness(self):
        """Test all output files have unique names."""
        filenames = list(OUTPUT_FILES.values())
        assert len(filenames) == len(set(filenames))


class TestGenomeExtensions:
    """Test genome file extensions."""

    def test_standard_extensions(self):
        """Test standard genome file extensions are included."""
        assert ".fasta" in GENOME_EXTENSIONS
        assert ".fna" in GENOME_EXTENSIONS
        assert ".fa" in GENOME_EXTENSIONS

    def test_extension_format(self):
        """Test extensions start with dot."""
        assert all(ext.startswith(".") for ext in GENOME_EXTENSIONS)


class TestCheckpointSettings:
    """Test pipeline checkpoint settings."""

    def test_checkpoint_interval(self):
        """Test checkpoint interval is reasonable."""
        assert CHECKPOINT_INTERVAL >= 1  # At least every genome
        assert CHECKPOINT_INTERVAL <= 100  # Not too infrequent