# PanGenomePlus API Reference

Complete technical documentation for developers using PanGenomePlus programmatically or contributing to the codebase.

## Table of Contents

1. [Core Architecture](#core-architecture)
2. [Main Pipeline API](#main-pipeline-api)
3. [Core Types and Data Structures](#core-types-and-data-structures)
4. [Module APIs](#module-apis)
5. [Configuration System](#configuration-system)
6. [Exception Handling](#exception-handling)
7. [Utility Functions](#utility-functions)
8. [Extension Points](#extension-points)
9. [Testing Framework](#testing-framework)
10. [Development Guidelines](#development-guidelines)

## Core Architecture

PanGenomePlus follows a modular pipeline architecture:

```
Input Genomes -> Annotation -> Feature Extraction -> Clustering -> Family Assignment -> Output Generation
```

### Package Structure
```
pangenomeplus/
|-- __main__.py              # Command-line interface
|-- pipeline.py              # Main pipeline orchestration
|-- core/
|   |-- types.py            # Data structures and enums
|   `-- exceptions.py       # Custom exception classes
|-- modules/
|   |-- annotation.py       # Genome annotation
|   |-- extraction.py       # Feature extraction
|   |-- clustering.py       # MMseqs2 clustering
|   |-- families.py         # Gene family assignment
|   `-- output.py           # Output format generation
`-- utils/
    `-- config.py           # Configuration management
```

## Main Pipeline API

### run_pangenome_pipeline()

The primary entry point for programmatic use.

```python
from pangenomeplus.pipeline import run_pangenome_pipeline
from pathlib import Path

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
    Execute complete pangenome analysis pipeline.

    Args:
        genomes_dir: Directory containing input genome FASTA files
        output_dir: Directory for pipeline outputs
        config: Complete pipeline configuration dictionary
        resume: Whether to resume from existing checkpoints
        stop_after: Optional stage to stop after execution
        validate_inputs: Whether to validate inputs before processing
        log_level: Logging verbosity level

    Returns:
        Dictionary containing pipeline results and metadata:
        {
            "start_time": float,
            "end_time": float,
            "runtime_seconds": float,
            "runtime_formatted": str,
            "pipeline_version": str,
            "config": Dict[str, Any],
            "n_genomes": int,
            "n_features": int,
            "n_families": int,
            "annotation_results": Dict[str, Any],
            "cluster_results": Dict[str, Any],
            "output_files": Dict[str, str]
        }

    Raises:
        ValidationError: If input validation fails
        PipelineError: If pipeline execution fails
        ToolExecutionError: If external tool execution fails
    """
```

### Example Usage

```python
from pangenomeplus.pipeline import run_pangenome_pipeline
from pangenomeplus.utils.config import create_default_configuration

# Create configuration
config = create_default_configuration()
config["clustering"]["protein"]["identity"] = 0.9
config["output"]["formats"] = ["transformer", "presence_absence"]

# Run pipeline
results = run_pangenome_pipeline(
    genomes_dir="./genomes",
    output_dir="./results",
    config=config,
    log_level="INFO"
)

print(f"Processed {results['n_genomes']} genomes")
print(f"Created {results['n_families']} gene families")
print(f"Runtime: {results['runtime_formatted']}")
```

### analyze_dataset_and_choose_strategy()

```python
def analyze_dataset_and_choose_strategy(
    genomes_dir: Union[str, Path],
    sample_size: int = 10
) -> Dict[str, Any]:
    """
    Analyze dataset characteristics for informational purposes.

    Args:
        genomes_dir: Directory containing genome files
        sample_size: Number of genomes to sample for analysis

    Returns:
        Dictionary with dataset information:
        {
            "n_genomes": int,
            "memory_available_gb": int,
            "description": str
        }
    """
```

### validate_pipeline_inputs()

```python
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
```

## Core Types and Data Structures

### FeatureType Enum

```python
from pangenomeplus.core.types import FeatureType

class FeatureType(Enum):
    """Genomic feature types supported by PanGenomePlus."""
    PROTEIN_CODING_GENE = "protein_coding_gene"
    TRNA = "trna"
    RRNA = "rrna"
    CRISPR = "crispr"
    INTERGENIC = "intergenic"
```

### GenomeFeature Dataclass

```python
@dataclass
class GenomeFeature:
    """
    Represents a genomic feature extracted from annotations.

    Attributes:
        pangenomeplus_id: Unique identifier within PanGenomePlus
        original_id: Original identifier from annotation tool
        genome_id: Identifier of the source genome
        feature_type: Type of genomic feature
        start: Start coordinate (1-based)
        end: End coordinate (1-based, inclusive)
        strand: Strand orientation ('+', '-', or '.')
        sequence: Nucleotide or amino acid sequence
        annotation: Additional annotation information
    """
    pangenomeplus_id: str
    original_id: str
    genome_id: str
    feature_type: FeatureType
    start: int
    end: int
    strand: str
    sequence: str
    annotation: Dict[str, Any] = field(default_factory=dict)
```

### ValidationResult

```python
@dataclass
class ValidationResult:
    """
    Result of input validation operations.

    Attributes:
        is_valid: Whether validation passed
        errors: List of error messages
        warnings: List of warning messages
        details: Additional validation details
    """
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)
```

### validate_genome_feature()

```python
def validate_genome_feature(feature: GenomeFeature) -> None:
    """
    Validate genome feature data integrity.

    Args:
        feature: GenomeFeature instance to validate

    Raises:
        ValueError: If feature data is invalid
    """
```

## Module APIs

### Annotation Module (pangenomeplus.modules.annotation)

#### annotate_genomes()

```python
def annotate_genomes(
    genome_paths: List[Path],
    output_dir: Path,
    tools: List[str] = None,
    threads: int = 4,
    config: Dict[str, Any] = None
) -> Dict[str, Any]:
    """
    Annotate multiple genomes using specified tools.

    Args:
        genome_paths: List of paths to genome FASTA files
        output_dir: Directory for annotation outputs
        tools: List of annotation tools to use
               ["prodigal", "trnascan", "barrnap", "minced"]
        threads: Number of parallel threads
        config: Tool-specific configuration options

    Returns:
        Dictionary mapping genome_id to annotation file paths:
        {
            "genome_id": {
                "prodigal": {"proteins": Path, "genes": Path},
                "trnascan": {"trnas": Path},
                "barrnap": {"rrnas": Path},
                "minced": {"crisprs": Path}
            }
        }
    """
```

#### Tool-specific functions

```python
def run_prodigal(
    genome_path: Path,
    output_dir: Path,
    genetic_code: int = 11,
    threads: int = 1
) -> Dict[str, Path]:
    """Run Prodigal gene prediction."""

def run_trnascan(
    genome_path: Path,
    output_dir: Path,
    threads: int = 1
) -> Dict[str, Path]:
    """Run tRNAscan-SE for tRNA detection."""

def run_barrnap(
    genome_path: Path,
    output_dir: Path,
    kingdom: str = "bac"
) -> Dict[str, Path]:
    """Run Barrnap for rRNA detection."""

def run_minced(
    genome_path: Path,
    output_dir: Path
) -> Dict[str, Path]:
    """Run MINCED for CRISPR detection."""
```

### Feature Extraction Module (pangenomeplus.modules.extraction)

#### extract_features_from_annotations()

```python
def extract_features_from_annotations(
    annotation_results: Dict[str, Any],
    output_dir: Path,
    extract_intergenic: bool = True,
    min_intergenic_length: int = 50
) -> Tuple[List[GenomeFeature], Dict[str, int]]:
    """
    Extract genomic features from annotation results.

    Args:
        annotation_results: Output from annotate_genomes()
        output_dir: Directory for feature extraction outputs
        extract_intergenic: Whether to extract intergenic regions
        min_intergenic_length: Minimum length for intergenic regions

    Returns:
        Tuple of (features_list, sequence_indices):
        - features_list: List of GenomeFeature objects
        - sequence_indices: Mapping of feature IDs to sequence indices
    """
```

#### Feature-specific extractors

```python
def extract_protein_features(
    annotation_results: Dict[str, Any]
) -> List[GenomeFeature]:
    """Extract protein-coding gene features."""

def extract_trna_features(
    annotation_results: Dict[str, Any]
) -> List[GenomeFeature]:
    """Extract tRNA gene features."""

def extract_rrna_features(
    annotation_results: Dict[str, Any]
) -> List[GenomeFeature]:
    """Extract rRNA gene features."""

def extract_crispr_features(
    annotation_results: Dict[str, Any]
) -> List[GenomeFeature]:
    """Extract CRISPR array features."""

def extract_intergenic_features(
    annotation_results: Dict[str, Any],
    min_length: int = 50
) -> List[GenomeFeature]:
    """Extract intergenic region features."""
```

### Clustering Module (pangenomeplus.modules.clustering)

#### run_clustering()

```python
def run_clustering(
    features_by_type: Dict[FeatureType, List[GenomeFeature]],
    output_dir: Path,
    config: Dict[str, Any]
) -> Dict[FeatureType, Dict[str, Any]]:
    """
    Cluster features by type using MMseqs2 easy-cluster.

    Args:
        features_by_type: Features grouped by genomic type
        output_dir: Directory for clustering outputs
        config: Clustering configuration parameters

    Returns:
        Dictionary mapping feature types to cluster results:
        {
            FeatureType.PROTEIN_CODING_GENE: {
                "cluster_file": Path,
                "representative_file": Path,
                "cluster_assignments": Dict[str, str]
            }
        }
    """
```

#### MMseqs2 wrapper functions

```python
def run_mmseqs_easy_cluster(
    input_fasta: Path,
    output_prefix: Path,
    identity: float = 0.8,
    coverage: float = 0.8,
    sensitivity: float = 7.5,
    threads: int = 4,
    use_gpu: bool = False,
    memory_limit: Optional[str] = None
) -> Dict[str, Path]:
    """Run MMseqs2 easy-cluster with specified parameters."""

def parse_mmseqs_clusters(
    cluster_file: Path
) -> Dict[str, str]:
    """Parse MMseqs2 cluster assignments."""
```

### Families Module (pangenomeplus.modules.families)

#### assign_gene_families()

```python
def assign_gene_families(
    cluster_results: Dict[FeatureType, Dict[str, Any]],
    id_mappings: Dict[str, str],
    output_dir: Path
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Assign gene family identifiers to clustered features.

    Args:
        cluster_results: Output from run_clustering()
        id_mappings: Mapping from PanGenomePlus IDs to original IDs
        output_dir: Directory for family assignment outputs

    Returns:
        Tuple of (gene_families, family_members):
        - gene_families: Mapping from feature ID to family ID
        - family_members: Mapping from family ID to member list
    """
```

#### Family processing functions

```python
def create_family_identifiers(
    cluster_assignments: Dict[str, str],
    feature_type: FeatureType
) -> Dict[str, str]:
    """Create standardized family identifiers."""

def save_family_assignments(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[str]],
    output_dir: Path
) -> Dict[str, Path]:
    """Save family assignments to files."""
```

### Output Module (pangenomeplus.modules.output)

#### generate_outputs()

```python
def generate_outputs(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[str]],
    features: List[GenomeFeature],
    id_mappings: Dict[str, Dict[str, str]],
    output_dir: Path,
    formats: List[str] = None,
    config: Dict[str, Any] = None,
    pipeline_results: Dict[str, Any] = None,
    annotation_results: Dict[str, Any] = None,
    cluster_results: Dict[str, Any] = None
) -> Dict[str, Path]:
    """
    Generate output files in specified formats.

    Args:
        gene_families: Feature to family mappings
        family_members: Family to members mappings
        features: List of all genomic features
        id_mappings: ID mapping dictionaries
        output_dir: Directory for output files
        formats: List of output formats to generate
        config: Output configuration options
        pipeline_results: Pipeline execution metadata
        annotation_results: Annotation stage results
        cluster_results: Clustering stage results

    Returns:
        Dictionary mapping format names to output file paths
    """
```

#### Format-specific generators

```python
def generate_transformer_format(
    gene_families: Dict[str, str],
    features: List[GenomeFeature],
    output_dir: Path,
    config: Dict[str, Any] = None
) -> Path:
    """Generate transformer format output for ML applications."""

def generate_presence_absence_matrix(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[str]],
    output_dir: Path
) -> Path:
    """Generate binary presence/absence matrix."""

def generate_roary_format(
    gene_families: Dict[str, str],
    family_members: Dict[str, List[str]],
    features: List[GenomeFeature],
    output_dir: Path
) -> Dict[str, Path]:
    """Generate Roary-compatible output files."""

def generate_fasta_outputs(
    family_members: Dict[str, List[str]],
    features: List[GenomeFeature],
    output_dir: Path
) -> Dict[str, Path]:
    """Generate FASTA files for gene families."""
```

## Configuration System

### create_default_configuration()

```python
from pangenomeplus.utils.config import create_default_configuration

def create_default_configuration() -> Dict[str, Any]:
    """
    Create default pipeline configuration.

    Returns:
        Complete default configuration dictionary with all sections:
        {
            "clustering": {
                "protein": {"identity": 0.8, "coverage": 0.8, "sensitivity": 7.5},
                "trna": {"identity": 0.95, "coverage": 0.9, "sensitivity": 7.5},
                "rrna": {"identity": 0.97, "coverage": 0.9, "sensitivity": 7.5},
                "crispr": {"identity": 0.9, "coverage": 0.8, "sensitivity": 6.0},
                "intergenic": {"identity": 0.7, "coverage": 0.6, "sensitivity": 6.0},
                "gpu": {"use_gpu": False, "memory_fraction": 0.8}
            },
            "annotation": {
                "tools": ["prodigal", "trnascan", "barrnap", "minced"],
                "prodigal": {"genetic_code": 11},
                "trnascan": {"kingdom": "bacteria"},
                "barrnap": {"kingdom": "bac"},
                "intergenic": {"extract": True, "min_length": 50}
            },
            "classification": {
                "hard_core_threshold": 1.0,
                "soft_core_threshold": 0.95,
                "shell_threshold": 0.15,
                "cloud_threshold": 0.15
            },
            "output": {
                "formats": ["transformer", "presence_absence"],
                "transformer": {"include_singletons": True, "feature_type_prefix": True},
                "compress": False
            },
            "resources": {
                "threads": 4,
                "memory_limit": None
            },
            "validation": {
                "input_files": True,
                "system_requirements": True
            }
        }
    """
```

### validate_configuration_schema()

```python
def validate_configuration_schema(config: Dict[str, Any]) -> ValidationResult:
    """
    Validate configuration against expected schema.

    Args:
        config: Configuration dictionary to validate

    Returns:
        ValidationResult indicating if configuration is valid
    """
```

### Configuration modification helpers

```python
def update_clustering_parameters(
    config: Dict[str, Any],
    identity: float = None,
    coverage: float = None,
    sensitivity: float = None,
    feature_type: str = "protein"
) -> Dict[str, Any]:
    """Update clustering parameters for specific feature type."""

def update_output_formats(
    config: Dict[str, Any],
    formats: List[str]
) -> Dict[str, Any]:
    """Update output format list in configuration."""

def update_resource_limits(
    config: Dict[str, Any],
    threads: int = None,
    memory_limit: str = None
) -> Dict[str, Any]:
    """Update resource allocation in configuration."""
```

## Exception Handling

### Custom Exceptions

```python
from pangenomeplus.core.exceptions import (
    PipelineError,
    ValidationError,
    ToolExecutionError,
    ConfigurationError
)

class PipelineError(Exception):
    """Base exception for pipeline-related errors."""
    pass

class ValidationError(PipelineError):
    """Raised when input validation fails."""
    pass

class ToolExecutionError(PipelineError):
    """Raised when external tool execution fails."""

    def __init__(self, message: str, tool: str = None, exit_code: int = None):
        super().__init__(message)
        self.tool = tool
        self.exit_code = exit_code

class ConfigurationError(PipelineError):
    """Raised when configuration is invalid."""
    pass
```

### Error handling patterns

```python
from pangenomeplus.core.exceptions import ToolExecutionError

try:
    result = run_pangenome_pipeline(genomes_dir, output_dir, config)
except ValidationError as e:
    print(f"Input validation failed: {e}")
except ToolExecutionError as e:
    print(f"Tool '{e.tool}' failed with exit code {e.exit_code}: {e}")
except PipelineError as e:
    print(f"Pipeline error: {e}")
```

## Utility Functions

### File handling utilities

```python
def list_genome_files(genomes_dir: Path) -> List[Path]:
    """List all genome FASTA files in directory."""

def setup_logging(level: str, log_file: Optional[Path] = None) -> None:
    """Setup logging configuration."""

def create_output_directories(base_dir: Path) -> Dict[str, Path]:
    """Create standard output directory structure."""
```

### Sequence utilities

```python
def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse FASTA file into sequence dictionary."""

def write_fasta(sequences: Dict[str, str], output_path: Path) -> None:
    """Write sequences to FASTA file."""

def translate_dna_sequence(dna_sequence: str, genetic_code: int = 11) -> str:
    """Translate DNA sequence to amino acid sequence."""
```

### ID management utilities

```python
def generate_pangenomeplus_id(
    genome_id: str,
    feature_type: FeatureType,
    index: int
) -> str:
    """Generate standardized PanGenomePlus feature ID."""

def parse_pangenomeplus_id(pgp_id: str) -> Dict[str, str]:
    """Parse PanGenomePlus ID into components."""

def create_family_id(
    feature_type: FeatureType,
    cluster_id: str
) -> str:
    """Create standardized gene family ID."""
```

## Extension Points

### Custom annotation tools

```python
class AnnotationTool:
    """Base class for custom annotation tools."""

    def __init__(self, name: str, executable: str):
        self.name = name
        self.executable = executable

    def run(
        self,
        genome_path: Path,
        output_dir: Path,
        **kwargs
    ) -> Dict[str, Path]:
        """Run annotation tool and return output file paths."""
        raise NotImplementedError

    def parse_output(self, output_file: Path) -> List[GenomeFeature]:
        """Parse tool output into GenomeFeature objects."""
        raise NotImplementedError
```

### Custom output formats

```python
class OutputFormat:
    """Base class for custom output formats."""

    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description

    def generate(
        self,
        gene_families: Dict[str, str],
        family_members: Dict[str, List[str]],
        features: List[GenomeFeature],
        output_dir: Path,
        **kwargs
    ) -> Path:
        """Generate output file in custom format."""
        raise NotImplementedError
```

### Plugin system

```python
def register_annotation_tool(tool: AnnotationTool) -> None:
    """Register custom annotation tool."""

def register_output_format(format: OutputFormat) -> None:
    """Register custom output format."""

def list_available_tools() -> List[str]:
    """List all available annotation tools."""

def list_available_formats() -> List[str]:
    """List all available output formats."""
```

## Testing Framework

### Test utilities

```python
from pangenomeplus.testing import (
    create_test_genome,
    create_test_config,
    run_pipeline_test
)

def create_test_genome(
    genome_id: str,
    n_genes: int = 100,
    output_path: Path = None
) -> Path:
    """Create synthetic genome for testing."""

def create_test_config(
    minimal: bool = True,
    **overrides
) -> Dict[str, Any]:
    """Create test configuration."""

def run_pipeline_test(
    test_genomes: List[Path],
    config: Dict[str, Any] = None,
    validate: bool = True
) -> Dict[str, Any]:
    """Run pipeline test with validation."""
```

### Test fixtures

```python
import pytest
from pangenomeplus.testing import TestFixtures

@pytest.fixture
def test_genomes():
    """Provide test genome files."""
    return TestFixtures.create_test_genomes(n_genomes=3)

@pytest.fixture
def default_config():
    """Provide default test configuration."""
    return TestFixtures.get_test_config()

def test_basic_pipeline(test_genomes, default_config, tmp_path):
    """Test basic pipeline functionality."""
    results = run_pangenome_pipeline(
        genomes_dir=test_genomes,
        output_dir=tmp_path,
        config=default_config
    )
    assert results["n_genomes"] == 3
    assert results["n_families"] > 0
```

## Development Guidelines

### Code style

- **PEP 8 compliance**: All code follows PEP 8 style guidelines
- **Type hints**: All functions include complete type annotations
- **Docstrings**: 100% function documentation coverage using Google style
- **Error handling**: Explicit exception handling with custom exception types

### Architecture principles

- **KISS (Keep It Simple, Stupid)**: Straightforward, uncomplicated solutions
- **YAGNI (You Aren't Gonna Need It)**: Implement only what's currently needed
- **SOLID principles**: Single responsibility, open-closed, etc.
- **Modular design**: Clear separation of concerns between modules

### Contributing workflow

1. **Fork and clone**: Fork repository and clone locally
2. **Branch**: Create feature branch from main
3. **Develop**: Write code following style guidelines
4. **Test**: Add tests for new functionality
5. **Document**: Update documentation and docstrings
6. **Submit**: Create pull request with detailed description

### Code review checklist

- [ ] Functions have complete type hints and docstrings
- [ ] Error cases are handled with appropriate exceptions
- [ ] New functionality includes tests
- [ ] Code follows PEP 8 style guidelines
- [ ] No unused imports or dead code
- [ ] Documentation is updated if needed

### Performance considerations

- **Memory usage**: Monitor memory consumption for large datasets
- **Thread safety**: Ensure thread-safe operations when using parallelism
- **Resource cleanup**: Properly clean up temporary files and processes
- **Scalability**: Design with large datasets (100,000+ genomes) in mind

---

## Example Integration

Here's a complete example of using the PanGenomePlus API:

```python
#!/usr/bin/env python3
"""Example PanGenomePlus integration script."""

import logging
from pathlib import Path
from pangenomeplus.pipeline import run_pangenome_pipeline
from pangenomeplus.utils.config import create_default_configuration
from pangenomeplus.core.exceptions import PipelineError

def main():
    # Setup logging
    logging.basicConfig(level=logging.INFO)

    # Configuration
    config = create_default_configuration()
    config["clustering"]["protein"]["identity"] = 0.9
    config["clustering"]["gpu"]["use_gpu"] = True
    config["output"]["formats"] = ["transformer", "presence_absence", "roary"]
    config["resources"]["threads"] = 8

    # Input/output paths
    genomes_dir = Path("./input_genomes")
    output_dir = Path("./pangenome_results")

    try:
        # Run pipeline
        results = run_pangenome_pipeline(
            genomes_dir=genomes_dir,
            output_dir=output_dir,
            config=config,
            log_level="INFO"
        )

        # Process results
        print(f"[OK] Analysis completed successfully!")
        print(f"  Genomes processed: {results['n_genomes']}")
        print(f"  Features extracted: {results['n_features']}")
        print(f"  Gene families: {results['n_families']}")
        print(f"  Runtime: {results['runtime_formatted']}")

        # Access output files
        output_files = results["output_files"]
        print(f"  Transformer format: {output_files['transformer_format']}")
        print(f"  Presence/absence matrix: {output_files['presence_absence_matrix']}")

    except PipelineError as e:
        print(f"[ERROR] Pipeline failed: {e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())
```

This API reference provides complete documentation for developers to use PanGenomePlus programmatically or contribute to its development.