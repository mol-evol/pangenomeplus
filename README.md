# PanGenomePlus

A high-performance pangenome analysis pipeline built on MMseqs2 with support for multiple genomic feature types and AI-ready output formats.

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Overview

PanGenomePlus is a streamlined pangenome analysis pipeline that leverages MMseqs2's powerful clustering capabilities to analyze protein-coding genes, tRNA, rRNA, CRISPR arrays, and intergenic regions. The pipeline follows KISS (Keep It Simple, Stupid) principles while maintaining professional code quality with 100% docstring coverage.

### Key Features

- **🚀 High Performance**: Built on MMseqs2 easy-cluster for optimal scaling from 2 to 500,000+ sequences
- **🧬 Multi-Feature Analysis**: Analyzes protein-coding genes, tRNA, rRNA, CRISPR arrays, and intergenic regions
- **🤖 AI-Ready Outputs**: Transformer format optimized for machine learning applications
- **📊 Multiple Output Formats**: Transformer, presence/absence matrices, Roary-compatible, and FASTA
- **⚡ GPU Acceleration**: Optional GPU support for large-scale analyses
- **🔧 Modular Architecture**: Clean separation of annotation → clustering → families → output
- **📋 Professional Quality**: 100% function documentation and PEP 8 compliance

## Quick Start

```bash
# Basic pangenome analysis
python -m pangenomeplus genomes/ output/

# High-sensitivity analysis with GPU acceleration
python -m pangenomeplus genomes/ output/ \
  --sensitivity 7.5 \
  --protein-identity 0.9 \
  --use-gpu \
  --threads 8

# Generate configuration template
python -m pangenomeplus --init-config my_config.yaml

# Analyze dataset characteristics only
python -m pangenomeplus genomes/ --analyze-only

# Run downstream feature-type analysis
python -m pangenomeplus --downstream-analysis results/ --output analysis/
```

## Installation

### Prerequisites

**External Tools** (must be installed and in PATH):
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (required)
- [Prodigal](https://github.com/hyattpd/Prodigal) (required)
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) (optional)
- [Barrnap](https://github.com/tseemann/barrnap) (optional)
- [MINCED](https://github.com/ctSkennerton/minced) (optional)

### Install PanGenomePlus

```bash
# Clone repository
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Quick Installation Check

```bash
# Verify external tools
mmseqs version
prodigal -v

# Test PanGenomePlus
python -m pangenomeplus --help
```

## Usage

### Basic Analysis

```bash
# Minimal command - analyzes genomes in 'genomes/' directory
python -m pangenomeplus genomes/ output/
```

**Input Requirements**:
- Directory containing genome FASTA files (`.fasta`, `.fa`, `.fna`)
- Each file should contain one genome
- Files can be gzipped (`.gz`)

### Advanced Options

#### Clustering Parameters
```bash
# High-sensitivity clustering
--sensitivity 7.5              # MMseqs2 sensitivity (1.0-8.0, default: 7.5)
--protein-identity 0.9         # Protein identity threshold (default: 0.8)
--protein-coverage 0.8         # Coverage threshold (default: 0.8)
--use-gpu                      # Enable GPU acceleration
--gpu-memory 16G              # GPU memory limit
```

#### Annotation Control
```bash
# Customize annotation tools
--annotation-tools prodigal,trnascan,barrnap,minced
--genetic-code 11             # Genetic code table (default: 11 = bacterial)
```

#### Output Formats
```bash
# Select output formats
--output-formats transformer,presence_absence,roary,fasta
--compress-output            # Compress large files
```

#### Pipeline Control
```bash
# Stop after specific stage
--stop-after annotation      # Options: annotation, clustering, families, analysis
--no-validate               # Skip input validation for faster startup
--resume                    # Resume from checkpoints (default: enabled)
--log-level DEBUG          # Logging verbosity
```

### Configuration Files

Generate a configuration template:
```bash
python -m pangenomeplus --init-config config.yaml
```

Example configuration:
```yaml
clustering:
  protein:
    identity: 0.8
    coverage: 0.8
    sensitivity: 7.5
  gpu:
    use_gpu: false
    memory_fraction: 0.8

annotation:
  tools: ["prodigal", "trnascan", "barrnap"]
  genetic_code: 11

output:
  formats: ["transformer", "presence_absence"]
  compress: false

resources:
  threads: 4
  memory_limit: null
```

Use configuration:
```bash
python -m pangenomeplus genomes/ output/ --config config.yaml
```

## Output Files

### Primary Outputs

| File | Description | Format |
|------|-------------|---------|
| `transformer_format.txt` | AI-ready sequence format | Tab-separated |
| `presence_absence_matrix.tsv` | Binary presence/absence matrix | TSV |
| `gene_to_family.tsv` | Gene-to-family mappings | TSV |
| `family_summary.tsv` | Family statistics | TSV |
| `pangenome_stats.json` | Analysis statistics | JSON |

### Transformer Format

The transformer format uses feature type prefixes for genomic element identification:

```
# Feature type prefixes:
P_  = Protein-coding genes
I_  = Intergenic regions
T_  = tRNA genes
R_  = rRNA genes
C_  = CRISPR arrays
S_  = Singleton genes

# Example output:
genome_A    P_FAM_000001 T_FAM_000002 I_INT_000001 R_FAM_000003
genome_B    P_FAM_000001 I_INT_000002 C_FAM_000004 P_FAM_000002 S_PGP_genome_B_001
```

### Roary-Compatible Output

Compatible with [Roary](https://sanger-pathogens.github.io/Roary/) for downstream analysis:
- `gene_presence_absence.csv` - Main Roary output format
- `summary_statistics.txt` - Pangenome statistics
- `clustered_proteins.fasta` - Representative sequences

## Performance & Scaling

### Dataset Size Guidelines

| Genomes | Memory (RAM) | Time* | Recommendation |
|---------|--------------|-------|----------------|
| 2-10 | 4-8 GB | 5-30 min | Standard settings |
| 10-100 | 8-32 GB | 30 min - 2 hours | Increase threads |
| 100-1000 | 32-128 GB | 2-8 hours | Consider GPU |
| 1000+ | 128+ GB | 8+ hours | GPU recommended |

*Approximate times on modern hardware

### Performance Optimization

```bash
# For large datasets (100+ genomes)
python -m pangenomeplus genomes/ output/ \
  --sensitivity 6.0 \
  --use-gpu \
  --threads 16 \
  --memory-limit 64G

# For maximum speed (lower sensitivity)
python -m pangenomeplus genomes/ output/ \
  --sensitivity 4.0 \
  --protein-identity 0.7

# For maximum sensitivity (slower but more thorough)
python -m pangenomeplus genomes/ output/ \
  --sensitivity 8.0 \
  --protein-identity 0.9
```

## Advanced Features

### Downstream Analysis

Analyze completed results by feature type:
```bash
python -m pangenomeplus --downstream-analysis results/ --output feature_analysis/
```

This generates:
- Feature-type separated statistics
- Core/shell/cloud classifications per feature type
- Visualization plots (if matplotlib available)

### Pangenome Classifications

- **Core genes**: Present in ≥95% of genomes
- **Shell genes**: Present in 15-95% of genomes
- **Cloud genes**: Present in <15% of genomes
- **Singletons**: Present in only one genome

Customize thresholds:
```bash
--core-threshold 0.99    # 99% for core
--cloud-threshold 0.05   # 5% for cloud
```

### Checkpointing

Pipeline automatically saves progress and can resume:
```bash
# Pipeline interrupted? Resume with:
python -m pangenomeplus genomes/ output/ --resume

# Force restart:
python -m pangenomeplus genomes/ output/ --no-resume
```

## Troubleshooting

### Common Issues

**1. MMseqs2 not found**
```bash
# Install MMseqs2
conda install -c conda-forge mmseqs2
# or
brew install mmseqs2
```

**2. Out of memory**
```bash
# Reduce memory usage
python -m pangenomeplus genomes/ output/ --memory-limit 16G --threads 4
```

**3. GPU errors**
```bash
# Disable GPU if causing issues
python -m pangenomeplus genomes/ output/ --no-gpu
```

**4. Large dataset slow**
```bash
# Use lower sensitivity for speed
python -m pangenomeplus genomes/ output/ --sensitivity 4.0
```

### Getting Help

- Check the [TUTORIAL.md](TUTORIAL.md) for step-by-step guidance
- See [INSTALLATION.md](INSTALLATION.md) for detailed setup instructions
- Review [EXAMPLES.md](EXAMPLES.md) for practical use cases
- Check [API.md](API.md) for technical documentation

## Architecture

```
Input Genomes
     ↓
Annotation (Prodigal, tRNAscan, Barrnap, MINCED)
     ↓
Feature Extraction (Proteins, tRNA, rRNA, CRISPR, Intergenic)
     ↓
Clustering (MMseqs2 easy-cluster)
     ↓
Family Assignment
     ↓
Output Generation (Transformer, Presence/Absence, Roary, FASTA)
```

### Code Quality

- **100% Function Documentation**: Every function has comprehensive docstrings
- **PEP 8 Compliance**: Consistent code formatting
- **Modular Design**: Clear separation of concerns
- **KISS Principles**: Simple, maintainable code
- **Type Hints**: Full type annotation support

## Citation

If you use PanGenomePlus in your research, please cite:

```
McInerney et al. (2024). PanGenomePlus: High-performance pangenome analysis
with multi-feature support. GitHub: https://github.com/mol-evol/pangenomeplus
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Acknowledgments

- **MMseqs2** team for the excellent clustering software
- **Prodigal** developers for gene prediction
- **tRNAscan-SE**, **Barrnap**, and **MINCED** teams for feature annotation tools
- Contributors and users of the PanGenomePlus pipeline

---

**Repository**: https://github.com/mol-evol/pangenomeplus
**Author**: James McInerney (james.mcinerney@nottingham.ac.uk)
**Version**: 1.0.0