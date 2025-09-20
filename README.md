# PanGenomePlus

An adaptive scale pangenome analysis pipeline that automatically adapts to dataset size and computational resources.

## Features

- **Adaptive Strategy Selection**: Automatically chooses optimal processing strategy based on dataset size
- **MMseqs2 Clustering**: High-performance protein clustering with GPU acceleration support
- **Multi-Feature Analysis**: Analyzes protein-coding genes, tRNA, rRNA, CRISPR arrays, and intergenic regions
- **Transformer Format Output**: AI-ready output format for machine learning applications
- **HPC Support**: GPU acceleration and memory optimization for large datasets
- **Feature-Type Separated Analysis**: Distinguishes between different genomic feature types

## Quick Start

```bash
# Basic analysis
pangenomeplus run genomes/ output/

# With custom parameters
pangenomeplus run genomes/ output/ \
  --protein-identity 0.9 \
  --sensitivity 7.5 \
  --use-gpu \
  --threads 8

# Analyze dataset characteristics
pangenomeplus analyze genomes/

# Generate configuration file
pangenomeplus init-config --output my_config.yaml
```

## CLI Parameters

### Clustering Parameters
- `--protein-identity`: Protein clustering identity threshold (0.0-1.0, default: 0.8)
- `--protein-coverage`: Protein clustering coverage threshold (0.0-1.0, default: 0.8)
- `--sensitivity`: MMseqs2 search sensitivity (1.0-8.0, default: 7.5)
- `--use-gpu`: Enable GPU acceleration for clustering
- `--gpu-memory`: GPU memory limit (e.g., 8G, 16G)

### Annotation Parameters
- `--genetic-code`: Genetic code table for Prodigal (1-31, default: 11)
- `--annotation-tools`: Comma-separated annotation tools (prodigal,trnascan,barrnap,minced)

### Classification Parameters
- `--core-threshold`: Core gene threshold (0.0-1.0, default: 0.95)
- `--cloud-threshold`: Cloud gene threshold (0.0-1.0, default: 0.15)

### Output Parameters
- `--output-formats`: Comma-separated output formats (transformer,presence_absence,roary,fasta)
- `--compress-output`: Compress large output files

## Transformer Format

The transformer format uses feature type prefixes for genomic element identification:

- **P_** = Protein-coding genes
- **I_** = Intergenic regions
- **T_** = tRNA genes
- **R_** = rRNA genes
- **C_** = CRISPR arrays

Example output:
```
genome_A P_FAM_001 T_FAM_002 I_INT_001 R_FAM_003
genome_B P_FAM_001 I_INT_002 C_FAM_004 P_FAM_002
```

## Installation

```bash
# Clone repository
git clone https://github.com/user/pangenomeplus.git
cd pangenomeplus

# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

## Requirements

- Python 3.8+
- MMseqs2
- Prodigal
- tRNAscan-SE (optional)
- Barrnap (optional)
- MINCED (optional)

## Strategy Selection

PanGenomePlus automatically selects processing strategies:

- **Direct** (≤100 genomes): Direct clustering approach
- **Batched** (≤1000 genomes): Batched hierarchical clustering
- **Hierarchical** (>1000 genomes): Distributed hierarchical clustering

## Configuration

Configuration files support YAML or JSON format. Generate a template:

```bash
pangenomeplus init-config --output config.yaml
```

Key configuration sections:
- `clustering`: MMseqs2 parameters and GPU settings
- `annotation`: Tool selection and parameters
- `classification`: Core/shell/cloud thresholds
- `output`: Format selection and compression
- `resources`: Thread and memory limits

## Output Files

- `transformer.txt`: AI-ready format for machine learning
- `presence_absence.tsv`: Binary presence/absence matrix
- `family_summary.tsv`: Gene family statistics
- `gene_to_family.tsv`: Gene-to-family mappings
- `pangenome_stats.json`: Analysis statistics

## Advanced Analysis

```bash
# Feature-type separated analysis
pangenomeplus feature-analysis results/ --output feature_analysis/

# Validate inputs before running
pangenomeplus validate genomes/ --config config.yaml
```