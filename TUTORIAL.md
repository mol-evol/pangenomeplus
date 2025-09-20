# PanGenomePlus Tutorial

A comprehensive step-by-step guide to performing pangenome analysis with PanGenomePlus.

## Table of Contents

1. [Prerequisites and Installation](#prerequisites-and-installation)
2. [First Analysis: Small Dataset](#first-analysis-small-dataset)
3. [Understanding the Output](#understanding-the-output)
4. [Advanced Configuration](#advanced-configuration)
5. [Large Dataset Analysis](#large-dataset-analysis)
6. [Downstream Analysis](#downstream-analysis)
7. [Troubleshooting Common Issues](#troubleshooting-common-issues)
8. [Best Practices](#best-practices)

## Prerequisites and Installation

### Step 1: Install External Dependencies

**Required tools that must be in your PATH:**

#### MMseqs2 (Required)
```bash
# Option 1: Using conda (recommended)
conda install -c conda-forge mmseqs2

# Option 2: Using Homebrew (macOS)
brew install mmseqs2

# Option 3: From source (see MMseqs2 documentation)
```

#### Prodigal (Required)
```bash
# Using conda
conda install -c bioconda prodigal

# Using Homebrew
brew install prodigal

# Verify installation
prodigal -v
```

#### Optional Tools
```bash
# tRNAscan-SE for tRNA detection
conda install -c bioconda trnascan-se

# Barrnap for rRNA detection
conda install -c bioconda barrnap

# MINCED for CRISPR detection
conda install -c bioconda minced
```

### Step 2: Install PanGenomePlus

```bash
# Clone the repository
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus

# Create and activate virtual environment
python -m venv pangenome_env
source pangenome_env/bin/activate  # On Windows: pangenome_env\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Install PanGenomePlus in development mode
pip install -e .
```

### Step 3: Verify Installation

```bash
# Test external tools
mmseqs version
prodigal -v

# Test PanGenomePlus
python -m pangenomeplus --help
```

## First Analysis: Small Dataset

Let's start with a small example using 3-5 bacterial genomes.

### Step 1: Prepare Your Data

Create a directory structure:
```bash
mkdir my_pangenome_analysis
cd my_pangenome_analysis
mkdir genomes
```

**Genome file requirements:**
- FASTA format (`.fasta`, `.fa`, `.fna`)
- One genome per file
- Files can be gzipped (`.gz`)
- Filenames will be used as genome identifiers

Example file structure:
```
genomes/
├── E_coli_K12.fasta
├── E_coli_O157.fasta
├── E_coli_CFT073.fasta
└── Salmonella_LT2.fasta
```

### Step 2: Run Basic Analysis

```bash
# Basic command - this will take 5-30 minutes depending on genome sizes
python -m pangenomeplus genomes/ results/
```

**What happens during analysis:**
1. **Validation**: Checks input files and external tools
2. **Annotation**: Runs Prodigal for gene prediction
3. **Feature Extraction**: Extracts proteins, tRNA, rRNA, intergenic regions
4. **Clustering**: Groups similar features using MMseqs2
5. **Family Assignment**: Creates gene families from clusters
6. **Output Generation**: Creates analysis files

### Step 3: Monitor Progress

PanGenomePlus provides detailed logging:
```bash
# For more verbose output
python -m pangenomeplus genomes/ results/ --verbose

# To see debug information
python -m pangenomeplus genomes/ results/ --log-level DEBUG
```

You'll see progress messages like:
```
2024-01-15 10:30:15 - INFO - Starting pangenome pipeline
2024-01-15 10:30:15 - INFO - Found 4 genome files
2024-01-15 10:30:16 - INFO - Starting genome annotation...
2024-01-15 10:30:45 - INFO - Annotation complete for 4 genomes
2024-01-15 10:30:45 - INFO - Starting clustering...
2024-01-15 10:31:20 - INFO - Created 3,456 gene families
2024-01-15 10:31:25 - INFO - Pipeline completed successfully
```

## Understanding the Output

After analysis completes, examine the results directory:

```bash
results/
├── annotation/           # Raw annotation files
│   ├── E_coli_K12/
│   ├── E_coli_O157/
│   └── ...
├── clustering/           # MMseqs2 clustering results
│   ├── protein_coding_gene/
│   ├── trna/
│   └── ...
├── families/            # Gene family assignments
│   ├── gene_to_family.tsv
│   ├── family_summary.tsv
│   └── family_members.json
├── features/            # Extracted genomic features
│   └── feature_summary.tsv
└── outputs/             # Final analysis results
    ├── transformer_format.txt
    ├── presence_absence_matrix.tsv
    └── ...
```

### Key Output Files

#### 1. Transformer Format (`outputs/transformer_format.txt`)
AI-ready format with feature type prefixes:
```
# Example content:
E_coli_K12    P_FAM_000001 P_FAM_000002 T_FAM_000001 I_INT_000001
E_coli_O157   P_FAM_000001 P_FAM_000003 T_FAM_000001 I_INT_000002
```

#### 2. Presence/Absence Matrix (`outputs/presence_absence_matrix.tsv`)
Binary matrix showing which families are present in each genome:
```
Family_ID       E_coli_K12    E_coli_O157    Salmonella_LT2
P_FAM_000001    1             1              1
P_FAM_000002    1             0              1
P_FAM_000003    0             1              0
```

#### 3. Family Summary (`families/family_summary.tsv`)
Statistics for each gene family:
```
family_id       size    feature_types          representative
P_FAM_000001    4       protein_coding_gene    PGP_E_coli_K12_001
P_FAM_000002    2       protein_coding_gene    PGP_E_coli_K12_045
```

### Step 4: Quick Analysis

```bash
# Count total families
wc -l results/families/family_summary.tsv

# Count core genes (present in all genomes)
python -c "
import pandas as pd
df = pd.read_csv('results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
core_genes = (df.sum(axis=1) == len(df.columns)).sum()
print(f'Core genes: {core_genes}')
print(f'Total families: {len(df)}')
print(f'Core percentage: {core_genes/len(df)*100:.1f}%')
"
```

## Advanced Configuration

### Step 1: Generate Configuration Template

```bash
python -m pangenomeplus --init-config my_config.yaml
```

### Step 2: Customize Configuration

Edit `my_config.yaml`:

```yaml
# Example advanced configuration
clustering:
  protein:
    identity: 0.9        # Higher stringency
    coverage: 0.9
    sensitivity: 7.5     # Higher sensitivity
  trna:
    identity: 0.95       # Very high for tRNA
  gpu:
    use_gpu: true        # Enable GPU if available
    memory_fraction: 0.8

annotation:
  tools: ["prodigal", "trnascan", "barrnap", "minced"]  # All tools
  prodigal:
    genetic_code: 11     # Bacterial genetic code
  intergenic:
    extract: true
    min_length: 100      # Longer intergenic regions

classification:
  hard_core_threshold: 1.0    # Must be in 100% genomes
  soft_core_threshold: 0.99   # 99% for soft core
  shell_threshold: 0.15
  cloud_threshold: 0.15

output:
  formats: ["transformer", "presence_absence", "roary"]
  transformer:
    include_singletons: true
    feature_type_prefix: true

resources:
  threads: 8
  memory_limit: "32G"

validation:
  input_files: true
  system_requirements: true
```

### Step 3: Run with Custom Configuration

```bash
python -m pangenomeplus genomes/ results_advanced/ --config my_config.yaml
```

## Large Dataset Analysis

For datasets with 100+ genomes, consider these optimizations:

### Performance Tuning

```bash
# High-performance configuration for large datasets
python -m pangenomeplus genomes/ results_large/ \
  --sensitivity 6.0 \
  --use-gpu \
  --threads 16 \
  --memory-limit 64G \
  --no-validate
```

### Memory Management

```bash
# For memory-constrained systems
python -m pangenomeplus genomes/ results/ \
  --sensitivity 4.0 \
  --threads 4 \
  --memory-limit 16G
```

### Batch Processing

For very large datasets (1000+ genomes), consider processing in batches:

```bash
# Process first 100 genomes
mkdir genomes_batch1
cp genomes/genome_{001..100}.fasta genomes_batch1/
python -m pangenomeplus genomes_batch1/ results_batch1/

# Process second 100 genomes
mkdir genomes_batch2
cp genomes/genome_{101..200}.fasta genomes_batch2/
python -m pangenomeplus genomes_batch2/ results_batch2/
```

### Checkpoint Recovery

If analysis is interrupted:
```bash
# Resume from where it left off
python -m pangenomeplus genomes/ results/ --resume

# Force restart
python -m pangenomeplus genomes/ results/ --no-resume
```

## Downstream Analysis

### Feature-Type Analysis

Analyze results by genomic feature type:

```bash
# Run downstream analysis
python -m pangenomeplus --downstream-analysis results/ --output feature_analysis/
```

This creates detailed analysis by feature type:
```
feature_analysis/
├── feature_type_analysis.json      # Statistics by feature type
├── core_shell_cloud_analysis.json  # Classifications
└── pangenome_composition_plot.png  # Visualization (if matplotlib available)
```

### Custom Analysis Scripts

#### Extract Core Genome
```python
import pandas as pd

# Load presence/absence matrix
matrix = pd.read_csv('results/outputs/presence_absence_matrix.tsv',
                     sep='\t', index_col=0)

# Find core genes (present in all genomes)
core_genes = matrix[matrix.sum(axis=1) == len(matrix.columns)]
print(f"Core genome: {len(core_genes)} families")

# Save core gene list
core_genes.index.to_series().to_csv('core_genes.txt', index=False, header=False)
```

#### Analyze Gene Family Sizes
```python
import pandas as pd
import matplotlib.pyplot as plt

# Load family summary
families = pd.read_csv('results/families/family_summary.tsv', sep='\t')

# Plot family size distribution
plt.figure(figsize=(10, 6))
plt.hist(families['size'], bins=50, alpha=0.7)
plt.xlabel('Family Size (number of genes)')
plt.ylabel('Number of Families')
plt.title('Gene Family Size Distribution')
plt.yscale('log')
plt.savefig('family_size_distribution.png', dpi=300, bbox_inches='tight')
```

## Troubleshooting Common Issues

### Issue 1: MMseqs2 Not Found
```
Error: MMseqs2 not found in PATH
```

**Solution:**
```bash
# Check if MMseqs2 is installed
which mmseqs

# If not found, install it
conda install -c conda-forge mmseqs2

# Add to PATH if needed
export PATH="/path/to/mmseqs:$PATH"
```

### Issue 2: Out of Memory
```
Error: Process killed (out of memory)
```

**Solutions:**
```bash
# Reduce memory usage
python -m pangenomeplus genomes/ results/ --memory-limit 8G --threads 2

# Use lower sensitivity (faster, less memory)
python -m pangenomeplus genomes/ results/ --sensitivity 4.0

# Process smaller batches
```

### Issue 3: GPU Errors
```
Error: CUDA out of memory
```

**Solutions:**
```bash
# Reduce GPU memory usage
python -m pangenomeplus genomes/ results/ --gpu-memory 8G

# Disable GPU
python -m pangenomeplus genomes/ results/ --no-gpu
```

### Issue 4: Long Runtime
```
Analysis taking too long...
```

**Solutions:**
```bash
# Use lower sensitivity for speed
python -m pangenomeplus genomes/ results/ --sensitivity 4.0

# Skip validation
python -m pangenomeplus genomes/ results/ --no-validate

# Use more threads
python -m pangenomeplus genomes/ results/ --threads 16
```

### Issue 5: Empty Results
```
Warning: No features extracted
```

**Check:**
1. Input files are valid FASTA format
2. Genomes contain protein-coding sequences
3. Prodigal is working: `prodigal -v`
4. Check log files for specific errors

## Best Practices

### 1. Data Preparation
- Use consistent naming conventions for genome files
- Ensure high-quality assemblies (N50 > 100kb recommended)
- Include representative genomes from your study group
- Consider genome completeness (>95% recommended)

### 2. Parameter Selection

**For closely related genomes** (same species):
```bash
--protein-identity 0.95 --sensitivity 7.5
```

**For distantly related genomes** (different genera):
```bash
--protein-identity 0.7 --sensitivity 6.0
```

**For exploratory analysis** (unknown diversity):
```bash
--protein-identity 0.8 --sensitivity 7.5  # Default values
```

### 3. Quality Control

Always examine:
- Number of genes per genome (should be consistent within species)
- Core genome size (typically 60-90% for bacterial species)
- Singleton rate (high rates may indicate contamination)

```bash
# Quick quality check
python -c "
import pandas as pd
matrix = pd.read_csv('results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
print('Genomes:', len(matrix.columns))
print('Total families:', len(matrix))
print('Core families:', (matrix.sum(axis=1) == len(matrix.columns)).sum())
print('Singletons:', (matrix.sum(axis=1) == 1).sum())
"
```

### 4. Performance Optimization

- Start with default parameters
- For large datasets (>100 genomes), enable GPU if available
- Monitor memory usage and adjust `--memory-limit` accordingly
- Use checkpointing for long-running analyses

### 5. Reproducibility

- Save configuration files used for analysis
- Document MMseqs2 and tool versions
- Keep detailed analysis logs

```bash
# Save environment info
python -m pangenomeplus --version > analysis_info.txt
mmseqs version >> analysis_info.txt
prodigal -v 2>> analysis_info.txt
```

## Next Steps

After completing this tutorial, you can:

1. **Explore the API**: See [API.md](API.md) for programmatic usage
2. **Review Examples**: Check [EXAMPLES.md](EXAMPLES.md) for specific use cases
3. **Advanced Installation**: See [INSTALLATION.md](INSTALLATION.md) for HPC setups
4. **Contribute**: Submit issues or improvements to the GitHub repository

---

**Need help?** Check the [main README](README.md) or open an issue on [GitHub](https://github.com/mol-evol/pangenomeplus).