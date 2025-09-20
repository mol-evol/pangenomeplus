# PanGenomePlus Examples

Practical examples demonstrating real-world usage scenarios for PanGenomePlus across different research contexts and data types.

## Table of Contents

1. [Basic Examples](#basic-examples)
2. [Research Scenarios](#research-scenarios)
3. [Parameter Optimization](#parameter-optimization)
4. [Integration Examples](#integration-examples)
5. [Large-Scale Analyses](#large-scale-analyses)
6. [Comparative Studies](#comparative-studies)
7. [Troubleshooting Examples](#troubleshooting-examples)
8. [Advanced Workflows](#advanced-workflows)

## Basic Examples

### Example 1: First Analysis with E. coli Genomes

```bash
# Download E. coli genomes (example URLs)
mkdir ecoli_genomes
cd ecoli_genomes

# Download representative E. coli strains
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2/GCF_000005825.2_ASM582v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/485/GCF_000010485.1_ASM1048v1/GCF_000010485.1_ASM1048v1_genomic.fna.gz

# Rename for clarity
mv GCF_000005825.2_ASM582v2_genomic.fna.gz E_coli_K12_MG1655.fasta.gz
mv GCF_000008865.2_ASM886v2_genomic.fna.gz E_coli_O157_H7_EDL933.fasta.gz
mv GCF_000010485.1_ASM1048v1_genomic.fna.gz E_coli_CFT073.fasta.gz

cd ..

# Basic analysis
python -m pangenomeplus ecoli_genomes/ ecoli_results/

# Check results
ls ecoli_results/outputs/
# Expected files:
# - transformer_format.txt
# - presence_absence_matrix.tsv
# - gene_to_family.tsv
# - family_summary.tsv
# - pangenome_stats.json
```

### Example 2: Custom Configuration

```bash
# Generate configuration template
python -m pangenomeplus --init-config ecoli_config.yaml

# Edit configuration for higher sensitivity
cat > ecoli_config.yaml << EOF
clustering:
  protein:
    identity: 0.9
    coverage: 0.9
    sensitivity: 7.5
  gpu:
    use_gpu: false

annotation:
  tools: ["prodigal", "trnascan", "barrnap"]
  genetic_code: 11

output:
  formats: ["transformer", "presence_absence", "roary"]
  compress: false

resources:
  threads: 8
  memory_limit: "32G"
EOF

# Run with custom configuration
python -m pangenomeplus ecoli_genomes/ ecoli_results_custom/ --config ecoli_config.yaml
```

### Example 3: GPU-Accelerated Analysis

```bash
# High-performance analysis with GPU
python -m pangenomeplus ecoli_genomes/ ecoli_results_gpu/ \
  --sensitivity 7.5 \
  --protein-identity 0.9 \
  --use-gpu \
  --threads 16 \
  --verbose
```

## Research Scenarios

### Scenario 1: Species-Level Pangenome Analysis

**Research Question**: What is the core and accessory genome of Salmonella enterica?

```bash
# Prepare Salmonella genomes
mkdir salmonella_analysis
cd salmonella_analysis

# Download 20+ Salmonella strains (example workflow)
# ... download commands ...

# Species-level analysis with optimized parameters
python -m pangenomeplus salmonella_genomes/ salmonella_pangenome/ \
  --sensitivity 7.5 \
  --protein-identity 0.8 \
  --protein-coverage 0.8 \
  --threads 8 \
  --log-level INFO

# Analyze results
cd salmonella_pangenome/outputs/

# Count core genes (present in all genomes)
python3 << EOF
import pandas as pd

# Load presence/absence matrix
matrix = pd.read_csv('presence_absence_matrix.tsv', sep='\t', index_col=0)
n_genomes = len(matrix.columns)

# Calculate pangenome statistics
core_genes = (matrix.sum(axis=1) == n_genomes).sum()
accessory_genes = (matrix.sum(axis=1) < n_genomes).sum()
singletons = (matrix.sum(axis=1) == 1).sum()

print(f"Salmonella enterica pangenome analysis:")
print(f"  Total genomes: {n_genomes}")
print(f"  Total gene families: {len(matrix)}")
print(f"  Core genes: {core_genes} ({core_genes/len(matrix)*100:.1f}%)")
print(f"  Accessory genes: {accessory_genes} ({accessory_genes/len(matrix)*100:.1f}%)")
print(f"  Singleton genes: {singletons} ({singletons/len(matrix)*100:.1f}%)")
EOF
```

### Scenario 2: Strain-Level Comparison within Species

**Research Question**: How do antibiotic-resistant and sensitive strains differ?

```bash
# Organize genomes by phenotype
mkdir resistance_study
cd resistance_study
mkdir resistant_strains sensitive_strains

# Copy genomes to appropriate directories
# cp path/to/resistant/*.fasta resistant_strains/
# cp path/to/sensitive/*.fasta sensitive_strains/

# Analyze each group separately
python -m pangenomeplus resistant_strains/ resistant_results/ \
  --sensitivity 8.0 \
  --protein-identity 0.95

python -m pangenomeplus sensitive_strains/ sensitive_results/ \
  --sensitivity 8.0 \
  --protein-identity 0.95

# Combined analysis
mkdir all_strains
cp resistant_strains/*.fasta all_strains/
cp sensitive_strains/*.fasta all_strains/

python -m pangenomeplus all_strains/ combined_results/ \
  --sensitivity 8.0 \
  --protein-identity 0.95

# Compare results
python3 << EOF
import pandas as pd
import numpy as np

# Load matrices
resistant = pd.read_csv('resistant_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
sensitive = pd.read_csv('sensitive_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
combined = pd.read_csv('combined_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)

# Find resistance-associated genes
resistant_genomes = [col for col in combined.columns if 'resistant' in col.lower()]
sensitive_genomes = [col for col in combined.columns if 'sensitive' in col.lower()]

# Calculate prevalence in each group
resistant_prevalence = combined[resistant_genomes].mean(axis=1)
sensitive_prevalence = combined[sensitive_genomes].mean(axis=1)

# Find genes enriched in resistant strains
enriched_in_resistant = combined[
    (resistant_prevalence > 0.8) & (sensitive_prevalence < 0.2)
].index

print(f"Genes enriched in resistant strains: {len(enriched_in_resistant)}")
for gene in enriched_in_resistant[:10]:  # Show first 10
    print(f"  {gene}: {resistant_prevalence[gene]:.2f} vs {sensitive_prevalence[gene]:.2f}")
EOF
```

### Scenario 3: Evolutionary Analysis Across Genera

**Research Question**: Compare pangenome evolution between related genera

```bash
# Multi-genus comparison
mkdir evolutionary_study
cd evolutionary_study

# Organize by genus
mkdir escherichia shigella salmonella

# Analyze each genus
for genus in escherichia shigella salmonella; do
    echo "Analyzing $genus..."
    python -m pangenomeplus ${genus}/ ${genus}_results/ \
      --sensitivity 6.0 \
      --protein-identity 0.7 \
      --threads 8
done

# Cross-genus analysis
mkdir all_genera
cp escherichia/*.fasta all_genera/
cp shigella/*.fasta all_genera/
cp salmonella/*.fasta all_genera/

python -m pangenomeplus all_genera/ cross_genus_results/ \
  --sensitivity 6.0 \
  --protein-identity 0.7 \
  --threads 8

# Analyze evolutionary patterns
python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt

# Load results
esc_matrix = pd.read_csv('escherichia_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
shi_matrix = pd.read_csv('shigella_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
sal_matrix = pd.read_csv('salmonella_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)

# Calculate pangenome sizes
genera_stats = {
    'Escherichia': {'genomes': len(esc_matrix.columns), 'families': len(esc_matrix)},
    'Shigella': {'genomes': len(shi_matrix.columns), 'families': len(shi_matrix)},
    'Salmonella': {'genomes': len(sal_matrix.columns), 'families': len(sal_matrix)}
}

for genus, stats in genera_stats.items():
    ratio = stats['families'] / stats['genomes']
    print(f"{genus}: {stats['families']} families / {stats['genomes']} genomes = {ratio:.1f} families per genome")
EOF
```

## Parameter Optimization

### High Sensitivity Analysis (Research Grade)

```bash
# Maximum sensitivity for publication-quality results
python -m pangenomeplus genomes/ high_sens_results/ \
  --sensitivity 8.0 \
  --protein-identity 0.95 \
  --protein-coverage 0.9 \
  --threads 16 \
  --memory-limit 64G \
  --output-formats transformer,presence_absence,roary,fasta \
  --log-level DEBUG
```

### Fast Preliminary Analysis

```bash
# Quick analysis for initial exploration
python -m pangenomeplus genomes/ quick_results/ \
  --sensitivity 4.0 \
  --protein-identity 0.7 \
  --protein-coverage 0.6 \
  --threads 4 \
  --no-validate \
  --output-formats transformer,presence_absence
```

### Memory-Efficient Analysis

```bash
# For memory-constrained systems
python -m pangenomeplus genomes/ memory_efficient/ \
  --sensitivity 5.0 \
  --protein-identity 0.8 \
  --threads 2 \
  --memory-limit 8G \
  --annotation-tools prodigal  # Only essential tools
```

## Integration Examples

### Example 1: Integration with Roary Workflow

```bash
# Generate Roary-compatible output
python -m pangenomeplus genomes/ pangenome_results/ \
  --output-formats roary \
  --threads 8

# Use Roary outputs for downstream analysis
cd pangenome_results/outputs/

# Roary files are compatible with existing Roary tools
roary_plots.py gene_presence_absence.csv

# Phylogenetic analysis with FastTree
FastTree -nt core_gene_alignment.aln > core_genome_tree.newick
```

### Example 2: Machine Learning Integration

```python
#!/usr/bin/env python3
"""Example: Using PanGenomePlus output for machine learning."""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# Load PanGenomePlus transformer format
data = []
genome_names = []

with open('pangenome_results/outputs/transformer_format.txt', 'r') as f:
    for line in f:
        if line.strip():
            parts = line.strip().split('\t')
            genome_names.append(parts[0])
            features = parts[1:] if len(parts) > 1 else []

            # Create binary feature vector
            all_features = set()
            for genome_line in open('pangenome_results/outputs/transformer_format.txt'):
                all_features.update(genome_line.strip().split('\t')[1:])

            feature_vector = [1 if feat in features else 0 for feat in sorted(all_features)]
            data.append(feature_vector)

# Convert to DataFrame
df = pd.DataFrame(data, index=genome_names, columns=sorted(all_features))

# PCA analysis
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df)

# Clustering
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(df)

# Visualization
plt.figure(figsize=(10, 8))
scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=clusters, cmap='viridis')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
plt.title('Pangenome-based Genome Clustering')
plt.colorbar(scatter)

# Annotate points
for i, genome in enumerate(genome_names):
    plt.annotate(genome, (pca_result[i, 0], pca_result[i, 1]),
                xytext=(5, 5), textcoords='offset points', fontsize=8)

plt.tight_layout()
plt.savefig('pangenome_pca.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Analysis complete:")
print(f"  Genomes: {len(genome_names)}")
print(f"  Features: {df.shape[1]}")
print(f"  PC1 variance explained: {pca.explained_variance_ratio_[0]:.1%}")
print(f"  PC2 variance explained: {pca.explained_variance_ratio_[1]:.1%}")
```

### Example 3: Integration with SLURM HPC

```bash
#!/bin/bash
#SBATCH --job-name=pangenome_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=compute

# Load modules
module load python/3.9
module load mmseqs2/13.45111
module load prodigal/2.6.3

# Activate environment
source ~/pangenome_env/bin/activate

# Set up parallel environment
export OMP_NUM_THREADS=16

# Run analysis
python -m pangenomeplus \
  $SLURM_SUBMIT_DIR/genomes/ \
  $SLURM_SUBMIT_DIR/results/ \
  --threads 16 \
  --memory-limit 120G \
  --sensitivity 7.5 \
  --protein-identity 0.9 \
  --log-level INFO

# Generate summary report
cd $SLURM_SUBMIT_DIR/results/outputs/
python3 << EOF
import pandas as pd
import json

# Load results
matrix = pd.read_csv('presence_absence_matrix.tsv', sep='\t', index_col=0)
stats = json.load(open('pangenome_stats.json'))

# Generate summary
n_genomes = len(matrix.columns)
n_families = len(matrix)
core_genes = (matrix.sum(axis=1) == n_genomes).sum()
accessory_genes = n_families - core_genes

summary = {
    'job_id': '$SLURM_JOB_ID',
    'genomes': n_genomes,
    'total_families': n_families,
    'core_genes': core_genes,
    'accessory_genes': accessory_genes,
    'core_percentage': round(core_genes/n_families*100, 1)
}

with open('slurm_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("SLURM job completed successfully")
print(f"Results saved to: $SLURM_SUBMIT_DIR/results/")
EOF
```

## Large-Scale Analyses

### Example 1: Processing 1000+ Genomes

```bash
# Batch processing strategy for very large datasets
mkdir large_scale_analysis
cd large_scale_analysis

# Organize genomes into manageable batches
mkdir -p batches/{batch_001..batch_010}

# Distribute genomes across batches (100 genomes per batch)
# ... distribution logic ...

# Process batches in parallel
for batch in batch_001 batch_002 batch_003; do
    echo "Processing $batch..."
    python -m pangenomeplus \
      batches/$batch/ \
      results_$batch/ \
      --sensitivity 6.0 \
      --protein-identity 0.8 \
      --use-gpu \
      --threads 8 \
      --memory-limit 32G &
done

wait  # Wait for all batches to complete

# Combine results (custom script)
python3 << EOF
import pandas as pd
from pathlib import Path

# Combine presence/absence matrices
all_matrices = []
for batch_dir in Path('.').glob('results_batch_*'):
    matrix_file = batch_dir / 'outputs' / 'presence_absence_matrix.tsv'
    if matrix_file.exists():
        matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
        all_matrices.append(matrix)

# Merge matrices
if all_matrices:
    combined_matrix = pd.concat(all_matrices, axis=1)
    combined_matrix = combined_matrix.fillna(0).astype(int)

    # Save combined results
    combined_matrix.to_csv('combined_presence_absence_matrix.tsv', sep='\t')

    print(f"Combined analysis:")
    print(f"  Total genomes: {len(combined_matrix.columns)}")
    print(f"  Total gene families: {len(combined_matrix)}")
    print(f"  Core genes: {(combined_matrix.sum(axis=1) == len(combined_matrix.columns)).sum()}")
EOF
```

### Example 2: Incremental Analysis Pipeline

```python
#!/usr/bin/env python3
"""Incremental pangenome analysis for continuously growing datasets."""

import os
import pandas as pd
from pathlib import Path
from pangenomeplus.pipeline import run_pangenome_pipeline
from pangenomeplus.utils.config import create_default_configuration

def incremental_analysis(genomes_dir, results_dir, checkpoint_dir):
    """Run incremental pangenome analysis."""

    genomes_dir = Path(genomes_dir)
    results_dir = Path(results_dir)
    checkpoint_dir = Path(checkpoint_dir)

    # Create directories
    results_dir.mkdir(exist_ok=True)
    checkpoint_dir.mkdir(exist_ok=True)

    # Check for existing analysis
    checkpoint_file = checkpoint_dir / 'processed_genomes.txt'
    processed_genomes = set()
    if checkpoint_file.exists():
        with open(checkpoint_file) as f:
            processed_genomes = set(line.strip() for line in f)

    # Find new genomes
    all_genomes = set(f.name for f in genomes_dir.glob('*.fasta*'))
    new_genomes = all_genomes - processed_genomes

    if not new_genomes:
        print("No new genomes to process")
        return

    print(f"Processing {len(new_genomes)} new genomes")

    # Create temporary directory with all genomes for analysis
    temp_dir = results_dir / 'temp_analysis'
    temp_dir.mkdir(exist_ok=True)

    # Copy all genomes (old + new) to temp directory
    for genome in all_genomes:
        src = genomes_dir / genome
        dst = temp_dir / genome
        if not dst.exists():
            os.link(src, dst)  # Hard link to save space

    # Configure analysis
    config = create_default_configuration()
    config["resources"]["threads"] = 8
    config["clustering"]["protein"]["identity"] = 0.8

    # Run complete analysis
    results = run_pangenome_pipeline(
        genomes_dir=temp_dir,
        output_dir=results_dir / 'latest',
        config=config
    )

    # Update checkpoint
    with open(checkpoint_file, 'w') as f:
        for genome in all_genomes:
            f.write(f"{genome}\n")

    # Clean up temp directory
    import shutil
    shutil.rmtree(temp_dir)

    print(f"Analysis updated:")
    print(f"  Total genomes: {results['n_genomes']}")
    print(f"  Total families: {results['n_families']}")

    return results

if __name__ == "__main__":
    results = incremental_analysis(
        genomes_dir="./incoming_genomes",
        results_dir="./pangenome_results",
        checkpoint_dir="./checkpoints"
    )
```

## Comparative Studies

### Example 1: Before/After Treatment Comparison

```bash
# Organize samples by treatment
mkdir treatment_study
cd treatment_study
mkdir pre_treatment post_treatment

# Analyze each timepoint
for timepoint in pre_treatment post_treatment; do
    python -m pangenomeplus \
      $timepoint/ \
      results_$timepoint/ \
      --sensitivity 8.0 \
      --protein-identity 0.95 \
      --output-formats transformer,presence_absence
done

# Compare results
python3 << EOF
import pandas as pd

# Load matrices
pre = pd.read_csv('results_pre_treatment/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
post = pd.read_csv('results_post_treatment/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)

# Find treatment-associated changes
# Genes lost after treatment
common_families = set(pre.index) & set(post.index)
lost_families = []
gained_families = []

for family in common_families:
    pre_prevalence = pre.loc[family].mean()
    post_prevalence = post.loc[family].mean()

    if pre_prevalence > 0.8 and post_prevalence < 0.2:
        lost_families.append((family, pre_prevalence, post_prevalence))
    elif pre_prevalence < 0.2 and post_prevalence > 0.8:
        gained_families.append((family, pre_prevalence, post_prevalence))

print(f"Treatment effect analysis:")
print(f"  Families lost: {len(lost_families)}")
print(f"  Families gained: {len(gained_families)}")

# Show top changes
print("\nTop families lost after treatment:")
for family, pre_prev, post_prev in lost_families[:5]:
    print(f"  {family}: {pre_prev:.2f} → {post_prev:.2f}")

print("\nTop families gained after treatment:")
for family, pre_prev, post_prev in gained_families[:5]:
    print(f"  {family}: {pre_prev:.2f} → {post_prev:.2f}")
EOF
```

### Example 2: Geographic Comparison

```bash
# Multi-location comparison
mkdir geographic_study
cd geographic_study

# Organize by location
mkdir north_america europe asia

# Analyze each region
for region in north_america europe asia; do
    echo "Analyzing $region..."
    python -m pangenomeplus \
      $region/ \
      results_$region/ \
      --sensitivity 7.0 \
      --protein-identity 0.85
done

# Cross-regional analysis
mkdir all_regions
cp north_america/*.fasta all_regions/
cp europe/*.fasta all_regions/
cp asia/*.fasta all_regions/

python -m pangenomeplus all_regions/ global_results/ \
  --sensitivity 7.0 \
  --protein-identity 0.85

# Analyze geographic patterns
python3 << EOF
import pandas as pd
import numpy as np

# Load regional results
na_matrix = pd.read_csv('results_north_america/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
eu_matrix = pd.read_csv('results_europe/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
as_matrix = pd.read_csv('results_asia/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)
global_matrix = pd.read_csv('global_results/outputs/presence_absence_matrix.tsv', sep='\t', index_col=0)

# Identify region-specific genes
regions = {
    'North America': [col for col in global_matrix.columns if 'NA_' in col],
    'Europe': [col for col in global_matrix.columns if 'EU_' in col],
    'Asia': [col for col in global_matrix.columns if 'AS_' in col]
}

print("Regional pangenome comparison:")
for region, genomes in regions.items():
    if genomes:
        region_prevalence = global_matrix[genomes].mean(axis=1)
        other_genomes = [g for g in global_matrix.columns if g not in genomes]
        other_prevalence = global_matrix[other_genomes].mean(axis=1)

        # Find region-specific families
        region_specific = global_matrix[
            (region_prevalence > 0.8) & (other_prevalence < 0.1)
        ].index

        print(f"  {region}: {len(genomes)} genomes, {len(region_specific)} region-specific families")
EOF
```

## Troubleshooting Examples

### Memory Issues

```bash
# Diagnose memory usage
echo "System memory information:"
free -h

# Run with memory monitoring
python -m pangenomeplus genomes/ results/ \
  --memory-limit 16G \
  --threads 4 \
  --sensitivity 5.0 \
  --log-level DEBUG 2>&1 | tee analysis.log

# If still failing, use minimal mode
python -m pangenomeplus genomes/ results_minimal/ \
  --annotation-tools prodigal \
  --sensitivity 4.0 \
  --protein-identity 0.7 \
  --threads 2 \
  --output-formats presence_absence
```

### Large Dataset Optimization

```bash
# Progressive optimization for large datasets
echo "Dataset size analysis:"
find genomes/ -name "*.fasta*" | wc -l

# Test with subset first
mkdir test_subset
cp genomes/*.fasta test_subset/ | head -10  # First 10 genomes

python -m pangenomeplus test_subset/ test_results/ \
  --sensitivity 6.0 \
  --log-level DEBUG

# If successful, scale up gradually
mkdir medium_subset
cp genomes/*.fasta medium_subset/ | head -50  # 50 genomes

python -m pangenomeplus medium_subset/ medium_results/ \
  --sensitivity 6.0 \
  --use-gpu \
  --threads 8

# Full analysis with optimized parameters
python -m pangenomeplus genomes/ full_results/ \
  --sensitivity 5.5 \
  --use-gpu \
  --threads 16 \
  --memory-limit 64G
```

### Tool Compatibility Issues

```bash
# Test external tool versions
echo "Checking tool compatibility:"

# MMseqs2 version check
mmseqs version
if [ $? -ne 0 ]; then
    echo "MMseqs2 not found or incompatible"
    exit 1
fi

# Prodigal check
prodigal -v
if [ $? -ne 0 ]; then
    echo "Prodigal not found or incompatible"
    exit 1
fi

# Test with minimal analysis
python -m pangenomeplus \
  --analyze-only \
  genomes/ \
  --log-level DEBUG

# Run with tool debugging
python -m pangenomeplus genomes/ debug_results/ \
  --log-level DEBUG \
  --threads 1 \
  --sensitivity 4.0 2>&1 | tee debug.log

# Check for specific error patterns
grep -i "error\|failed\|exception" debug.log
```

## Advanced Workflows

### Example 1: Quality Control Pipeline

```python
#!/usr/bin/env python3
"""Quality control pipeline for pangenome analysis."""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def quality_control_analysis(results_dir):
    """Perform quality control on pangenome results."""

    results_dir = Path(results_dir)

    # Load results
    matrix = pd.read_csv(results_dir / 'outputs' / 'presence_absence_matrix.tsv',
                        sep='\t', index_col=0)
    family_summary = pd.read_csv(results_dir / 'families' / 'family_summary.tsv',
                                sep='\t')

    # Quality metrics
    n_genomes = len(matrix.columns)
    n_families = len(matrix)

    # Genome-level statistics
    genome_stats = pd.DataFrame({
        'genome': matrix.columns,
        'n_families': matrix.sum(axis=0),
        'unique_families': (matrix == 1).sum(axis=0)
    })

    # Family-level statistics
    family_stats = pd.DataFrame({
        'family': matrix.index,
        'prevalence': matrix.sum(axis=1),
        'prevalence_pct': matrix.sum(axis=1) / n_genomes * 100
    })

    # Quality checks
    print("Quality Control Report")
    print("=" * 50)

    # Check for outlier genomes
    mean_families = genome_stats['n_families'].mean()
    std_families = genome_stats['n_families'].std()
    outliers = genome_stats[
        (genome_stats['n_families'] < mean_families - 2*std_families) |
        (genome_stats['n_families'] > mean_families + 2*std_families)
    ]

    print(f"Dataset overview:")
    print(f"  Genomes: {n_genomes}")
    print(f"  Gene families: {n_families}")
    print(f"  Average families per genome: {mean_families:.1f} ± {std_families:.1f}")

    if len(outliers) > 0:
        print(f"\nPotential outlier genomes ({len(outliers)}):")
        for _, row in outliers.iterrows():
            print(f"  {row['genome']}: {row['n_families']} families")

    # Core genome analysis
    core_threshold = 0.95
    core_families = family_stats[family_stats['prevalence_pct'] >= core_threshold*100]
    accessory_families = family_stats[family_stats['prevalence_pct'] < core_threshold*100]
    singleton_families = family_stats[family_stats['prevalence'] == 1]

    print(f"\nPangenome composition:")
    print(f"  Core genes (≥{core_threshold*100}%): {len(core_families)} ({len(core_families)/n_families*100:.1f}%)")
    print(f"  Accessory genes: {len(accessory_families)} ({len(accessory_families)/n_families*100:.1f}%)")
    print(f"  Singleton genes: {len(singleton_families)} ({len(singleton_families)/n_families*100:.1f}%)")

    # Generate plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Families per genome
    axes[0,0].hist(genome_stats['n_families'], bins=20, alpha=0.7, edgecolor='black')
    axes[0,0].axvline(mean_families, color='red', linestyle='--', label=f'Mean: {mean_families:.0f}')
    axes[0,0].set_xlabel('Families per genome')
    axes[0,0].set_ylabel('Number of genomes')
    axes[0,0].set_title('Distribution of families per genome')
    axes[0,0].legend()

    # Family prevalence
    axes[0,1].hist(family_stats['prevalence_pct'], bins=20, alpha=0.7, edgecolor='black')
    axes[0,1].axvline(core_threshold*100, color='red', linestyle='--', label=f'Core threshold: {core_threshold*100}%')
    axes[0,1].set_xlabel('Family prevalence (%)')
    axes[0,1].set_ylabel('Number of families')
    axes[0,1].set_title('Distribution of family prevalence')
    axes[0,1].legend()

    # Unique families per genome
    axes[1,0].scatter(genome_stats['n_families'], genome_stats['unique_families'], alpha=0.6)
    axes[1,0].set_xlabel('Total families')
    axes[1,0].set_ylabel('Unique families')
    axes[1,0].set_title('Unique vs total families per genome')

    # Cumulative pangenome curve
    prevalence_counts = family_stats['prevalence'].value_counts().sort_index()
    cumulative = prevalence_counts.cumsum()
    axes[1,1].plot(prevalence_counts.index, cumulative.values, 'o-')
    axes[1,1].set_xlabel('Minimum prevalence')
    axes[1,1].set_ylabel('Cumulative families')
    axes[1,1].set_title('Cumulative pangenome curve')

    plt.tight_layout()
    plt.savefig(results_dir / 'quality_control_plots.png', dpi=300, bbox_inches='tight')

    # Save QC report
    qc_report = {
        'dataset_stats': {
            'n_genomes': n_genomes,
            'n_families': n_families,
            'mean_families_per_genome': float(mean_families),
            'std_families_per_genome': float(std_families)
        },
        'pangenome_composition': {
            'core_families': len(core_families),
            'accessory_families': len(accessory_families),
            'singleton_families': len(singleton_families),
            'core_percentage': float(len(core_families)/n_families*100)
        },
        'outlier_genomes': outliers['genome'].tolist()
    }

    import json
    with open(results_dir / 'quality_control_report.json', 'w') as f:
        json.dump(qc_report, f, indent=2)

    print(f"\nQuality control complete. Report saved to {results_dir}/quality_control_report.json")

    return qc_report

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python qc_pipeline.py <results_directory>")
        sys.exit(1)

    quality_control_analysis(sys.argv[1])
```

### Example 2: Automated Report Generation

```python
#!/usr/bin/env python3
"""Automated pangenome analysis report generator."""

import pandas as pd
import json
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

def generate_analysis_report(results_dir, output_file="pangenome_report.html"):
    """Generate comprehensive HTML report from PanGenomePlus results."""

    results_dir = Path(results_dir)

    # Load data
    matrix = pd.read_csv(results_dir / 'outputs' / 'presence_absence_matrix.tsv',
                        sep='\t', index_col=0)

    # Calculate statistics
    n_genomes = len(matrix.columns)
    n_families = len(matrix)
    core_genes = (matrix.sum(axis=1) == n_genomes).sum()
    accessory_genes = n_families - core_genes
    singletons = (matrix.sum(axis=1) == 1).sum()

    # Generate plots
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Plot 1: Pangenome composition pie chart
    sizes = [core_genes, accessory_genes - singletons, singletons]
    labels = ['Core', 'Accessory', 'Singletons']
    colors = ['#ff9999', '#66b3ff', '#99ff99']

    axes[0,0].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[0,0].set_title('Pangenome Composition')

    # Plot 2: Family prevalence distribution
    prevalence = matrix.sum(axis=1) / n_genomes * 100
    axes[0,1].hist(prevalence, bins=20, alpha=0.7, edgecolor='black')
    axes[0,1].set_xlabel('Prevalence (%)')
    axes[0,1].set_ylabel('Number of families')
    axes[0,1].set_title('Family Prevalence Distribution')

    # Plot 3: Families per genome
    families_per_genome = matrix.sum(axis=0)
    axes[1,0].hist(families_per_genome, bins=15, alpha=0.7, edgecolor='black')
    axes[1,0].set_xlabel('Families per genome')
    axes[1,0].set_ylabel('Number of genomes')
    axes[1,0].set_title('Families per Genome Distribution')

    # Plot 4: Heatmap of core genes (if manageable size)
    if core_genes <= 50:  # Only show if reasonable number
        core_matrix = matrix[matrix.sum(axis=1) == n_genomes]
        sns.heatmap(core_matrix, cmap='Blues', cbar=False, ax=axes[1,1])
        axes[1,1].set_title(f'Core Gene Presence ({core_genes} genes)')
    else:
        axes[1,1].text(0.5, 0.5, f'Core genes: {core_genes}\n(Too many to display)',
                      ha='center', va='center', transform=axes[1,1].transAxes)
        axes[1,1].set_title('Core Genes Summary')

    plt.tight_layout()
    plt.savefig(results_dir / 'analysis_plots.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Generate HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PanGenomePlus Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
            .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
            .stat-box {{ background-color: #e8f4fd; padding: 15px; border-radius: 5px; text-align: center; }}
            .stat-number {{ font-size: 24px; font-weight: bold; color: #2c5aa0; }}
            .stat-label {{ font-size: 14px; color: #666; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .plot {{ text-align: center; margin: 20px 0; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>PanGenomePlus Analysis Report</h1>
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Dataset:</strong> {results_dir.name}</p>
        </div>

        <h2>Summary Statistics</h2>
        <div class="stats">
            <div class="stat-box">
                <div class="stat-number">{n_genomes}</div>
                <div class="stat-label">Genomes</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{n_families:,}</div>
                <div class="stat-label">Gene Families</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{core_genes:,}</div>
                <div class="stat-label">Core Genes</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{core_genes/n_families*100:.1f}%</div>
                <div class="stat-label">Core Percentage</div>
            </div>
        </div>

        <h2>Pangenome Composition</h2>
        <table>
            <tr><th>Category</th><th>Count</th><th>Percentage</th></tr>
            <tr><td>Core genes (100%)</td><td>{core_genes:,}</td><td>{core_genes/n_families*100:.2f}%</td></tr>
            <tr><td>Accessory genes</td><td>{accessory_genes-singletons:,}</td><td>{(accessory_genes-singletons)/n_families*100:.2f}%</td></tr>
            <tr><td>Singleton genes</td><td>{singletons:,}</td><td>{singletons/n_families*100:.2f}%</td></tr>
            <tr><td><strong>Total</strong></td><td><strong>{n_families:,}</strong></td><td><strong>100.00%</strong></td></tr>
        </table>

        <h2>Analysis Plots</h2>
        <div class="plot">
            <img src="analysis_plots.png" alt="Analysis Plots" style="max-width: 100%;">
        </div>

        <h2>Genome Statistics</h2>
        <table>
            <tr><th>Genome</th><th>Gene Families</th><th>Unique Families</th></tr>
    """

    # Add genome statistics
    for genome in matrix.columns:
        total_families = matrix[genome].sum()
        unique_families = ((matrix[genome] == 1) & (matrix.sum(axis=1) == 1)).sum()
        html_content += f"<tr><td>{genome}</td><td>{total_families}</td><td>{unique_families}</td></tr>"

    html_content += """
        </table>

        <h2>Files Generated</h2>
        <ul>
            <li><strong>presence_absence_matrix.tsv:</strong> Binary presence/absence matrix</li>
            <li><strong>transformer_format.txt:</strong> AI-ready sequence format</li>
            <li><strong>gene_to_family.tsv:</strong> Gene-to-family mappings</li>
            <li><strong>family_summary.tsv:</strong> Family statistics</li>
            <li><strong>pangenome_stats.json:</strong> Analysis metadata</li>
        </ul>

        <hr>
        <p><em>Report generated by PanGenomePlus automated analysis pipeline</em></p>
    </body>
    </html>
    """

    # Save report
    with open(results_dir / output_file, 'w') as f:
        f.write(html_content)

    print(f"Report generated: {results_dir / output_file}")
    return results_dir / output_file

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python generate_report.py <results_directory>")
        sys.exit(1)

    generate_analysis_report(sys.argv[1])
```

---

These examples demonstrate the flexibility and power of PanGenomePlus across various research scenarios. Each example can be adapted to specific research needs by adjusting parameters and analysis approaches.

For additional examples and use cases, check the [TUTORIAL.md](TUTORIAL.md) for step-by-step guidance or the [API.md](API.md) for programmatic usage.