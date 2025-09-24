# PanGenomePlus Project Objectives

## Primary Objective
Create a comprehensive pangenome analysis pipeline, in the python programming language,  that analyzes ALL genomic features across multiple genomes, including:

- Protein-coding genes
- Intergenic regions
- tRNAs
- rRNAs
- CRISPR arrays
- Other non-coding regions

## Core Philosophy

### KISS Principle (Keep It Simple, Stupid)
- Heavy lifting performed by specialized external tools
- Python code orchestrates tools and processes outputs
- No reinventing the wheel for complex bioinformatics tasks

### YAGNI Principle (You Aren't Gonna Need It)
- Implement only what is needed for core functionality
- Avoid over-engineering and premature optimization

### SOLID Principles
- Single Responsibility: Each module has one clear purpose
- Open/Closed: Extensible without modification
- Interface Segregation: Clear, minimal interfaces
- Dependency Inversion: Depend on abstractions, not concretions

## External Tool Specifications

### 1. Prodigal (Gene Prediction)
**Purpose**: Predict protein-coding genes in prokaryotic genomes

**Input**:
- FASTA genome files
- Command: `prodigal -i genome.fasta -o genes.gff -a proteins.faa -d nucleotides.fna -f gff`

**Output**:
- GFF3 annotation file with gene coordinates
- Protein sequences (.faa)
- Nucleotide sequences (.fna)

**Processing Requirements**:
- Parse GFF3 to extract gene coordinates (single-pass, O(n) algorithm)
- Extract sequences from external tool outputs and assign compact IDs universally across all feature types with coordinate preservation
- Validate sequence lengths match coordinates
- Handle incomplete genes at contig boundaries

**Linear Time Extraction Algorithm**:
```python
# Optimal O(n) approach - features are already coordinate-sorted
def extract_features_linear(gff_file, fasta_file):
    features = []
    # Single pass through pre-sorted GFF3 - no additional sorting needed
    for line in gff_file:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        start, end = int(fields[3]), int(fields[4])
        features.append(GenomeFeature(start=start, end=end, ...))

    # Direct sequence extraction using coordinates
    genome_seq = SeqIO.read(fasta_file, "fasta")
    for feature in features:
        feature.sequence = str(genome_seq.seq[feature.start-1:feature.end])

    return features  # Total time: O(n) where n = number of features
```

**Performance Targets**:
- Processing rate: >1,000 features/second per genome
- Memory usage: <100MB per genome during extraction
- Linear scaling: 200 genomes should take ~200x single genome time

### 2. MMseqs2 (Protein Clustering)
**Purpose**: Cluster protein sequences by homology to identify gene families

**Input**:
- Protein FASTA file (all genomes combined)
- Parameters: --min-seq-id 0.8 --coverage 0.8 --cov-mode 1

**Output**:
- Cluster assignments (representative -> members)
- Cluster representatives FASTA
- Alignment information

**Processing Requirements**:
- Parse cluster TSV format (representative<tab>member)
- Map cluster assignments back to original gene IDs
- Create gene family mappings
- Generate family size statistics
- Classify families as core/accessory/cloud

### 3. tRNAscan-SE (tRNA Detection)
**Purpose**: Identify and annotate tRNA genes

**Input**:
- FASTA genome files
- Command: `tRNAscan-SE -B genome.fasta`

**Output**:
- Tab-delimited file with tRNA coordinates and anticodons
- Predicted secondary structures

**Processing Requirements**:
- Parse coordinate format (genome, start, end, anticodon)
- Extract tRNA sequences from genome
- Convert to GFF3 format for consistency
- Handle overlapping predictions

### 4. Barrnap (rRNA Detection)
**Purpose**: Identify ribosomal RNA genes

**Input**:
- FASTA genome files
- Command: `barrnap --kingdom bac genome.fasta`

**Output**:
- GFF3 format with rRNA annotations
- Subunit classifications (16S, 23S, 5S)

**Processing Requirements**:
- Parse GFF3 rRNA annotations
- Extract rRNA sequences
- Classify by subunit type
- Validate length ranges for each subunit

### 5. MINCED (CRISPR Detection)
**Purpose**: Identify CRISPR arrays and spacers

**Input**:
- FASTA genome files
- Command: `minced genome.fasta output.txt`

**Output**:
- Text format with CRISPR coordinates
- Spacer sequences and repeats

**Processing Requirements**:
- Parse MINCED text output format
- Extract spacer sequences
- Convert to GFF3 annotations
- Handle variable spacer lengths

## Complete Pipeline Flow

### Overview
**Raw Genomes** → **Annotation Tools** → **Feature Extraction** → **Compact ID Assignment** → **Clustering** → **Family Assignment** → **Core/Accessory/Cloud Classification** → **Multi-Format Output Generation**

### Detailed Pipeline Stages

### 1. Annotation Stage
**Inputs**: Genome FASTA files
**Processes** (parallel execution):
- Run Prodigal for protein-coding genes → protein sequences (.faa) + coordinates (.gff)
- Run tRNAscan-SE for tRNAs → tRNA coordinates (.tab)
- Run Barrnap for rRNAs → rRNA annotations (.gff)
- Run MINCED for CRISPR elements → CRISPR spacers (.txt)
**Outputs**: Native output formats from each external tool

### 2. Feature Extraction & Compact ID Assignment Stage
**Inputs**: Native tool outputs + genome FASTA
**Processes**:
- Parse each tool's native output format
- Extract sequences using genomic coordinates
- Generate intergenic regions (pipeline-calculated between features)
- **Universal compact ID assignment across ALL feature types**:
  - Proteins: P1, P2, P3...
  - tRNAs: T1, T2, T3...
  - rRNAs: R1, R2, R3...
  - CRISPR: C1, C2, C3...
  - Intergenic: I1, I2, I3...
- **Create mapping tables**: compact_id ↔ full_metadata (including coordinates)
**Outputs**:
- Feature-type FASTA files with compact headers:
  - `protein_coding_gene.fasta` (>P1, >P2, >P3...)
  - `tRNA.fasta` (>T1, >T2, >T3...)
  - `rRNA.fasta` (>R1, >R2, >R3...)
  - `CRISPR.fasta` (>C1, >C2, >C3...)
  - `intergenic_region.fasta` (>I1, >I2, >I3...)
- Mapping tables with coordinate preservation

### 3. Clustering Stage
**Inputs**: Feature-type FASTA files with compact headers
**Processes**:
- **Protein sequences**: MMseqs2 clustering (unidirectional - trust Prodigal orientations)
- **rRNA sequences**: MMseqs2 clustering (unidirectional - trust Barrnap orientations)
- **tRNA sequences**: MMseqs2 clustering (unidirectional - trust tRNAscan-SE orientations)
- **Bidirectional nucleotide clustering**: MMseqs2 clustering with bidirectional search (combined FASTA approach):
  - Process: Create combined FASTA with forward + reverse-complement sequences
  - Parameters: --search-type 1 (nucleotide mode)
  - Post-processing: Map results back to original sequence orientations
  - Feature types:
    - Intergenic regions (relaxed parameters: 70% identity, 50% coverage)
    - CRISPR spacers (horizontal transfer affects orientation)
**Rationale**: Trust specialized tools' strand predictions; use bidirectional only for features without inherent orientation
**Outputs**: Cluster assignment files per feature type using compact IDs

### 4. Family Assignment Stage
**Inputs**: Cluster assignments + compact ID mappings
**Processes**:
- **Create family IDs with feature type prefixes**:
  - Family IDs: FAM_P1, FAM_T42, FAM_I5K, FAM_R7, FAM_C12
  - Singletons: SING_P2B7K, SING_T15 (tracked but not counted as families)
- Only clustered groups count as true gene families
- Calculate family size distributions
- Classify families by genome presence (core/accessory/cloud)
**Outputs**:
- `gene_to_family.tsv` (using compact IDs)
- `family_summary.tsv`

### 5. Pangenome Classification Stage
**Inputs**: Family assignments + genome counts
**Processes**:
- **Core families**: ≥95% genomes
- **Accessory families**: 15-95% genomes
- **Cloud families**: <15% genomes
- Classification applied separately per feature type
**Outputs**: Classified family lists per feature type

### 6. Multi-Format Output Generation Stage
**Inputs**: Family assignments + mapping tables + coordinate information
**Processes**:
- **Coordinate-based reconstruction**: Use mapping tables to recover genomic positions
- **Transformer format generation**:
  - Format: `genome_name P1 I5 P2 T42 P3 I7 R2 P4`
  - **Critical requirement**: Families ordered by genomic coordinate (NOT alphabetical)
  - Process: Retrieve families present in genome → Sort by start coordinate → Output space-separated tokens
- **Traditional format generation**:
  - Presence/absence matrices with meaningful biological headers
  - Roary-compatible output with original gene names
  - Summary statistics with biological context
**Outputs**:
- **Transformer format**: Coordinate-ordered family sequences per genome
- **Presence/absence matrices**: Binary matrices (CSV, TSV)
- **Roary-compatible format**: Traditional pangenome analysis format
- **Summary statistics**: Family counts, core/accessory breakdown per feature type
- **Visualizations**: Plots and pangenome charts

### Critical Technical Requirements

#### MMseqs2 Search Directionality
- **Proteins**: Unidirectional - trust Prodigal strand predictions
- **rRNAs**: Unidirectional - trust Barrnap strand predictions
- **tRNAs**: Unidirectional - trust tRNAscan-SE strand predictions
- **Intergenic regions**: Bidirectional - no inherent orientation
- **CRISPR spacers**: Bidirectional - horizontal transfer affects orientation

#### Coordinate Preservation
- All compact IDs linked to genomic coordinates in mapping tables
- Transformer output must preserve exact genomic order of families
- No alphabetical or frequency-based sorting in final output

## Feature Type Specifications

### Protein-Coding Genes
- **Clustering Method**: MMseqs2 (80% identity, 80% coverage)
- **Search Direction**: Forward strand only (pre-oriented by Prodigal)
- **Family Definition**: Homologous groups from clustering
- **Expected Range**: 15,000-70,000 families for E. coli
- **Extraction**: Direct from Prodigal .faa files using sequence headers
- **Compact ID Assignment**: P1, P2, P3, ... (sequential assignment)
- **Performance**: Should process >5,000 proteins/second

### Intergenic Regions
- **Definition**: Non-annotated regions ≥50bp between features
- **Clustering Method**: MMseqs2 (relaxed: 70% identity, 50% coverage)
- **Search Direction**: Bidirectional using combined FASTA (forward + reverse-complement sequences)
- **Rationale**: Recombination events affect orientation; no inherent directionality
- **Special Handling**: High variability expected, many singletons normal
- **Compact ID Assignment**: I1, I2, I3, ... (sequential assignment)
- **Extraction Algorithm**:
  ```python
  # Linear time intergenic extraction - O(n) where n = number of features
  def extract_intergenic_linear(features, genome_seq, min_length=50):
      intergenic = []
      prev_end = 1
      for feature in sorted_features:  # Already sorted from GFF3
          if feature.start - prev_end >= min_length:
              intergenic.append(GenomeFeature(
                  start=prev_end,
                  end=feature.start-1,
                  sequence=genome_seq[prev_end-1:feature.start-1]
              ))
          prev_end = max(prev_end, feature.end + 1)
      return intergenic
  ```

### tRNAs
- **Clustering Method**: MMseqs2 sequence similarity clustering
- **Search Direction**: Unidirectional (trust tRNAscan-SE strand predictions)
- **Rationale**: tRNAscan-SE provides accurate strand orientations using specialized models
- **Expected Count**: ~40-80 tRNAs per genome (20 amino acids × 1-4 isoacceptors)
- **Validation**: Check anticodon predictions and secondary structure
- **Compact ID Assignment**: T1, T2, T3, ... (sequential assignment)
- **Extraction**: Parse tRNAscan-SE output coordinates and extract sequences with strand information

### rRNAs
- **Clustering Method**: MMseqs2 clustering (separate by subunit type: 16S, 23S, 5S)
- **Search Direction**: Unidirectional (trust Barrnap strand predictions)
- **Rationale**: Barrnap provides accurate strand orientations using HMM-based detection
- **Expected Count**: 1-7 rRNA operons per genome
- **Validation**: Length ranges (16S: ~1,500bp, 23S: ~2,900bp, 5S: ~120bp)
- **Compact ID Assignment**: R1, R2, R3, ... (sequential assignment)
- **Extraction**: Direct from Barrnap GFF3 output coordinates with strand information

### CRISPR Elements
- **Clustering Method**: MMseqs2 spacer sequence similarity
- **Search Direction**: Bidirectional using combined FASTA (forward + reverse-complement sequences)
- **Rationale**: CRISPR spacers can have variable orientations; horizontal transfer
- **Variable Presence**: Not all genomes contain CRISPR systems
- **Special Handling**: Account for horizontal transfer of spacers
- **Compact ID Assignment**: C1, C2, C3, ... (sequential assignment)
- **Extraction**: Parse MINCED text output format for spacer coordinates

## Quality Control Requirements

### Input Validation
- Verify genome file integrity (FASTA format, sequence lengths)
- Check annotation completeness (minimum gene count thresholds)
- Validate external tool outputs before processing

### Process Monitoring
- Track clustering statistics (success rates, singleton counts)
- Monitor memory usage and processing times
- Log warnings for unusual patterns (too many/few features)

### Output Validation
- Verify family count ranges against literature expectations
- Check core genome size consistency
- Validate presence/absence matrix completeness

## Error Handling Strategy

### Tool Failures
- Retry failed tool executions with adjusted parameters
- Graceful degradation (continue without failed feature types)
- Clear error reporting with suggested fixes

### Data Quality Issues
- Skip problematic genomes with warnings
- Report statistics on filtered/excluded data
- Provide recommendations for data quality improvement

### Performance Issues
- Implement timeouts for long-running operations
- Add progress monitoring for large datasets
- Optimize algorithms before adding complexity

## Performance Anti-Patterns to Avoid

### Feature Extraction Anti-Patterns
1. **Redundant Sorting**: GFF3 files from Prodigal are already coordinate-sorted
   ```python
   # WRONG - O(n log n) unnecessary sorting
   features = sorted(features, key=lambda x: x.start)

   # RIGHT - O(n) direct processing, features already sorted
   for feature in features:  # Process in order
   ```

2. **Repeated File I/O**: Load sequences once, reuse in memory
   ```python
   # WRONG - O(n) file reads for n features
   for feature in features:
       seq = SeqIO.read(genome_file, "fasta")
       extract_sequence(seq, feature)

   # RIGHT - O(1) file read, O(n) sequence extraction
   genome_seq = SeqIO.read(genome_file, "fasta")
   for feature in features:
       extract_sequence(genome_seq, feature)
   ```

3. **Complex Data Structures**: Use simple lists when possible
   ```python
   # WRONG - unnecessary overhead for simple operations
   feature_dict = {f.id: f for f in features}
   sorted_features = sorted(feature_dict.values())

   # RIGHT - direct list processing
   for feature in features:  # Already sorted from GFF3
   ```

### Expected Linear Time Operations
- **GFF3 Parsing**: O(n) where n = lines in file
- **Sequence Extraction**: O(m) where m = total sequence length
- **Intergenic Detection**: O(f) where f = number of features
- **Total Extraction**: O(n + m + f) = linear in input size

### Memory Efficiency Guidelines
- Process one genome at a time to limit memory usage
- Use generators/iterators for large datasets
- Release sequence objects after processing each genome
- Avoid storing full dataset in memory simultaneously

## Ultra-Compact Sequence ID System

### Purpose and Scale Requirements
Design ultra-compact sequence identifiers optimized for:
- **Scale**: 9,999,999+ genomes
- **Sequence count**: Potentially 1+ trillion sequences across all feature types
- **MMseqs2 optimization**: Minimal memory usage and optimal clustering performance
- **Storage efficiency**: 95% reduction compared to current verbose IDs

### Current Problem
Current IDs are extremely verbose and inefficient:
```
PGP_E_coli_001_GCF_000597845_protein_coding_gene_NZ_CP007265.1_1_NZ_CP007265.1_47_1450
```
- **Length**: ~90+ characters
- **Memory impact**: Massive overhead for billion+ sequences
- **MMseqs2 performance**: Suboptimal due to long headers

### Compact ID Format Specification

#### Core Format: `{TYPE}{BASE36_ID}`
- **Feature type prefix** (1 character):
  - `P`: Protein-coding genes
  - `I`: Intergenic regions
  - `T`: tRNAs
  - `R`: rRNAs
  - `C`: CRISPR elements

- **Base36 sequence ID** (variable length):
  - **Character set**: 0-9, A-Z (36 total characters)
  - **Sequential assignment**: IDs assigned incrementally during processing
  - **Examples**: `1`, `2B7K`, `9ZZZ`, `1A2B3C`

#### Format Examples
```
P1          → First protein sequence
P2B7K       → Protein sequence #2,847,432
I5          → Fifth intergenic region
T42         → 42nd tRNA sequence
R7          → Seventh rRNA sequence
C1A         → CRISPR sequence #1,296
```

#### Scalability Capacity
- **6 characters**: 36^6 = 2.1 billion unique IDs per feature type
- **7 characters**: 36^7 = 78.3 billion unique IDs per feature type
- **8 characters**: 36^8 = 2.8 trillion unique IDs per feature type

### Mapping System Design

#### Forward Lookup: Compact ID → Full Metadata
**Purpose**: Retrieve complete information from compact ID
**Performance**: O(1) constant time lookup
**Data structure**:
```
compact_to_full = {
    "P1": {
        "genome_id": "E_coli_001",
        "original_id": "gene_1",
        "contig": "NZ_CP007265.1",
        "start": 47,
        "end": 1450,
        "strand": "+",
        "product": "DNA polymerase III subunit alpha"
    },
    "I5": {
        "genome_id": "E_coli_001",
        "original_id": "intergenic_1",
        "contig": "NZ_CP007265.1",
        "start": 2556,
        "end": 2800,
        "strand": ".",
        "product": "intergenic_region"
    }
}
```

#### Reverse Lookup: Genome/Position → Compact ID
**Purpose**: Find compact ID from genomic coordinates
**Performance**: O(1) constant time lookup
**Data structure**:
```
location_to_compact = {
    ("E_coli_001", 47, 1450): "P1",
    ("E_coli_001", 2556, 2800): "I5",
    ("E_coli_002", 100, 500): "P25A"
}
```

#### Genome-Based Lookup: Genome → Feature Lists
**Purpose**: Retrieve all features for a specific genome
**Performance**: O(1) to get list, O(k) to iterate k features
**Data structure**:
```
genome_features = {
    "E_coli_001": {
        "proteins": ["P1", "P2", "P3", "P7B"],
        "intergenic": ["I5", "I6", "I12"],
        "tRNAs": ["T42", "T43"],
        "rRNAs": ["R7"],
        "CRISPR": []
    }
}
```

### Storage and Performance Strategy

#### Small to Medium Scale (up to 1M genomes)
- **Implementation**: In-memory hash tables/dictionaries
- **Memory usage**: ~100MB to 100GB depending on genome count
- **Lookup performance**:
  - Forward: <1 microsecond
  - Reverse: <10 microseconds
- **Persistence**: JSON/Pickle files for mapping tables

#### Large Scale (1M to 10M genomes)
- **Implementation**: Memory-mapped files or embedded database
- **Memory usage**: ~100GB with intelligent caching
- **Lookup performance**:
  - Forward: <10 microseconds
  - Reverse: <100 microseconds
- **Persistence**: SQLite with proper indexing

#### Ultra Scale (10M+ genomes)
- **Implementation**: Full database backend (PostgreSQL/SQLite)
- **Storage**: Terabytes, distributed if necessary
- **Lookup performance**:
  - Forward: <1 millisecond
  - Reverse: <5 milliseconds
- **Database schema**:
```sql
CREATE TABLE sequences (
    compact_id TEXT PRIMARY KEY,
    feature_type CHAR(1),
    genome_id TEXT,
    original_id TEXT,
    contig_id TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    strand CHAR(1),
    product TEXT
);
CREATE INDEX idx_genome_pos ON sequences(genome_id, start_pos);
CREATE INDEX idx_feature_type ON sequences(feature_type);
```

### MMseqs2 Integration Requirements

#### FASTA Output Format
- **Sequence headers**: Use compact IDs exclusively
- **Example FASTA**:
```
>P1
MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLD...
>P2
MKFTVEREHLLKPLQQVSGPLGGRPTLPILGNLLLQVADGTLSLTGTDLE...
>I5
ATCGATCGATCGATCGATCGATCGATCG...
```

#### Mapping File Generation
- **Purpose**: Preserve complete traceability alongside compact FASTA
- **Format**: JSON/TSV files containing full mapping tables
- **Co-location**: Mapping files stored with FASTA files for each feature type

#### Clustering Output Processing
- **MMseqs2 input**: Compact ID FASTA files
- **MMseqs2 output**: Cluster assignments using compact IDs
- **Post-processing**: Map cluster results back to full metadata using lookup tables

### Quality Assurance Requirements

#### ID Uniqueness Validation
- **Global uniqueness**: No duplicate compact IDs across all feature types
- **Sequence counters**: Maintain separate counters per feature type
- **Collision detection**: Verify uniqueness during ID assignment

#### Mapping Integrity Validation
- **Bidirectional consistency**: Forward and reverse lookups must be consistent
- **Completeness**: Every sequence must have complete mapping entry
- **Coordinate validation**: Genomic coordinates must be valid and non-overlapping where expected

#### Performance Monitoring
- **Lookup timing**: Monitor query performance as scale increases
- **Memory usage**: Track memory consumption of mapping structures
- **Database optimization**: Ensure indexes remain effective at scale

### Benefits and Expected Outcomes

#### Storage Efficiency
- **Header reduction**: 90+ characters → 3-8 characters (95% reduction)
- **Memory savings**: Massive reduction in RAM usage for large-scale analysis
- **Disk space**: Significantly smaller FASTA files and faster I/O

#### MMseqs2 Performance
- **Reduced memory usage**: Less RAM required for sequence databases
- **Faster processing**: Shorter IDs improve clustering speed
- **Better scalability**: MMseqs2 can handle larger datasets efficiently

#### Pipeline Scalability
- **Billion-scale ready**: System designed for 1+ trillion sequences
- **Future-proof**: Base36 encoding provides massive ID space
- **Maintenance efficiency**: Simple format easy to debug and validate

## Success Metrics

### Biological Accuracy
- Core gene families match literature expectations (2,000-4,000 for E. coli)
- Total protein families in reasonable range (15,000-70,000)
- Appropriate singleton rates (<20% for proteins)

### Performance Targets
- Process 200 genomes in <4 hours on standard hardware
- Memory usage <32GB for typical datasets
- Linear scaling with genome count

### Output Quality
- Complete presence/absence matrices (no missing values)
- Consistent family naming across output formats
- Comprehensive summary statistics and visualizations

## Software Development Methodology

### Core Development Philosophy
PanGenomePlus development follows collaborative, iterative principles with unwavering commitment to clean, maintainable code that serves current needs without premature optimization.

### Requirement Validation Protocol
Before implementing any solution, automatically:

#### Identify Core Requirements
- **Functionality**: What specific biological analysis capability is needed?
- **Use cases**: How will researchers actually use this feature?
- **Constraints**: What are the essential technical and biological limitations?

#### Question Ambiguous Specifications
- **Unclear requirements**: Seek clarification before coding
- **Speculative features**: Challenge "nice to have" additions
- **Premature optimization**: Resist complexity without demonstrated need
- **Mixed responsibilities**: Ensure single-purpose components

### SOLID Principles Enforcement

#### Single Responsibility Principle
- Each pipeline stage handles exactly one biological or computational concern
- External tools (Prodigal, MMseqs2) remain responsible for their specialized tasks
- Python code orchestrates, does not replicate optimized external functionality

#### Open/Closed Principle
- Pipeline components extensible for new feature types without modifying existing code
- New output formats addable without changing core clustering logic
- External tool integration expandable without touching extraction algorithms

#### Liskov Substitution Principle
- Feature extractors interchangeable regardless of source tool
- Output formatters substitutable without breaking downstream analysis
- Clustering components replaceable while maintaining interface contracts

#### Interface Segregation Principle
- Compact ID system provides minimal, specific interfaces
- External tool wrappers expose only necessary functionality
- Output generators depend only on required data structures

#### Dependency Inversion Principle
- Core pipeline depends on abstractions, not concrete tool implementations
- External tools referenced through stable interfaces
- Database/storage backends abstracted from business logic

### Solution Generation Standards

#### Code Quality Hierarchy
1. **Clarity over cleverness**: Readable code preferred over sophisticated solutions
2. **Simplicity over flexibility**: Current needs addressed without future speculation
3. **Current needs over future possibilities**: YAGNI principle strictly enforced
4. **Explicit over implicit**: Clear interfaces and error handling required

#### Implementation Validation Checklist
Before presenting any solution, verify:
- **Simplicity**: Is this the simplest possible solution for current requirements?
- **Necessity**: Is every component and abstraction immediately needed?
- **Responsibility**: Are concerns properly separated across components?
- **Extensibility**: Can this be extended without modifying existing code?
- **Dependencies**: Are external dependencies properly abstracted?

### Development Workflow

#### Phase 1: Requirements Clarification
- **Business context**: Why is this analysis capability needed?
- **User scenarios**: How will researchers interact with this feature?
- **Technical constraints**: What are the performance and scalability requirements?
- **Integration needs**: How does this fit with existing pipeline stages?

#### Phase 2: Solution Design
- **Propose simplest viable approach**: Start with minimal implementation
- **Identify challenges**: Surface technical and biological complications early
- **Highlight trade-offs**: Make performance vs. complexity decisions explicit
- **Challenge assumptions**: Question whether requirements are correctly understood

#### Phase 3: Test-Driven Implementation
- **Write failing tests first**: Define expected behavior before implementation
- **Implement minimal code**: Satisfy tests with simplest possible solution
- **Verify correctness**: Ensure biological accuracy and technical performance
- **Refactor for clarity**: Improve code readability without changing behavior

### Forbidden Development Patterns

#### Prohibited Practices
- **"Just in case" features**: No functionality without immediate use case
- **Premature abstractions**: No generic solutions without concrete requirements
- **Mixed responsibilities**: No components handling multiple concerns
- **Future requirements**: No implementation of anticipated but unconfirmed needs
- **Premature optimization**: No performance tuning without demonstrated bottlenecks

#### Code Quality Standards
- **Single responsibility per function/class**: Each unit serves one clear purpose
- **Clear interface boundaries**: Minimal, well-defined component interactions
- **Minimal dependencies**: Reduce coupling between pipeline stages
- **Explicit error handling**: Clear failure modes and recovery strategies
- **Biological accuracy validation**: Verify results match literature expectations

### External Tool Integration Philosophy
Following KISS principles, external tools remain unmodified while Python code provides minimal orchestration:
- **Tool expertise**: Leverage specialized software (Prodigal, MMseqs2) for complex tasks
- **Interface simplicity**: Minimal wrappers around external tool execution
- **Format conversion only**: Python handles only input/output format standardization
- **No reinvention**: Never replicate functionality available in optimized external tools

### Development Communication Protocol

#### Structured Response Format
All development communications must follow this standardized structure:

1. **Requirement Clarification**
   - Restate understanding of biological and technical requirements
   - Identify assumptions and constraints
   - Clarify ambiguous specifications before proceeding

2. **Core Solution Design**
   - Present simplest viable approach for current needs
   - Identify key architectural decisions
   - Highlight integration points with existing pipeline stages

3. **Implementation Details**
   - Specify concrete technical approach
   - Define interfaces and data structures
   - Address error handling and edge cases

4. **Key Design Decisions**
   - Explain trade-offs and alternatives considered
   - Justify complexity when unavoidable
   - Document decisions for future reference

5. **Validation Results**
   - Confirm biological accuracy against literature
   - Verify performance meets scalability requirements
   - Demonstrate adherence to SOLID principles

### Collaborative Execution Mode

#### Team Member Roles and Behaviors

**Team Member**
- Proactively engage in iterative development process
- Take ownership of assigned pipeline components
- Contribute to collective code quality and biological accuracy

**Critical Thinker**
- Challenge assumptions about requirements and implementation approaches
- Suggest simpler alternatives when complexity is unnecessary
- Question whether features address real researcher needs

**Quality Guardian**
- Maintain high standards through test-driven development
- Enforce adherence to KISS, YAGNI, and SOLID principles
- Prevent technical debt accumulation

#### Core Principle Maintenance
Continuously demonstrate commitment to:

**KISS (Keep It Simple, Stupid)**
- Choose readable solutions over clever optimizations
- Minimize cognitive load in code review and maintenance
- Prefer explicit over implicit implementation approaches

**YAGNI (You Aren't Gonna Need It)**
- Implement only features with immediate, confirmed use cases
- Resist adding flexibility for hypothetical future requirements
- Remove unused code and abstractions

**SOLID Principles**
- Single responsibility per component and function
- Open for extension, closed for modification
- Substitutable implementations without breaking contracts
- Segregated interfaces for specific concerns
- Abstract dependencies rather than concrete implementations

**DRY (Don't Repeat Yourself)**
- Extract common functionality into reusable components
- Avoid duplicating business logic across pipeline stages
- Centralize configuration and shared constants

#### Professional Behavior Standards

**Ownership**
- Take full responsibility for code quality and biological accuracy
- Follow through on commitments and deadlines
- Proactively address technical debt and maintainability issues

**Initiative**
- Identify potential issues before they become problems
- Suggest improvements to existing pipeline components
- Contribute solutions beyond immediate assigned tasks

**Collaboration**
- Engage in constructive dialogue during design discussions
- Provide specific, actionable feedback during code review
- Share knowledge and mentor others on team practices

### Error Handling and Correction Protocol

#### Principle Violation Detection
When identifying violations of development principles:

1. **Identify Specific Principle Breach**
   - Name the violated principle (KISS, YAGNI, SOLID, DRY)
   - Point to exact code location or design decision
   - Explain why current approach violates principle

2. **Explain Violation Clearly**
   - Describe negative consequences of current approach
   - Identify maintenance, testing, or complexity issues
   - Reference specific examples where possible

3. **Provide Simplest Correction**
   - Propose minimal change to address violation
   - Avoid introducing additional complexity in fix
   - Maintain existing functionality and interfaces

4. **Verify Correction Maintains Requirements**
   - Confirm biological accuracy unchanged
   - Validate performance characteristics preserved
   - Ensure integration points remain stable

### Continuous Validation Framework

#### Ongoing Monitoring
Throughout all development activities, continuously monitor for:

**Scope Creep**
- Features beyond confirmed researcher requirements
- Functionality for hypothetical future use cases
- Complexity not justified by immediate needs

**Unnecessary Complexity**
- Abstractions without concrete implementation variations
- Generalized solutions for specific problems
- Over-engineered interfaces and data structures

**Mixed Responsibilities**
- Components handling multiple unrelated concerns
- Business logic mixed with I/O operations
- External tool integration coupled with data processing

**Premature Optimization**
- Performance tuning without demonstrated bottlenecks
- Complex algorithms for small-scale data
- Memory optimizations before memory pressure identified

#### Correction Strategies
When violations detected, correct by:

**Returning to Core Requirements**
- Review original biological analysis objectives
- Eliminate functionality not serving researcher needs
- Focus implementation on confirmed use cases

**Simplifying Design**
- Reduce number of abstractions and interfaces
- Combine components with single, clear responsibility
- Choose explicit over configurable approaches

**Separating Concerns**
- Extract mixed responsibilities into focused components
- Isolate external tool integration from business logic
- Separate configuration from implementation logic

**Focusing on Immediate Needs**
- Implement minimum viable solution for current requirements
- Defer optimization until performance issues demonstrated
- Remove speculative features and flexibility

## Programming Architecture Guidelines

### Core Programming Principles

Following KISS and YAGNI principles, PanGenomePlus implementation must prioritize simplicity and immediate needs over theoretical flexibility.

#### 1. Fewer Classes, More Functions
- **Approach**: Function-first architecture for bioinformatics processing
- **Rationale**: Data transformation pipelines benefit from functional composition
- **Implementation**: Use module-level functions for extraction, clustering, output generation
- **Example**: `extract_proteins(genome_file)` instead of `ProteinExtractor().extract()`

#### 2. Flatter Structure Initially
- **Approach**: Avoid deep module hierarchies and inheritance chains
- **Rationale**: Complexity should emerge from actual requirements, not anticipation
- **Implementation**: Group related functions in modules, avoid sub-packages until necessary
- **Growth strategy**: Refactor into deeper structure only when modules become unwieldy

#### 3. Composition Over Inheritance
- **Approach**: Combine simple components rather than building class hierarchies
- **Rationale**: Genomic features share data but have different behaviors
- **Implementation**: Use simple dataclasses + specialized functions
- **Example**: `Feature` dataclass + `validate_protein()`, `validate_trna()` functions

#### 4. Simple Data Structures
- **Approach**: Dataclasses and dictionaries over abstract base classes
- **Rationale**: Genomic data is naturally structured and well-understood
- **Implementation**: `@dataclass` for entities, `Dict` for mappings, `List` for collections
- **Avoid**: Abstract base classes unless you have 3+ implementations with different behaviors

#### 5. Add Abstraction Only When You Have Concrete Need
- **Principle**: No "future-proofing" or "what if" abstractions
- **Trigger**: Create abstractions only when existing code breaks or becomes unmaintainable
- **Test**: Can you name 2+ concrete implementations that would use this abstraction today?
- **Examples of premature abstraction**: Strategy patterns, plugin systems, configuration frameworks

### Recommended Module Organization

#### Simple, Function-Based Structure
```
pangenomeplus/
├── extraction.py          # Feature extraction functions
├── clustering.py          # MMseqs2 wrapper functions
├── compact_ids.py         # ID management functions
├── outputs.py             # Output format functions
├── pipeline.py            # Main orchestration
├── utils.py               # File I/O, validation
└── cli.py                 # Command-line interface
```

#### Data Modeling Strategy
```python
# Simple, clear data structures
@dataclass
class Feature:
    compact_id: str
    genome_id: str
    contig: str
    start: int
    end: int
    strand: str
    sequence: str
    feature_type: str  # P, T, R, I, C

# Function-based processing
def extract_proteins(genome_file: str) -> List[Feature]: ...
def cluster_sequences(features: List[Feature]) -> ClusterResult: ...
```

### Implementation Anti-Patterns to Avoid

#### ❌ Over-Abstraction
```python
# WRONG: Complex hierarchy for simple data
class GenomeFeature(ABC):
    @abstractmethod
    def validate(self) -> bool: ...

class ProteinFeature(GenomeFeature): ...
class tRNAFeature(GenomeFeature): ...
```

#### ❌ Strategy Patterns for Single Use Cases
```python
# WRONG: Strategy for known, fixed requirement
class IDGenerator(ABC): ...
class Base36IDGenerator(IDGenerator): ...
```

#### ❌ Plugin Systems for Known Tools
```python
# WRONG: Plugin system for 5 known external tools
class ToolPlugin(ABC): ...
class ProdigalPlugin(ToolPlugin): ...
```

#### ✅ KISS-Compliant Alternatives
```python
# RIGHT: Simple dataclass + functions
@dataclass
class Feature: ...

def validate_protein(feature: Feature) -> bool: ...
def generate_compact_id(feature_type: str, counter: int) -> str: ...
def run_prodigal(genome_file: str) -> List[Feature]: ...
```

### Growth Strategy

#### Phase 1: Start Simple
- Single `Feature` dataclass for all genomic features
- Module-level functions for each external tool integration
- Dictionary-based compact ID management
- Direct MMseqs2 function calls

#### Phase 2: Refactor When Necessary
- **Trigger**: When current approach becomes difficult to maintain
- **Evidence**: Code duplication, unclear responsibilities, performance issues
- **Approach**: Extract classes only for components with complex state management

#### Phase 3: Add Complexity Only for Concrete Benefits
- **Inheritance**: Only when you have genuine behavioral polymorphism
- **Abstraction**: Only when you have multiple implementations with shared interface
- **Configuration**: Only when users need different behaviors

### External Tool Integration Philosophy

#### KISS Approach: Direct Function Calls
```python
def run_prodigal(genome_file: str, mode: str = "normal") -> List[Feature]:
    """Direct tool execution with minimal wrapping"""
    cmd = f"prodigal -i {genome_file} -f gff -o output.gff"
    subprocess.run(cmd, shell=True, check=True)
    return parse_prodigal_gff("output.gff")
```

#### Avoid Over-Engineering
- No tool abstraction layers unless you support multiple tools for same task
- No configuration frameworks unless users need different tool behaviors
- No plugin systems unless tools are truly interchangeable

### Code Organization Principles

#### Module Responsibility
- **extraction.py**: All external tool integrations
- **clustering.py**: All MMseqs2 operations
- **outputs.py**: All output format generation
- **pipeline.py**: Orchestration and workflow control

#### Function Design
- Pure functions where possible (sequence operations, validation)
- Clear input/output types with type hints
- Single responsibility per function
- Minimal side effects (file I/O isolated where possible)

## System Requirements

### Operating Environment
- **Operating System**: Linux (any modern distribution)
- **Memory**: 8GB minimum, scales linearly with dataset size
- **Storage**: Local filesystem with sufficient space for intermediate files
- **Python**: Version 3.8 or higher

### External Tool Dependencies
- **Prodigal**: Gene prediction (any recent version)
- **MMseqs2**: Sequence clustering (version 13+)
- **tRNAscan-SE**: tRNA detection (version 2.0+)
- **Barrnap**: rRNA detection (version 0.9+)
- **MINCED**: CRISPR detection (version 0.4+)
- **Tools must be accessible via PATH or specified via environment variables**

## Input Data Requirements

### Genome File Specifications
- **Format**: FASTA files (.fasta, .fna, .fa extensions accepted)
- **Quality expectation**: Complete, high-quality genomes (user responsibility to verify)
- **Content**: Assembled genomic sequences (NO raw reads accepted)
- **File structure**: Single directory containing all genome FASTA files

### Optional Pre-existing Annotations
- **GFF3 files**: Can bypass Prodigal if complete annotations already available
- **Naming convention**: GFF files must match genome FASTA filenames
- **Format**: Standard GFF3 with gene, CDS, and other feature annotations

### Data Quality Assumptions
- Genomes are complete or near-complete assemblies
- Minimal contamination and assembly errors
- Suitable for production pangenome analysis
- Quality control performed upstream of PanGenomePlus

## Configuration and Parameter Management

### Complete Command-Line Interface Specification

**Total Parameters**: 24 (18 core analysis + 6 feature selection)

#### Core Analysis Parameters (18)

**Prodigal Gene Prediction (2 parameters)**
- `--prodigal-mode {normal,meta}` - Organism type (default: normal)
- `--translation-table {auto,11,4,1-25}` - Genetic code (default: auto)

**MMseqs2 Clustering (6 parameters)**
- `--clustering-identity FLOAT` - Sequence identity threshold (default: 0.8)
- `--clustering-coverage FLOAT` - Coverage threshold (default: 0.8)
- `--clustering-sensitivity FLOAT` - Sensitivity level (default: 4.0)
- `--coverage-mode {query,target,bidirectional}` - Coverage calculation method (default: bidirectional)
- `--cluster-mode {set-cover,connected-component,greedy}` - Clustering algorithm (default: set-cover)
- `--max-seqs INT` - Maximum sequences per query for memory management (default: 10000)

**tRNA Detection (3 parameters)**
- `--trna-model {bacteria,archaea,eukaryota,mitochondrial}` - Organism-specific model (default: bacteria)
- `--trna-score-cutoff FLOAT` - Score threshold for tRNA detection (default: 20.0)
- `--trna-search-mode {relaxed,normal,strict}` - Search stringency (default: normal)

**rRNA Detection (2 parameters)**
- `--rrna-kingdom {bac,arc,euk,mito}` - Kingdom-specific model (default: bac)
- `--rrna-evalue FLOAT` - E-value threshold (default: 1e-6)

**CRISPR Detection (3 parameters)**
- `--crispr-min-repeats INT` - Minimum repeat count for CRISPR array (default: 3)
- `--crispr-repeat-length-range MIN,MAX` - Repeat length bounds (default: 23,47)
- `--crispr-spacer-length-range MIN,MAX` - Spacer length bounds (default: 26,50)

**Performance Control (2 parameters)**
- `--threads INT` - Number of parallel threads (default: auto-detect)
- `--memory-limit SIZE` - Memory limit for MMseqs2 operations (default: auto)

#### Feature Selection Parameters (6)

**Skip Flags**
- `--skip-trna` - Don't analyze tRNAs (skips tRNAscan-SE)
- `--skip-rrna` - Don't analyze rRNAs (skips Barrnap)
- `--skip-crispr` - Don't analyze CRISPR elements (skips MINCED)
- `--skip-intergenic` - Don't extract intergenic regions

**Analysis Subset Flags**
- `--protein-only` - Analyze only protein-coding genes (skips all non-coding)
- `--non-coding-only` - Analyze only non-coding features (skips proteins)

#### Organism Preset Configurations

**`--preset bacteria`** (default for prokaryotes)
```bash
--prodigal-mode normal --translation-table 11 --trna-model bacteria --rrna-kingdom bac
```

**`--preset archaea`** (for archaeal genomes)
```bash
--prodigal-mode normal --translation-table 11 --trna-model archaea --rrna-kingdom arc
```

**`--preset eukaryota`** (for eukaryotic genomes)
```bash
--prodigal-mode normal --translation-table 1 --trna-model eukaryota --rrna-kingdom euk
```

**`--preset metagenomic`** (for metagenomic assemblies)
```bash
--prodigal-mode meta --translation-table auto --trna-search-mode relaxed
```

### Pipeline Behavior with Feature Selection

#### Feature Selection Impact on Pipeline Stages

**Annotation Stage Modifications**
- Skip flags prevent execution of corresponding external tools
- Tool output directories created only for requested feature types
- Conditional dependency checking based on selected features

**Extraction Stage Adaptations**
- Feature type filtering applied during GFF3 parsing
- Compact ID assignment maintains type prefixes for selected features only
- Coordinate ordering preserved regardless of feature selection

**Clustering Stage Variations**
- FASTA files generated only for requested feature types
- MMseqs2 execution limited to relevant sequence types
- Bidirectional clustering applied only to selected nucleotide features

**Output Format Consistency**
- Transformer output includes only analyzed feature types in coordinate order
- Presence/absence matrices reflect selected features only
- Summary statistics calculated on analyzed subset
- Clear labeling indicates which features were processed

#### Feature Selection Examples

**Protein-only analysis**: `--protein-only`
- Pipeline: Annotation → Protein extraction → Protein clustering → Family assignment → Output
- Tools used: Prodigal only
- Output: Core/accessory/cloud protein families

**Non-coding comprehensive**: `--non-coding-only`
- Pipeline: Full annotation → Non-coding extraction → Multi-type clustering → Output
- Tools used: tRNAscan-SE, Barrnap, MINCED (Prodigal still needed for intergenic boundaries)
- Output: tRNA, rRNA, CRISPR, and intergenic families

**Selective analysis**: `--skip-crispr --skip-rrna`
- Pipeline: Partial annotation → Protein + tRNA + intergenic → Clustering → Output
- Tools used: Prodigal, tRNAscan-SE
- Output: Protein, tRNA, and intergenic families

### Parameter Interaction and Validation Rules

#### Mutually Exclusive Parameters
- `--protein-only` and `--non-coding-only` cannot be used together
- `--skip-*` flags ignored when using subset flags (`--protein-only`, `--non-coding-only`)

#### Dependency Requirements
- `--skip-intergenic` requires at least one other feature type
- `--non-coding-only` requires Prodigal for intergenic region boundary detection
- Organism presets override individual tool-specific parameters

#### Parameter Range Validation
- Identity and coverage thresholds: 0.0-1.0
- Sensitivity values: 1.0-8.0 (MMseqs2 range)
- Thread count: 1 to system maximum
- Repeat/spacer length ranges: positive integers, min < max

### External Tool Configuration
- **Auto-detection**: Attempt to find tools via PATH
- **Manual specification**: Override via environment variables (PRODIGAL_PATH, MMSEQS_PATH, etc.)
- **Version checking**: Warn if tool versions are outdated or incompatible
- **Conditional execution**: Only required tools executed based on feature selection

## Checkpoint and Recovery System

### Resume Functionality
- **Checkpoint stages**: After annotation, extraction, clustering, families
- **Implementation**: Detect existing intermediate files and resume from last completed stage
- **State preservation**: Save progress indicators and intermediate results
- **Manual restart**: Force restart from specific stage via command-line flag

### Error Handling Strategy
- **Failed genomes**: Log errors and continue processing remaining genomes
- **Tool failures**: Report specific external tool errors with diagnostic information
- **Graceful degradation**: Continue analysis with subset of successfully processed genomes
- **No timeouts**: Allow external tools to complete without artificial time limits

### Recovery Protocols
- **Partial completion**: Preserve successfully processed results
- **Incremental processing**: Add new genomes to existing analysis
- **Failed stage recovery**: Restart from specific pipeline stage
- **Resource exhaustion**: Clear guidance on memory and disk space requirements

## Output Specifications

### File Naming Conventions
- **Transformer format**: `pangenome_transformer.txt`
- **Presence/absence matrix**: `presence_absence_matrix.csv`
- **Family summary**: `family_summary.tsv`
- **Mapping tables**: `compact_id_mappings.json`
- **Log files**: `pangenomeplus_YYYYMMDD_HHMMSS.log`

### Directory Structure
```
output_directory/
├── transformer/
│   └── pangenome_transformer.txt
├── matrices/
│   ├── presence_absence_matrix.csv
│   └── family_summary.tsv
├── mappings/
│   └── compact_id_mappings.json
└── logs/
    └── pangenomeplus_YYYYMMDD_HHMMSS.log
```

### Output Format Validation
- **Transformer format**: Verify coordinate ordering and family completeness
- **Matrix format**: Validate genome and family name consistency
- **Mapping integrity**: Ensure bidirectional compact ID mapping accuracy

## Installation and Deployment

### Installation Methods
- **Python package**: `pip install pangenomeplus`
- **Development install**: `pip install -e .` from source directory
- **Conda environment**: Optional conda package with dependencies
- **Container deployment**: Docker container for isolated execution

### External Tool Setup
- **User responsibility**: Install and configure external tools independently
- **PATH requirements**: Ensure all tools accessible via system PATH
- **Version compatibility**: Use recent versions of all external tools
- **Installation verification**: Test tool execution before running pipeline

### Environment Setup
- **Python environment**: Virtual environment recommended
- **System dependencies**: Standard build tools for potential compilation
- **Permission requirements**: Read/write access to input and output directories

## Logging and Monitoring

### Logging Levels
- **ERROR**: Pipeline failures and critical issues requiring intervention
- **WARNING**: Quality concerns, missing features, non-fatal issues
- **INFO**: Progress updates, stage completion, major milestones
- **DEBUG**: Detailed processing information for troubleshooting

### Progress Reporting
- **Stage completion**: Clear indication of pipeline stage progress
- **Genome processing**: Count of completed vs. total genomes
- **Feature statistics**: Counts of extracted and clustered features
- **Time estimates**: Simple completion time projections for long stages

### Log Format
- **Standard Python logging**: Consistent timestamp and level formatting
- **Structured information**: Key metrics in parseable format
- **Error context**: Sufficient detail for troubleshooting failures
- **Performance metrics**: Basic timing and resource usage information

## Implementation Priorities

### Phase 1: Core Functionality
1. Fix tRNA and intergenic region clustering
2. Implement proper sequence extraction for all feature types
3. Validate clustering parameters per feature type

### Phase 2: Output Formats
1. Generate all required output formats
2. Add comprehensive summary statistics
3. Create visualization components

### Phase 3: Optimization
1. Performance profiling and optimization
2. Memory usage improvements
3. Parallel processing where beneficial

### Phase 4: Validation
1. Literature comparison studies
2. Benchmark against existing tools
3. User acceptance testing
