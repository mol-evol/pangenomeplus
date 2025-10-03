# PanGenomePlus Technical Debt and Code Quality Issues

## Priority Levels
- 🔴 **CRITICAL**: Blocking issues or major bugs
- 🟠 **HIGH**: Significant technical debt affecting maintainability
- 🟡 **MEDIUM**: Code quality issues that should be addressed
- 🟢 **LOW**: Minor improvements and optimizations

---

## 🔴 CRITICAL Issues

### 1. Massive Functions Exceeding Maintainability Limits
These functions are too large and complex, making them difficult to test, debug, and maintain.

#### `pangenome_analysis.py:618-1035` - generate_comprehensive_markdown_report()
- **Size**: 417 lines (4x recommended maximum)
- **Complexity**: Multiple nested sections, hardcoded HTML/CSS
- **Fix**: Extract each report section into separate functions:
  ```python
  # Suggested refactor:
  def generate_summary_section(stats: Dict) -> str
  def generate_diversity_section(stats: Dict) -> str
  def generate_visualization_section(figures: List) -> str
  def generate_comprehensive_markdown_report(...):
      sections = [
          generate_summary_section(stats),
          generate_diversity_section(stats),
          generate_visualization_section(figures)
      ]
      return "\n\n".join(sections)
  ```

#### `pipeline.py:833-1209` - process_genomes()
- **Size**: 376 lines
- **Complexity**: Handles all genome processing in single function
- **Fix**: Break into logical stages:
  ```python
  def validate_genome_files(genome_dir: Path) -> List[Path]
  def run_annotation_tools(genome_file: Path, config: Config) -> AnnotationResult
  def extract_features(annotations: AnnotationResult) -> List[Feature]
  def process_genomes(genome_dir: Path, config: Config):
      genomes = validate_genome_files(genome_dir)
      for genome in genomes:
          annotations = run_annotation_tools(genome, config)
          features = extract_features(annotations)
  ```

#### `cli.py:22-289` - create_parser()
- **Size**: 267 lines
- **Complexity**: All arguments defined in single function
- **Fix**: Group related arguments:
  ```python
  def add_clustering_args(parser: ArgumentParser)
  def add_tool_args(parser: ArgumentParser)
  def add_output_args(parser: ArgumentParser)
  ```

---

## 🟠 HIGH Priority Issues

### 2. Generic Exception Handling Masking Errors
Found 14 instances of `except Exception as e:` that hide specific error types.

#### Locations:
- `extraction.py:178, 240, 302, 364, 426` - Tool execution catches
- `pipeline.py:512, 623, 734, 845, 956` - Processing catches
- `pangenome_analysis.py:234, 345, 456, 567` - Analysis catches

#### Example Problem:
```python
# Current (BAD):
try:
    result = run_prodigal(genome)
except Exception as e:
    logger.error(f"Error: {e}")
    return None

# Should be (GOOD):
try:
    result = run_prodigal(genome)
except FileNotFoundError as e:
    logger.error(f"Genome file not found: {e}")
    raise
except subprocess.CalledProcessError as e:
    logger.error(f"Prodigal failed: {e}")
    return None
```

### 3. Code Duplication - Genome Loading Pattern
Same genome loading logic repeated 5 times across modules.

#### Locations:
- `extraction.py:45-67` - load_genome()
- `pipeline.py:123-145` - _load_genome_sequence()
- `pangenome_analysis.py:234-256` - read_genome()
- `clustering.py:78-100` - get_genome_sequence()
- `outputs.py:345-367` - load_genome_file()

#### Fix: Create single utility function:
```python
# utils/genome_io.py
def load_genome(filepath: Path) -> SeqRecord:
    """Single source of truth for genome loading"""
    if not filepath.exists():
        raise FileNotFoundError(f"Genome not found: {filepath}")
    return SeqIO.read(filepath, "fasta")
```

### 4. Magic Numbers Throughout Codebase
Hard-coded values that should be constants.

#### Translation Tables (`extraction.py`):
```python
# Current:
translation_table = 11  # Line 234
if organism == "archaea":
    translation_table = 11  # Line 245
elif organism == "mitochondria":
    translation_table = 2  # Line 247

# Should be:
from constants import TRANSLATION_TABLES
translation_table = TRANSLATION_TABLES[organism]
```

#### Clustering Thresholds (`clustering.py`):
```python
# Current:
identity = 0.8  # Line 123
coverage = 0.8  # Line 124
sensitivity = 4.0  # Line 125

# Should be:
from constants import DEFAULT_CLUSTERING_PARAMS
identity = config.get('identity', DEFAULT_CLUSTERING_PARAMS['identity'])
```

#### Visualization Parameters (`pangenome_analysis.py`):
```python
# Current:
figsize=(12, 6)  # Line 456
top_n = 50  # Line 567
dpi = 300  # Line 678

# Should be:
from constants import VIZ_PARAMS
figsize = VIZ_PARAMS['default_figsize']
```

---

## 🟡 MEDIUM Priority Issues

### 5. Resource Cleanup Issues
Files and processes not properly cleaned up.

#### `extraction.py:145` - Temporary files not deleted
```python
# Missing cleanup:
temp_file = Path(f"/tmp/prodigal_{uuid4()}.gff")
# ... processing ...
# Never deleted!

# Fix with context manager:
with tempfile.NamedTemporaryFile(suffix='.gff') as temp_file:
    # ... processing ...
    # Automatically cleaned up
```

#### `clustering.py:234` - MMseqs2 databases not removed
```python
# Current: Creates DB but no cleanup
mmseqs_db = output_dir / "mmseqs_db"

# Should use try/finally or context manager
```

### 6. Inconsistent Path Handling
Mix of string and Path operations.

#### Locations:
- `pipeline.py`: Uses Path objects
- `extraction.py`: Uses strings
- `cli.py`: Mixed usage

#### Fix: Standardize on pathlib.Path:
```python
# Bad:
genome_file = args.genome_dir + "/" + filename

# Good:
genome_file = Path(args.genome_dir) / filename
```

### 7. Missing Input Validation
No validation before processing.

#### `pipeline.py:456` - No genome format check
```python
# Current: Assumes valid FASTA
genome = load_genome(genome_file)

# Should validate:
def validate_fasta(filepath: Path) -> bool:
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
        return len(records) > 0
    except Exception:
        return False
```

### 8. Hardcoded Output Paths
Output filenames scattered through code.

#### Examples:
- `outputs.py:123`: `"presence_absence_matrix.csv"`
- `outputs.py:234`: `"family_summary.tsv"`
- `outputs.py:345`: `"pangenome_transformer.txt"`

#### Fix: Centralize in constants:
```python
# constants.py
OUTPUT_FILES = {
    'presence_absence': 'presence_absence_matrix.csv',
    'summary': 'family_summary.tsv',
    'transformer': 'pangenome_transformer.txt'
}
```

---

## 🟢 LOW Priority Issues

### 9. Missing Type Hints
Functions lack proper type annotations.

#### Priority Functions for Type Hints:
- `extraction.py`: All run_* functions
- `clustering.py`: cluster_sequences()
- `pipeline.py`: process_genomes()

### 10. Logging Inconsistencies
Mixed logging approaches.

#### Issues:
- Some modules use logger, others use print()
- No consistent log format
- Missing debug logging in key areas

### 11. Missing Docstrings
Many functions lack documentation.

#### Priority for Documentation:
- All public API functions
- Complex internal functions
- Functions with non-obvious parameters

### 12. No Progress Bars for Long Operations
User has no feedback during long runs.

#### Add Progress Bars To:
- Genome processing loop
- MMseqs2 clustering
- Feature extraction

### 13. Configuration Not Centralized
Settings scattered across modules.

#### Create Config System:
```python
# config.py
@dataclass
class PanGenomeConfig:
    clustering_identity: float = 0.8
    clustering_coverage: float = 0.8
    min_intergenic_length: int = 50
    core_threshold: float = 0.95
    accessory_lower: float = 0.15
```

---

## Implementation Priority Order

### Phase 1: Critical Fixes (Week 1)
1. Refactor massive functions
2. Fix generic exception handling
3. Eliminate code duplication

### Phase 2: High Priority (Week 2)
4. Replace magic numbers with constants
5. Fix resource cleanup
6. Standardize path handling

### Phase 3: Medium Priority (Week 3)
7. Add input validation
8. Centralize output paths
9. Add missing type hints

### Phase 4: Polish (Week 4)
10. Improve logging
11. Add docstrings
12. Add progress bars
13. Create configuration system

---

## Quick Wins (Can be done immediately)

1. **Add constants.py module** - Central location for all magic numbers
2. **Fix tRNAscan-SE -Q flag** - Already done, verify in all code paths
3. **Remove print() statements** - Replace with proper logging
4. **Delete commented-out code** - Clean up codebase
5. **Fix obvious typos** - "succesfully" → "successfully", etc.

---

## Testing Requirements

After addressing these issues, add tests for:
- Each refactored function
- Error handling paths
- Edge cases (empty genomes, missing files)
- Configuration validation
- Output format verification

---

## Notes

- Many issues are interconnected (e.g., refactoring large functions will naturally improve exception handling)
- Start with high-impact, low-effort fixes (constants, cleanup)
- Consider using automated tools (black, isort, mypy) for consistent formatting
- Document changes in CHANGELOG as you go
