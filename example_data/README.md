# Example Dataset

This directory contains 3 *Escherichia coli* genomes for testing PanGenomePlus.

## Files

- `ecoli_example_001.fasta` - *E. coli* genome 1 (~4.6 MB)
- `ecoli_example_002.fasta` - *E. coli* genome 2 (~4.5 MB)
- `ecoli_example_003.fasta` - *E. coli* genome 3 (~4.6 MB)

## Quick Test

Run PanGenomePlus on these genomes:

```bash
pangenomeplus --genome-dir example_data/ --output-dir example_output/
```

## Expected Results

**Runtime**: 2-5 minutes (depending on system)

**Expected output**:
- Core protein families: ~3,500-4,000
- Accessory protein families: ~500-1,000
- tRNA families: ~40-50
- rRNA families: 3-6
- CRISPR families: Variable (may be 0)

## Purpose

Use this small dataset to:
- Verify PanGenomePlus installation
- Test external tool dependencies
- Familiarize yourself with output formats
- Benchmark performance on your system
