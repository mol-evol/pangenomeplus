"""Simplified CLI interface for PanGenomePlus."""

import argparse
import sys
from pathlib import Path
from typing import Optional, Literal

from pangenomeplus.utils.config import (
    load_configuration, create_default_configuration, save_configuration
)
from pangenomeplus.core.exceptions import ConfigurationError, PipelineError
def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for PanGenomePlus CLI."""
    parser = argparse.ArgumentParser(
        prog='pangenomeplus',
        description='PanGenomePlus: Adaptive scale pangenome analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python -m pangenomeplus genomes/ output/

  # With custom parameters
  python -m pangenomeplus genomes/ output/ --sensitivity 6.0 --use-gpu --threads 8

  # Generate config template
  python -m pangenomeplus --init-config config.yaml

  # Analyze dataset only
  python -m pangenomeplus genomes/ --analyze-only

  # Run downstream analysis
  python -m pangenomeplus --downstream-analysis results/ --output analysis/
        """.strip()
    )

    # Positional arguments (optional)
    parser.add_argument('genomes_dir', nargs='?', type=Path,
                       help='Directory containing genome FASTA files')
    parser.add_argument('output_dir', nargs='?', type=Path,
                       help='Output directory for pipeline results')

    # Special modes
    special_group = parser.add_argument_group('Special modes')
    special_group.add_argument('--init-config', type=Path, metavar='FILE',
                              help='Create default configuration file and exit')
    special_group.add_argument('--analyze-only', action='store_true',
                              help='Only analyze dataset characteristics and exit')
    special_group.add_argument('--downstream-analysis', type=Path, metavar='RESULTS_DIR',
                              help='Run downstream analysis on existing pipeline results')

    # Core options
    core_group = parser.add_argument_group('Core options')
    core_group.add_argument('--config', '-c', type=Path,
                           help='Configuration file (YAML or JSON)')
    core_group.add_argument('--output', '-o', type=Path,
                           help='Output directory for downstream analysis')
    core_group.add_argument('--no-validate', action='store_true',
                           help='Skip input validation for faster startup')
    core_group.add_argument('--verbose', '-v', action='store_true',
                           help='Enable verbose logging')
    core_group.add_argument('--resume', action='store_true', default=True,
                           help='Resume from existing checkpoints (default: True)')
    core_group.add_argument('--no-resume', dest='resume', action='store_false',
                           help='Do not resume from existing checkpoints')
    core_group.add_argument(
        '--stop-after', choices=['annotation', 'clustering', 'families', 'analysis'],
                           help='Stop pipeline after specified stage')
    core_group.add_argument(
        '--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO',
                           help='Logging verbosity level (default: INFO)')

    # MMseqs2 clustering parameters
    clustering_group = parser.add_argument_group('Clustering parameters')
    clustering_group.add_argument('--protein-identity', type=float, metavar='FLOAT',
                                 help='Protein clustering identity threshold (0.0-1.0, default: 0.8)'
                                )
    clustering_group.add_argument('--protein-coverage', type=float, metavar='FLOAT',
                                 help='Protein clustering coverage threshold (0.0-1.0, default: 0.8)')
    clustering_group.add_argument('--sensitivity', type=float, metavar='FLOAT',
                                 help='MMseqs2 search sensitivity (1.0-8.0, default: 7.5)')
    clustering_group.add_argument('--use-gpu', action='store_true',
                                 help='Enable GPU acceleration for clustering')
    clustering_group.add_argument('--gpu-memory', type=str, metavar='SIZE',
                                 help='GPU memory limit (e.g., 8G, 16G)')

    # Annotation parameters
    annotation_group = parser.add_argument_group('Annotation parameters')
    annotation_group.add_argument('--genetic-code', type=int, metavar='INT',
                                 help='Genetic code table for Prodigal (1-31, default: 11)')
    annotation_group.add_argument('--annotation-tools', type=str, metavar='TOOLS',
                                 help='Comma-separated annotation tools (prodigal,trnascan,barrnap,minced)')

    # Classification parameters
    classification_group = parser.add_argument_group('Classification parameters')
    classification_group.add_argument('--core-threshold', type=float, metavar='FLOAT',
                                     help='Core gene threshold (0.0-1.0, default: 0.95)')
    classification_group.add_argument('--cloud-threshold', type=float, metavar='FLOAT',
                                     help='Cloud gene threshold (0.0-1.0, default: 0.15)')

    # Resource parameters
    resource_group = parser.add_argument_group('Resource parameters')
    resource_group.add_argument('--threads', '-t', type=int, metavar='INT',
                               help='Number of threads to use (default: 4)')
    resource_group.add_argument('--memory-limit', type=str, metavar='SIZE',
                               help='Memory limit (e.g., 32G, 64G)')

    # Output parameters
    output_group = parser.add_argument_group('Output parameters')
    output_group.add_argument('--output-formats', type=str, metavar='FORMATS',
                             help='Comma-separated output formats (transformer,presence_absence,roary,fasta)')
    output_group.add_argument('--compress-output', action='store_true',
                             help='Compress large output files')

    # Downstream analysis options
    downstream_group = parser.add_argument_group('Downstream analysis options')
    downstream_group.add_argument('--genome-pattern', default=r'(E_coli_\d+)',
                                 help='Regex pattern to extract genome IDs from gene IDs (default: %(default)s)')
    downstream_group.add_argument('--sample-size', type=int, default=10,
                                 help='Number of genomes to sample for analysis (default: %(default)s)')

    # Version
    parser.add_argument('--version', action='version', version='PanGenomePlus 1.0.0')

    return parser


def cli() -> None:
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()

    try:
        # Handle special modes first
        if args.init_config:
            _handle_init_config(args.init_config)
            return

        if args.analyze_only:
            if not args.genomes_dir:
                print("Error: GENOMES_DIR required for --analyze-only", file=sys.stderr)
                sys.exit(1)
            _handle_analyze_only(args.genomes_dir, args.sample_size)
            return

        if args.downstream_analysis:
            output_path = args.output or Path('feature_analysis_results')
            _handle_downstream_analysis(args.downstream_analysis, output_path, args.genome_pattern)
            return

        # Main pipeline mode - require both arguments
        if not args.genomes_dir or not args.output_dir:
            print("Error: Both GENOMES_DIR and OUTPUT_DIR are required for pipeline execution", file=sys.stderr)
            print("Use --help to see all available options", file=sys.stderr)
            sys.exit(1)

        # Validate genomes_dir exists
        if not args.genomes_dir.exists():
            print(f"Error: Genomes directory does not exist: {args.genomes_dir}", file=sys.stderr)
            sys.exit(1)
        if not args.genomes_dir.is_dir():
            print(f"Error: Genomes path is not a directory: {args.genomes_dir}", file=sys.stderr)
            sys.exit(1)

        # Load or create configuration
        if args.config:
            if not args.config.exists():
                print(f"Error: Configuration file does not exist: {args.config}", file=sys.stderr)
                sys.exit(1)
            pipeline_config = load_configuration(args.config)
            print(f"Loaded configuration from {args.config}")
        else:
            pipeline_config = create_default_configuration()
            if args.verbose:
                print("Using default configuration")

        # Override configuration with CLI parameters
        _apply_cli_overrides(pipeline_config, {
            'protein_identity': args.protein_identity,
            'protein_coverage': args.protein_coverage,
            'sensitivity': args.sensitivity,
            'use_gpu': args.use_gpu,
            'gpu_memory': args.gpu_memory,
            'genetic_code': args.genetic_code,
            'annotation_tools': args.annotation_tools,
            'core_threshold': args.core_threshold,
            'cloud_threshold': args.cloud_threshold,
            'threads': args.threads,
            'memory_limit': args.memory_limit,
            'output_formats': args.output_formats,
            'compress_output': args.compress_output
        })

        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)

        # Import and run pipeline
        from pangenomeplus.pipeline import run_pangenome_pipeline

        print("Starting pangenome analysis...")
        print(f"  Input genomes: {args.genomes_dir}")
        print(f"  Output directory: {args.output_dir}")
        if args.verbose:
            print(f"  Log level: {args.log_level}")
            print(f"  Validation: {'disabled' if args.no_validate else 'enabled'}")

        results = run_pangenome_pipeline(
            genomes_dir=args.genomes_dir,
            output_dir=args.output_dir,
            config=pipeline_config,
            resume=args.resume,
            stop_after=args.stop_after,
            validate_inputs=not args.no_validate,
            log_level=args.log_level
        )

        print("Pipeline completed successfully!")
        print(f"Processed {results.get('n_genomes', 0)} genomes")
        print(f"Generated {results.get('n_families', 0)} gene families")

    except (ConfigurationError, PipelineError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def _handle_init_config(output_path: Path) -> None:
    """Handle --init-config mode."""
    try:
        config = create_default_configuration()
        save_configuration(config, output_path)
        print(f"Created default configuration: {output_path}")
    except Exception as e:
        print(f"Error creating configuration: {e}", file=sys.stderr)
        sys.exit(1)


def _handle_analyze_only(genomes_dir: Path, sample_size: int) -> None:
    """Handle --analyze-only mode."""
    try:
        from pangenomeplus.pipeline import analyze_dataset_and_choose_strategy

        print(f"Analyzing dataset in {genomes_dir}...")

        dataset_info = analyze_dataset_and_choose_strategy(
            genomes_dir=genomes_dir,
            sample_size=sample_size
        )

        print("\nDataset Analysis Results:")
        print(f"  Number of genomes: {dataset_info['n_genomes']}")
        print(f"  Available memory: {dataset_info['memory_available_gb']} GB")
        print(f"  Description: {dataset_info['description']}")

    except Exception as e:
        print(f"Error analyzing dataset: {e}", file=sys.stderr)
        sys.exit(1)


def _handle_downstream_analysis(pipeline_results_dir: Path, output_path: Path, genome_pattern: str) -> None:
    """Handle --downstream-analysis mode."""
    try:
        from pangenomeplus.modules.downstream_analysis import run_feature_type_analysis

        print(f"Running feature-type analysis on results in {pipeline_results_dir}...")

        # Construct paths to required files
        family_summary_file = pipeline_results_dir / "families" / "family_summary.tsv"
        gene_to_family_file = pipeline_results_dir / "families" / "gene_to_family.tsv"

        # Validate required files exist
        if not family_summary_file.exists():
            print(f"Error: Required file not found: {family_summary_file}", file=sys.stderr)
            sys.exit(1)

        if not gene_to_family_file.exists():
            print(f"Error: Required file not found: {gene_to_family_file}", file=sys.stderr)
            sys.exit(1)

        # Run analysis
        stats = run_feature_type_analysis(
            family_summary_file=family_summary_file,
            gene_to_family_file=gene_to_family_file,
            output_dir=output_path,
            genome_id_pattern=genome_pattern
        )

        print("\nFeature-type analysis completed!")
        print(f"Results saved to: {output_path}")

        # Display summary statistics
        print("\nSummary:")
        for feature_type, type_stats in stats.items():
            if feature_type in ['protein', 'intergenic']:
                name = 'Protein-coding' if feature_type == 'protein' else 'Intergenic'
                print(f"  {name}: {type_stats['total_families']} families")
                print(f"    Core: {type_stats['core_families']} ({type_stats['core_percentage']:.1f}%)")
                print(f"    Shell: {type_stats['shell_families']} ({type_stats['shell_percentage']:.1f}%)")
                print(f"    Cloud: {type_stats['cloud_families']} ({type_stats['cloud_percentage']:.1f}%)")

    except Exception as e:
        print(f"Error running feature analysis: {e}", file=sys.stderr)
        sys.exit(1)


def _apply_cli_overrides(config: dict, cli_params: dict) -> None:
    """Apply CLI parameter overrides to configuration."""

    # MMseqs2 clustering parameters
    if cli_params['protein_identity'] is not None:
        config.setdefault('clustering', {}).setdefault('protein', {})['identity'] = cli_params['protein_identity']
    if cli_params['protein_coverage'] is not None:
        config.setdefault('clustering', {}).setdefault('protein', {})['coverage'] = cli_params['protein_coverage']
    if cli_params['sensitivity'] is not None:
        for feature_type in ['protein', 'trna', 'rrna', 'intergenic']:
            config.setdefault('clustering', {}).setdefault(feature_type, {})['sensitivity'] = cli_params['sensitivity']

    # GPU parameters
    if cli_params['use_gpu']:
        config.setdefault('clustering', {}).setdefault('gpu', {})['use_gpu'] = True
    if cli_params['gpu_memory'] is not None:
        config.setdefault('clustering', {}).setdefault('gpu', {})['memory_limit'] = cli_params['gpu_memory']

    # Annotation parameters
    if cli_params['genetic_code'] is not None:
        config.setdefault('annotation', {}).setdefault('prodigal', {})['genetic_code'] = cli_params['genetic_code']
    if cli_params['annotation_tools'] is not None:
        tools_list = [tool.strip() for tool in cli_params['annotation_tools'].split(',')]
        config.setdefault('annotation', {})['tools'] = tools_list

    # Classification thresholds
    if cli_params['core_threshold'] is not None:
        config.setdefault('classification', {})['soft_core_threshold'] = cli_params['core_threshold']
    if cli_params['cloud_threshold'] is not None:
        config.setdefault('classification', {})['cloud_threshold'] = cli_params['cloud_threshold']

    # Resource parameters
    if cli_params['threads'] is not None:
        config.setdefault('resources', {})['threads'] = cli_params['threads']
    if cli_params['memory_limit'] is not None:
        config.setdefault('resources', {})['memory_limit'] = cli_params['memory_limit']

    # Output parameters
    if cli_params['output_formats'] is not None:
        formats_list = [fmt.strip() for fmt in cli_params['output_formats'].split(',')]
        config.setdefault('output', {})['formats'] = formats_list
    if cli_params['compress_output']:
        config.setdefault('output', {}).setdefault('transformer', {})['compress'] = True


if __name__ == '__main__':
    cli()