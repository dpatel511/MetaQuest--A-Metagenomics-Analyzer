"""
MetaQuest CLI - Command Line Interface (FASTQ-only version)
============================================================

A professional command-line interface for the MetaQuest metagenomics pipeline.
Now exclusively supports FASTQ input for streamlined workflow.

Author: MetaQuest Development Team
"""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path
from textwrap import dedent
from .exceptions import MetaQuestError
from .io.output_formatter import Colors, OutputFormatter, set_formatter
from .config import __version__

# ============================================================================
# APP METADATA
# ============================================================================

__app_name__ = "MetaQuest"
__tagline__ = "An Integrated Metagenomics Analysis Pipeline"
__author__ = "MetaQuest Development Team"


def _normalize_global_options(argv):
    """Allow display/config options on either side of the subcommand."""
    switches = {'--debug', '--quiet', '--no-color', '--verbose', '--low-memory'}
    global_args = []
    command_args = []
    index = 0
    while index < len(argv):
        value = argv[index]
        if value in switches:
            global_args.append(value)
        elif value == '--config' and index + 1 < len(argv):
            global_args.extend(argv[index:index + 2])
            index += 1
        elif value.startswith('--config='):
            global_args.append(value)
        else:
            command_args.append(value)
        index += 1
    return global_args + command_args


def _prepare_output_directory(args, parser):
    """Enforce explicit fresh, forced, or resumed output behavior."""
    if args.command != "run":
        return
    output = Path(args.output)
    if output.exists() and not output.is_dir():
        parser.error(f"output path is not a directory: {output}")
    exists_with_content = output.exists() and any(output.iterdir())
    if args.resume:
        if not exists_with_content or not (output / "analysis_metadata.json").is_file():
            parser.error("--resume requires an existing MetaQuest output directory")
        return
    if exists_with_content and not args.force:
        parser.error(
            f"output directory is not empty: {output}. "
            "Choose a new directory, use --resume, or use --force."
        )
    if exists_with_content:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
        backup = output.with_name(f"{output.name}.metaquest-backup-{stamp}")
        output.rename(backup)
    output.mkdir(parents=True, exist_ok=True)


# ============================================================================
# ARGUMENT GROUP SETUP FUNCTIONS
# ============================================================================

def setup_fastq_validation_args(parser):
    """
    Add FASTQ-specific validation arguments to a parser.
    
    Configures quality control thresholds for FASTQ validation including:
    - Minimum quality scores
    - Sequence count requirements
    - Overrepresentation detection
    """
    val_group = parser.add_argument_group(
        'FASTQ Validation Options',
        'Configure quality control parameters for FASTQ input files'
    )
    
    val_group.add_argument(
        '-q', '--min-quality',
        type=int,
        default=20,
        metavar='SCORE',
        help="Minimum mean Phred quality score threshold (default: 20)"
    )
    
    val_group.add_argument(
        '-n', '--min-sequences',
        type=int,
        default=100,
        metavar='NUM',
        help="Minimum number of sequences required for analysis (default: 100)"
    )
    
    val_group.add_argument(
        '--overrep-threshold',
        type=float,
        default=0.1,
        metavar='PCT',
        help="Fraction threshold for flagging overrepresented sequences (default: 0.1)"
    )
    val_group.add_argument(
        '--strict-validation',
        action='store_true',
        help="Treat quality warnings as fatal validation errors",
    )

def setup_annotation_args(parser):
    """
    Add functional annotation control arguments to a parser.
    
    Configures the thread allocation for provisional functional annotation.
    
    Pyrodigal runs in metagenomic mode on contigs of at least 200 bp.
    """
    annot_group = parser.add_argument_group(
        'Annotation Control Options',
        'Control Pyrodigal gene prediction and functional annotation'
    )
    
    annot_group.add_argument(
        '--annotation-threads',
        type=int,
        default=8,
        metavar='NUM',
        help="Number of CPU threads for functional annotation (default: 8)"
    )

# ============================================================================
# VALIDATION HANDLER
# ============================================================================

def handle_validation(file_paths, args, formatter):
    """
    Run comprehensive FASTQ file validation and return validation status.
    
    Args:
        file_paths: List of FASTQ file paths to validate
        args: Command line arguments containing validation parameters
        formatter: OutputFormatter instance for consistent logging
        
    Returns:
        bool: True if all files are valid, False if any file fails validation
    """
    formatter.section_header("FILE VALIDATION")
    
    from .io.file_validator import FileValidator

    validator = FileValidator()
    
    # Configure FASTQ validation parameters
    validator.quality_threshold = args.min_quality
    validator.min_sequences = args.min_sequences
    validator.overrep_threshold = args.overrep_threshold
    validator.strict = args.strict_validation
    
    formatter.info("Validation Configuration:")
    formatter.result({
        'Quality threshold': f"Q{args.min_quality}",
        'Minimum sequences': f"{args.min_sequences:,}",
        'Overrepresentation cutoff': f"{args.overrep_threshold:.1%}"
    }, indent=2)

    all_valid = True
    reports = []
    total_files = len(file_paths)
    
    # Validate each file individually
    for i, file_path in enumerate(file_paths, 1):
        if total_files > 1:
            formatter.info(f"Validating file {i}/{total_files}: {Path(file_path).name}")
        else:
            formatter.info(f"Validating file: {Path(file_path).name}")
        
        is_valid, report = validator.validate_and_analyze(file_path, 'fastq')
        if report:
            reports.append(report)
        
        if not is_valid:
            all_valid = False
            formatter.info(f"Validation failed: {file_path}")
        else:
            formatter.success(f"Validation passed: {file_path}")
    
    if len(file_paths) == 2:
        paired_valid, paired_error = validator.validate_pair(*file_paths)
        if not paired_valid:
            all_valid = False
            formatter.error(paired_error)
        else:
            formatter.success("Paired FASTQ identifiers and counts are synchronized")

    args.validation_results = reports
    args.validation_status = "passed" if all_valid else "failed"
    print()  # Add spacing
    
    if total_files > 1:
        formatter.info(f"Validation Summary: {total_files} file(s) processed")
    
    if all_valid:
        formatter.success("All files validated successfully")
        formatter.info("Files meet quality standards and are ready for analysis")
    else:
        formatter.info("Validation failed for one or more files")
        formatter.info("Please address quality issues before proceeding")
    
    return all_valid


# ============================================================================
# CUSTOM HELP FORMATTER
# ============================================================================

class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Enhanced help text formatter with improved readability.
    
    Features:
    - Extended help position for better alignment
    - Consistent column widths
    - Cleaner option string formatting
    """
    
    def __init__(self, prog):
        super().__init__(prog, max_help_position=35, width=100)
    
    def _format_action_invocation(self, action):
        """Format command-line option strings with consistent styling."""
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(f'{option_string}')
                parts[-1] += f' {args_string}'
            return ', '.join(parts)


# ============================================================================
# PARSER CREATION
# ============================================================================

def create_parser():
    """
    Create and configure the main argument parser with all subcommands.
    
    Returns:
        argparse.ArgumentParser: Fully configured parser with all commands and options
    """
    
    description = dedent(f"""
        MetaQuest {__version__} · research-use short-read metagenomics

        Validate reads, profile taxonomy, assemble contigs, predict genes,
        and produce stable descriptive reports.

        Quick start:
          metaquest check
          metaquest run --single reads.fastq.gz --output results/
    """)
    
    parser = argparse.ArgumentParser(
        description=description,
        prog="metaquest",
        formatter_class=CustomHelpFormatter,
        epilog="Run 'metaquest COMMAND --help' for command-specific options."
    )
    
    # Global options
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'MetaQuest version {__version__}'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug mode with full diagnostic output (includes all commands and stderr)'
    )
    output_mode = parser.add_mutually_exclusive_group()
    output_mode.add_argument(
        '--quiet',
        action='store_true',
        help='Show errors only'
    )
    parser.add_argument(
        '--no-color',
        action='store_true',
        help='Disable ANSI colors and terminal animation'
    )
    output_mode.add_argument(
        '--verbose',
        action='store_true',
        help='Show additional progress and diagnostic context'
    )

    parser.add_argument(
        '--config',
        metavar='FILE',
        help='Path to YAML configuration file (overrides defaults). Generate one with: metaquest init-config'
    )
    parser.add_argument(
        '--low-memory',
        action='store_true',
        help='Use Kraken2 memory mapping to reduce RAM use'
    )
    
    # Create subcommand parsers
    subparsers = parser.add_subparsers(
        dest='command',
        required=True,
        title='Commands',
        metavar='COMMAND'
    )

    # ========================================================================
    # CHECK COMMAND
    # ========================================================================
    parser_check = subparsers.add_parser(
        'check',
        help='Verify system dependencies and database status',
        description=dedent("""
            System Dependency Verification
            ==============================
            
            Performs comprehensive validation of:
                • Stable command-line tools (Kraken2, Bracken, MEGAHIT,
                  Pyrodigal, and DIAMOND)
                • Stable Python package dependencies
                • Reference database files and indices
            
            Run this before your first analysis to ensure proper configuration.
            Identifies missing dependencies with installation suggestions.
        """),
        formatter_class=CustomHelpFormatter
    )
    parser_check.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root override'
    )

    # ========================================================================
    # VALIDATE COMMAND
    # ========================================================================
    parser_validate = subparsers.add_parser(
        'validate',
        help='Validate FASTQ input files without running analysis',
        description=dedent("""
            FASTQ File Validation
            =====================
            
            Quality control checks for FASTQ files:
                • Format validation (proper headers, sequences)
                • Quality score analysis
                • Sequence count verification
                • Overrepresentation detection
                • File integrity checks
            
            Recommended to run before analysis to catch issues early.
        """),
        formatter_class=CustomHelpFormatter
    )
    
    validate_mode = parser_validate.add_mutually_exclusive_group(required=True)
    validate_mode.add_argument(
        '--single',
        metavar='READS.fastq',
        help="Single-end FASTQ file"
    )
    validate_mode.add_argument(
        '--paired',
        nargs=2,
        metavar=('R1.fastq', 'R2.fastq'),
        help="Paired-end FASTQ files"
    )
    validate_mode.add_argument(
        '--interleaved',
        metavar='INTERLEAVED.fastq',
        help="Interleaved paired-end FASTQ file"
    )
    setup_fastq_validation_args(parser_validate)

    # ========================================================================
    # ANALYZE COMMAND
    # ========================================================================
    parser_analyze = subparsers.add_parser(
        'run',
        aliases=['analyze'],
        help='Run the short-read analysis pipeline',
        description=dedent("""
            Complete Metagenomics Analysis Pipeline
            ========================================
            
            Executes the full MetaQuest pipeline on FASTQ reads:
            
            Pipeline Steps:
                1. Taxonomic Classification (Kraken2 + Bracken)
                2. Metagenomic Assembly (MEGAHIT)
                3. Gene Prediction (Pyrodigal metagenomic mode)
                4. Functional Annotation (eggNOG-mapper)
                5. Gene Abundance (BBMap read mapping, TPM)
                6. Descriptive Report Generation
            
            Use --taxonomy-only for rapid taxonomic classification, or
            --skip-functional to stop after gene prediction, or
            --skip-abundance to skip read mapping and gene abundance.
        """),
        formatter_class=CustomHelpFormatter
    )
    
    parser_analyze.set_defaults(command='run')

    # Add all standard arguments
    parser_analyze.add_argument(
        '-o', '--output',
        default='results',
        metavar='DIR',
        help="Output directory for results (default: results)"
    )
    parser_analyze.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root override'
    )
    parser_analyze.add_argument(
        '--skip-validation',
        action='store_true',
        help="Skip input file validation (not recommended for production use)"
    )
    output_mode = parser_analyze.add_mutually_exclusive_group()
    output_mode.add_argument(
        '--force',
        action='store_true',
        help="Move an existing output directory to a timestamped backup and start fresh",
    )
    output_mode.add_argument(
        '--resume',
        action='store_true',
        help="Resume from an existing MetaQuest output directory",
    )
    parser_analyze.add_argument(
        '--taxonomy-only',
        action='store_true',
        help="Run only Kraken2/Bracken taxonomy and descriptive reporting"
    )
    parser_analyze.add_argument(
        '--skip-annotation',
        action='store_true',
        help="Deprecated alias for --taxonomy-only"
    )
    parser_analyze.add_argument(
        '--skip-functional',
        action='store_true',
        help="Run assembly and Pyrodigal but skip eggNOG functional annotation"
    )
    parser_analyze.add_argument(
        '--skip-abundance',
        action='store_true',
        help="Skip read mapping and gene/functional abundance estimation"
    )
    
    # Add argument groups
    setup_annotation_args(parser_analyze)
    setup_fastq_validation_args(parser_analyze)
    parser_analyze.add_argument(
        '--plot-formats', nargs='+', choices=('svg', 'png', 'pdf'),
        help='Figure formats (default from configuration: svg png)',
    )
    
    fastq_mode = parser_analyze.add_mutually_exclusive_group(required=True)
    fastq_mode.add_argument(
        '--single',
        metavar='READS.fastq',
        help="Single-end FASTQ file"
    )
    fastq_mode.add_argument(
        '--paired',
        nargs=2,
        metavar=('R1.fastq', 'R2.fastq'),
        help="Paired-end FASTQ files (forward and reverse reads)"
    )
    fastq_mode.add_argument(
        '--interleaved',
        metavar='INTERLEAVED.fastq',
        help="Interleaved paired-end FASTQ file"
    )

    # ========================================================================
    # COMPARE COMMAND
    # ========================================================================
    subparsers.add_parser(
        'compare',
        help='Reserved for a future validated comparative workflow',
        description=(
            'The experimental comparison implementation was removed. '
            'This command name is retained for compatibility.'
        ),
        formatter_class=CustomHelpFormatter,
    )
    # ========================================================================
    # SETUP-DB COMMAND
    # ========================================================================
    parser_setup_db = subparsers.add_parser(
        'databases',
        aliases=['setup-db'],
        help='Inspect or install reference databases',
        description=dedent("""
            Database Setup Utility
            ======================
            
            Installs versioned reference data outside the source repository.
            Use --list to inspect availability and size before downloading.

            Examples:
                metaquest databases --list
                metaquest databases --db-dir /data/metaquest --database taxonomy
                metaquest databases --db-dir /data/metaquest --database functional
                metaquest check --db-dir /data/metaquest

            Use these commands instead of running Kraken database builders or
            eggNOG-mapper database-download scripts directly.
        """),
        formatter_class=CustomHelpFormatter
    )

    parser_setup_db.set_defaults(command='databases')

    # ========================================================================
    # INIT-CONFIG COMMAND
    # ========================================================================
    parser_init_config = subparsers.add_parser(
        'init-config',
        help='Generate a default configuration YAML file',
        description=dedent("""
            Generate Default Configuration
            ===============================

            Creates a metaquest.yaml configuration file with all default values.
            Edit this file to customize pipeline behavior without CLI flags.

            Usage:
                metaquest init-config
                metaquest init-config -o my_config.yaml
        """),
        formatter_class=CustomHelpFormatter
    )
    parser_init_config.add_argument(
        '-o', '--output',
        default='metaquest.yaml',
        metavar='FILE',
        help='Output path for the config file (default: metaquest.yaml)'
    )

    db_options = parser_setup_db.add_argument_group('Database Selection')
    db_options.add_argument(
        '--db-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Database root (default: METAQUEST_DB_DIR, config, or ./databases)'
    )
    db_options.add_argument(
        '--database',
        choices=('taxonomy', 'functional', 'amr', 'virulence'),
        metavar='NAME',
        help="Install one database using MetaQuest's managed downloader"
    )
    db_options.add_argument(
        '--list',
        action='store_true',
        help='List database releases, sizes, and installation status'
    )
    db_options.add_argument(
        '--force',
        action='store_true',
        help='Replace an existing invalid or older installation'
    )

    return parser

# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """
    Main entry point for the MetaQuest CLI.
    
    Orchestrates:
        • Argument parsing and validation
        • Output formatter initialization
        • Command execution routing
        • Error handling and user feedback
        • Runtime tracking and reporting
    """
    
    parser = create_parser()
    args = parser.parse_args(_normalize_global_options(sys.argv[1:]))
    _prepare_output_directory(args, parser)
    
    # Determine verbosity level from command-line flags
    verbosity = (
        'debug' if args.debug else
        'verbose' if args.verbose else
        'silent' if args.quiet else
        'standard'
    )
    
    # Setup log file path if output directory specified (not for init-config)
    log_file = None
    if hasattr(args, 'output') and args.output and args.command not in ('init-config', 'check'):
        log_file = Path(args.output) / 'metaquest.log'
    
    # Initialize global output formatter
    formatter = OutputFormatter(
        verbosity=verbosity,
        log_file=log_file,
        append_log=bool(getattr(args, 'resume', False)),
    )
    if args.no_color:
        formatter.colors_enabled = False
        Colors.disable()
    set_formatter(formatter)
    
    # Show banner (unless version or help requested)
    if not any(arg in sys.argv for arg in ['-v', '--version', '-h', '--help']):
        formatter.banner(__app_name__, __version__, __tagline__)

    start_time = time.monotonic()

    try:
        # ====================================================================
        # COMMAND: CHECK
        # ====================================================================
        if args.command == 'check':
            from .io.utils import run_system_check

            formatter.section_header("SYSTEM DEPENDENCY CHECK")
            formatter.info("Verifying installation of required tools and databases...")
            
            from .settings import load_config
            run_system_check(
                formatter,
                config=load_config(
                    Path(args.config) if args.config else None,
                    db_dir=args.db_dir,
                ),
                require_functional=True,
                require_abundance=True,
            )
            
            formatter.success("System check completed successfully")
            formatter.info("All required dependencies are properly configured")
            sys.exit(0)

        # ====================================================================
        # FILE PATH VALIDATION (for analyze and validate commands)
        # ====================================================================
        if args.command in ['run', 'validate']:
            file_paths = []
            
            if args.single:
                file_paths = [args.single]
            elif args.paired:
                file_paths = args.paired
            elif args.interleaved:
                file_paths = [args.interleaved]

            # Check file existence before proceeding
            for f in file_paths:
                if not Path(f).exists():
                    formatter.error(
                        f"Input file not found: {f}",
                        solutions=[
                            "Check the file path for typos",
                            "Verify the file exists in the current directory",
                            "Use absolute path instead of relative path",
                            "Check file permissions"
                        ]
                    )
                    sys.exit(1)

        # ====================================================================
        # COMMAND: VALIDATE
        # ====================================================================
        if args.command == 'validate':
            is_valid = handle_validation(file_paths, args, formatter)
            sys.exit(0 if is_valid else 1)
        
        # ====================================================================
        # COMMAND: ANALYZE
        # ====================================================================
        elif args.command == 'run':
            from .core.analysis import run_analysis
            from .io.utils import run_system_check

            formatter.section_header("Preflight")
            taxonomy_only = bool(args.taxonomy_only or args.skip_annotation)
            if args.skip_annotation:
                formatter.warning("--skip-annotation is deprecated; use --taxonomy-only")
            if taxonomy_only and args.skip_functional:
                formatter.warning("--skip-functional has no effect with --taxonomy-only")

            # System check
            with formatter.spinner("Verifying required tools and databases"):
                from .settings import load_config
                check_config = load_config(
                    Path(args.config) if args.config else None,
                    db_dir=args.db_dir,
                )
                run_system_check(
                    formatter,
                    config=check_config,
                    taxonomy_only=taxonomy_only,
                    require_interleaved=bool(args.interleaved),
                    require_functional=not taxonomy_only and not args.skip_functional,
                    require_abundance=not taxonomy_only and not args.skip_abundance,
                )
            formatter.success("Tools and databases are ready")

            # File validation (unless skipped)
            if not args.skip_validation:
                is_valid = handle_validation(file_paths, args, formatter)
                if not is_valid:
                    formatter.error(
                        "Analysis aborted due to validation failure",
                        solutions=[
                            "Fix validation errors reported above",
                            "Improve input file quality if needed",
                            "Use --skip-validation flag to bypass (not recommended)"
                        ]
                    )
                    sys.exit(1)
                formatter.success("Input validation completed - files meet quality standards")
            else:
                args.validation_status = "skipped"
                args.validation_results = []
                formatter.warning("Skipping file validation (--skip-validation flag detected)")
                formatter.info("Proceeding with unvalidated data")
            
            formatter.section_header("Run configuration")
            formatter.result({
                'Input': ', '.join(file_paths),
                'Read type': 'Single-end' if args.single else 'Paired-end' if args.paired else 'Interleaved',
                'Database': str(check_config.databases.base_dir),
                'Workflow': (
                    'Taxonomy only'
                    if taxonomy_only
                    else 'Taxonomy + assembly + gene prediction'
                    if args.skip_functional
                    else 'Taxonomy + assembly + gene prediction + eggNOG'
                ),
                'Gene abundance': (
                    'Skipped' if taxonomy_only or args.skip_abundance else 'BBMap + TPM'
                ),
                'Memory mode': 'Low memory' if args.low_memory else 'Default',
                'Threads': str(args.annotation_threads),
                'Output': args.output,
            })

            run_analysis(file_paths, args.output, args)

            elapsed = time.monotonic() - start_time
            summary_path = Path(args.output) / 'analysis_summary.json'
            summary = json.loads(summary_path.read_text(encoding='utf-8'))
            metrics = {}
            taxonomy = summary.get('taxonomy', {})
            assembly = summary.get('assembly', {})
            annotation = summary.get('annotation', {})
            if taxonomy:
                metrics['Reported taxa'] = taxonomy.get('reported_taxa', 0)
                metrics['Kraken classification rate'] = (
                    f"{taxonomy.get('classification_rate', 0):.2%}"
                )
                unit = 'reads' if args.single else 'fragments'
                metrics[f'Species-assigned {unit}'] = (
                    f"{taxonomy.get('species_assigned_reads', 0):,}"
                )
            preprocessing = summary.get('preprocessing', {})
            if preprocessing.get('after'):
                metrics['Reads retained after fastp'] = (
                    f"{preprocessing.get('retained_fraction', 0):.2%}"
                )
            if assembly:
                metrics['Assembly contigs'] = f"{assembly.get('total_contigs', 0):,}"
                metrics['Assembly N50'] = f"{assembly.get('n50', 0):,} bp"
            if annotation:
                metrics['Predicted genes'] = f"{annotation.get('predicted_genes', 0):,}"
                if not args.skip_functional:
                    metrics['eggNOG annotated genes'] = (
                        f"{annotation.get('annotated_genes', 0):,}"
                    )
            abundance = summary.get('abundance', {})
            if abundance:
                metrics[f"Mapped {abundance.get('count_unit', 'reads')}"] = (
                    f"{abundance.get('mapping_rate', 0):.2%}"
                )
                metrics['Assigned to genes'] = (
                    f"{abundance.get('gene_assignment_rate', 0):.2%}"
                )
            formatter.section_header("Results")
            formatter.result(metrics)
            formatter.success(f"Completed in {formatter._format_time(elapsed)}")
            formatter.info(f"Results → {summary_path.resolve()}")
            formatter.info(f"HTML report → {(Path(args.output) / 'report.html').resolve()}")
            formatter.info(f"Figures → {(Path(args.output) / 'plots').resolve()}")

        # ====================================================================
        # COMMAND: COMPARE
        # ====================================================================
        elif args.command == 'compare':
            formatter.section_header("COMPARATIVE ANALYSIS")
            formatter.error(
                "Comparative analysis is not available in the stable runtime",
                solutions=[
                    "Run taxonomy and functional analysis for each sample first",
                    "Comparative analysis will return after its methods are validated",
                ],
            )
            sys.exit(1)

        # ====================================================================
        # COMMAND: SETUP-DB
        # ====================================================================
        elif args.command == 'databases':
            from .database_manager import DATABASES, database_rows, install_database
            from .settings import load_config

            formatter.section_header("DATABASE SETUP")

            configured = load_config(Path(args.config)) if args.config else load_config()
            db_root = args.db_dir or configured.databases.base_dir
            formatter.result({
                'Database root': str(Path(db_root).expanduser().resolve()),
            })

            if args.list:
                for row in database_rows(db_root):
                    formatter.info(
                        f"{row['Database']}: {row['Status']} | "
                        f"{row['Release']} | {row['Size']}"
                    )
                sys.exit(0)

            if not args.database:
                formatter.error(
                    "No database selected",
                    solutions=[
                        "Run 'metaquest databases --list' to inspect available databases",
                        "Run 'metaquest databases --database taxonomy' for taxonomy data",
                        "Run 'metaquest databases --database functional' for functional data",
                    ]
                )
                sys.exit(1)

            selected = DATABASES[args.database]
            if selected.download_gb is not None:
                formatter.warning(
                    f"The {args.database} download is approximately "
                    f"{selected.download_gb:g} GB"
                )
            target = install_database(
                args.database,
                db_root,
                force=args.force,
                notify=formatter.info,
            )
            elapsed = time.time() - start_time
            formatter.success(
                f"Database setup completed in {formatter._format_time(elapsed)}"
            )
            formatter.info(f"Installed at: {target}")

        # ====================================================================
        # COMMAND: INIT-CONFIG
        # ====================================================================
        elif args.command == 'init-config':
            import shutil
            from .settings import _DEFAULT_CONFIG_PATH

            output_path = Path(args.output)
            if output_path.exists():
                formatter.warning(f"File already exists: {output_path}")
                formatter.info("Overwriting with default configuration")

            shutil.copy2(_DEFAULT_CONFIG_PATH, output_path)
            formatter.success(f"Configuration file created: {output_path}")
            formatter.info("Edit this file to customize pipeline behavior.")
            formatter.info("Usage: metaquest run --config metaquest.yaml --single reads.fq -o results/")

    except KeyboardInterrupt:
        formatter.warning("Operation interrupted by user (Ctrl+C)")
        sys.exit(130)

    except MetaQuestError as e:
        formatter.error(str(e))
        if e.suggestions:
            formatter.info("Suggestions:")
            for s in e.suggestions:
                formatter.info(f"  - {s}", indent=2)
        sys.exit(1)

    except Exception as e:
        formatter.error(
            f"Unexpected error: {e}",
            solutions=[
                "Check the log file for details",
                "Verify input file integrity",
                "Ensure sufficient disk space and memory",
                "Report issue at GitHub if problem persists"
            ],
            doc_link="https://github.com/Karudhoru/MetaQuest--A-Metagenomics-Analyzer/issues"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
