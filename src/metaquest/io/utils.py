"""
MetaQuest Utilities Module - CLEANED v2.0
==========================================

Core utility functions for assembly, validation, and system checks.
Progress parsing removed - using spinners instead.
"""

import subprocess
from pathlib import Path
import importlib

def _check_command(name, version_cmd="--version"):
    """Helper to check for a single command-line tool."""
    try:
        if name == 'seqtk':
            subprocess.run([name], capture_output=True, text=True, check=False)
            return True
        
        cmd = name.split() + [version_cmd] if version_cmd else name.split()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_system_check(
    formatter,
    *,
    config=None,
    taxonomy_only: bool = False,
    require_interleaved: bool = False,
    require_functional: bool = False,
    require_abundance: bool = False,
):
    """
    Runs a comprehensive check of all dependencies using the formatter for output.
    """
    errors = []
    warnings = []
    
    formatter.info("Performing system-wide checks...")

    # Command-Line Tools
    formatter.substep("Checking command-line tools...")
    tools = {'fastp': '--version', 'kraken2': '--version', 'bracken': '--version'}
    if not taxonomy_only:
        tools['megahit'] = '--version'
    if require_functional:
        tools.update({'diamond': 'version', 'emapper.py': '--version'})
    if require_abundance and not taxonomy_only:
        tools['bbmap.sh'] = '--version'
    if require_interleaved:
        tools['reformat.sh'] = None
    for tool, cmd in tools.items():
        if not _check_command(tool, cmd):
            errors.append(f"Tool not found: '{tool}'. Please install it, e.g., via 'conda install -c bioconda {tool}'.")

    # Python Packages
    formatter.substep("Checking Python packages...")
    packages = {
        'Bio': 'biopython',
        'pandas': 'pandas',
        'numpy': 'numpy',
        'matplotlib': 'matplotlib',
        'seaborn': 'seaborn',
    }
    if not taxonomy_only:
        packages['pyrodigal'] = 'pyrodigal'
    for pkg, install_name in packages.items():
        try:
            importlib.import_module(pkg)
        except ImportError:
            errors.append(f"Python package not found: '{install_name}'. Please install it, e.g., 'pip install {install_name}'.")

    # Databases
    formatter.substep("Checking database files...")
    if config is None:
        from metaquest.settings import get_config
        config = get_config()
    required_dbs = {
        "taxonomy": (config.databases.kraken_db / "hash.k2d",),
    }
    if require_functional:
        required_dbs["functional"] = (
            config.databases.functional_dir / "eggnog.db",
            config.databases.functional_dir / "eggnog.taxa.db",
            config.databases.functional_dir / "eggnog_proteins.dmnd",
        )
    for database, paths in required_dbs.items():
        missing = [path for path in paths if not path.exists()]
        if missing:
            errors.append(
                f"Required {database} database files are missing: "
                f"{', '.join(str(path) for path in missing)}. Install them with "
                f"'metaquest databases --db-dir {config.databases.base_dir} "
                f"--database {database}'."
            )

    # Final Report
    if warnings:
        for warning in warnings:
            formatter.warning(warning)
    
    if errors:
        formatter.error(
            "System check failed due to missing critical dependencies.",
            solutions=errors + ["Please resolve the issues above before running an analysis."]
        )
        raise SystemExit(1)

def split_interleaved(interleaved_fastq: str, output_dir: Path, formatter) -> list:
    """
    Splits an interleaved FASTQ into two files using reformat.sh.
    """
    r1 = output_dir / "split_R1.fastq.gz"
    r2 = output_dir / "split_R2.fastq.gz"

    cmd = [
        "reformat.sh", f"in={interleaved_fastq}",
        f"out1={r1}", f"out2={r2}", "qin=33", "qout=33", "overwrite=true",
    ]
    
    returncode, _, stderr = formatter.run_subprocess(
        cmd, "Splitting interleaved file", show_command=True
    )

    if returncode != 0:
        formatter.error(
            "reformat.sh failed to split the file.",
            solutions=["Ensure BBTools is installed and the input file is a valid FASTQ."]
        )
        raise RuntimeError(f"reformat.sh failed. Stderr: {stderr}")
    
    return [str(r1), str(r2)]
