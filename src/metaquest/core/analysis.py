"""
MetaQuest Core Analysis Module — FASTQ-only v6.0.0
===================================================
Delegates to pipeline.runner for actual orchestration.
This module provides backward-compatible entry points.
"""

import hashlib
import json
from pathlib import Path
from dataclasses import asdict, replace

from ..exceptions import ConfigError
from ..settings import get_config, load_config
from ..pipeline.runner import build_default_pipeline
from ..pipeline.context import PipelineContext


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def _validate_resume(output_dir, input_files, config, workflow):
    metadata_path = output_dir / "analysis_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    recorded_inputs = metadata.get("input_files", [])
    current_inputs = [
        {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in input_files
    ]
    if recorded_inputs != current_inputs:
        raise ConfigError("--resume input files do not match the recorded run")
    if metadata.get("effective_config") != _jsonable(asdict(config)):
        raise ConfigError("--resume configuration does not match the recorded run")
    recorded_workflow = metadata.get("workflow", {})
    for key, value in workflow.items():
        if recorded_workflow.get(key) != value:
            raise ConfigError(f"--resume workflow option does not match: {key}")


def run_analysis(input_file, output_dir, cli_args=None):
    """
    Main analysis entry point — backward compatible.

    Args:
        input_file: FASTQ file path(s).
        output_dir: Output directory path.
        cli_args: Command-line arguments namespace.
    """
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Load config (from cli_args.config if available, else defaults)
    config_path = getattr(cli_args, "config", None)
    db_dir = getattr(cli_args, "db_dir", None)
    if config_path or db_dir:
        config = load_config(
            Path(config_path) if config_path else None,
            db_dir=Path(db_dir) if db_dir else None,
        )
    else:
        config = get_config()

    # Preserve the existing CLI while ensuring command-line values actually
    # reach the immutable configuration used by pipeline stages.
    annotation_threads = getattr(cli_args, "annotation_threads", None)
    if annotation_threads is not None:
        config = replace(
            config,
            annotation=replace(config.annotation, threads=annotation_threads),
        )
    plot_formats = getattr(cli_args, "plot_formats", None)
    if plot_formats:
        config = replace(
            config,
            reporting=replace(config.reporting, plot_formats=tuple(plot_formats)),
        )

    # Determine read mode
    read_mode = "paired"
    if cli_args:
        if getattr(cli_args, "single", False):
            read_mode = "single"
        elif getattr(cli_args, "interleaved", False):
            read_mode = "interleaved"

    # Normalize input files
    if isinstance(input_file, str):
        input_files = [Path(input_file)]
    elif isinstance(input_file, (list, tuple)):
        input_files = [Path(f) for f in input_file]
    else:
        input_files = [Path(input_file)]

    skip_annotation = bool(
        getattr(cli_args, "taxonomy_only", False)
        or getattr(cli_args, "skip_annotation", False)
    )
    skip_functional = getattr(cli_args, "skip_functional", False)
    skip_abundance = bool(getattr(cli_args, "skip_abundance", False))
    resume = getattr(cli_args, "resume", False)

    if resume:
        _validate_resume(
            output_dir_path,
            input_files,
            config,
            {
                "read_mode": read_mode,
                "taxonomy_only": skip_annotation,
                "skip_functional": skip_functional,
                "skip_abundance": skip_abundance,
                "low_memory": bool(getattr(cli_args, "low_memory", False)),
            },
        )

    # Build and run pipeline
    ctx = PipelineContext(
        config=config,
        input_files=input_files,
        output_dir=output_dir_path,
        read_mode=read_mode,
        skip_annotation=skip_annotation,
        skip_functional=skip_functional,
        skip_abundance=skip_abundance,
        low_memory=getattr(cli_args, "low_memory", False),
        resume=resume,
    )
    ctx.metadata["validation"] = {
        "status": getattr(cli_args, "validation_status", "not_run"),
        "strict": bool(getattr(cli_args, "strict_validation", False)),
        "bypassed": bool(getattr(cli_args, "skip_validation", False)),
        "reports": _jsonable(getattr(cli_args, "validation_results", [])),
    }

    pipeline = build_default_pipeline(
        config,
        skip_annotation=skip_annotation,
        skip_functional=skip_functional,
        skip_abundance=skip_abundance,
    )

    ctx = pipeline.run(ctx)
