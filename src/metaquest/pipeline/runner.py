"""
Pipeline Runner — orchestrates stages in sequence.
"""

from __future__ import annotations

import gc
import hashlib
import json
import logging
import subprocess
import time
from datetime import datetime
from dataclasses import asdict
from pathlib import Path
from typing import Any

from metaquest.exceptions import MetaQuestError, PipelineStageError
from metaquest.settings import MetaQuestConfig
from metaquest.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


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


def _tool_version(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=False, timeout=15
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unavailable"
    output = f"{result.stdout}\n{result.stderr}".strip()
    return output.splitlines()[0] if output else "unknown"


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


class PipelineRunner:
    """
    Runs the MetaQuest analysis pipeline as a sequence of stages.

    Each stage receives the PipelineContext and returns it (possibly mutated).
    Stages are functions: (ctx: PipelineContext) -> PipelineContext
    """

    def __init__(self, config: MetaQuestConfig):
        self.config = config
        self._stages: list[tuple[str, Any]] = []

    def add_stage(self, name: str, stage_fn):
        """Register a stage function."""
        self._stages.append((name, stage_fn))
        return self

    def run(self, ctx: PipelineContext) -> PipelineContext:
        """Execute all registered stages in order."""
        from metaquest.io.output_formatter import get_formatter

        total = len(self._stages)
        formatter = get_formatter()
        logger.info("Pipeline starting with %d stage(s)", total)
        formatter.section_header("Pipeline")

        attempt_started = datetime.now().isoformat()
        history = []
        if ctx.resume:
            previous = _read_json(ctx.output_dir / "analysis_metadata.json")
            history = list(previous.get("attempt_history", []))
            if not history and previous:
                prior_attempt = {
                    "started_at": previous.get("started_at"),
                    "finished_at": previous.get("timestamp"),
                    "status": previous.get("status", "unknown"),
                    "stage_timings_seconds": previous.get(
                        "stage_timings_seconds", {}
                    ),
                }
                if previous.get("failed_stage"):
                    prior_attempt["failed_stage"] = previous["failed_stage"]
                    prior_attempt["failure_reason"] = previous.get("failure_reason")
                history.append(prior_attempt)
            original_started = previous.get(
                "original_started_at", previous.get("started_at", attempt_started)
            )
            original_timings = previous.get(
                "original_stage_timings_seconds",
                previous.get("stage_timings_seconds", {}),
            )
        else:
            original_started = attempt_started
            original_timings = {}

        current_attempt = {
            "started_at": attempt_started,
            "status": "running",
            "stage_timings_seconds": {},
        }
        history.append(current_attempt)
        ctx.metadata.update(
            {
                "status": "running",
                "started_at": original_started,
                "original_started_at": original_started,
                "last_attempt_started_at": attempt_started,
                "original_stage_timings_seconds": original_timings,
                "attempt_history": history,
                "stage_timings_seconds": original_timings,
            }
        )
        self._save_metadata(ctx)

        for i, (name, stage_fn) in enumerate(self._stages, 1):
            logger.info("[%d/%d] %s", i, total, name)
            formatter.info(f"[{i}/{total}] {name}")
            started = time.monotonic()
            try:
                with formatter.spinner(name):
                    ctx = stage_fn(ctx)
                ctx.completed_stages.append(name)
                elapsed = time.monotonic() - started
                elapsed = round(elapsed, 3)
                current_attempt["stage_timings_seconds"][name] = elapsed
                if not ctx.resume:
                    ctx.metadata.setdefault("stage_timings_seconds", {})[name] = elapsed
                    ctx.metadata["original_stage_timings_seconds"][name] = elapsed
                else:
                    ctx.metadata.setdefault("resume_stage_timings_seconds", {})[
                        name
                    ] = elapsed
                self._save_metadata(ctx)
                formatter.success(
                    f"{name.ljust(28)} {formatter._format_time(elapsed)}"
                )
            except PipelineStageError as exc:
                self._record_failure(ctx, name, exc)
                raise
            except MetaQuestError as exc:
                self._record_failure(ctx, name, exc)
                raise
            except Exception as e:
                self._record_failure(ctx, name, e)
                raise PipelineStageError(name, str(e), cause=e) from e
            finally:
                gc.collect()

        ctx.metadata["status"] = "completed"
        current_attempt["status"] = "completed"
        current_attempt["finished_at"] = datetime.now().isoformat()
        self._save_metadata(ctx)
        logger.info("Pipeline complete: %d stage(s) finished", len(ctx.completed_stages))
        return ctx

    def _record_failure(self, ctx: PipelineContext, stage: str, exc: Exception) -> None:
        """Persist completed stages and failure provenance before propagating."""
        ctx.metadata.update(
            {
                "status": "failed",
                "failed_stage": stage,
                "failure_reason": str(exc),
            }
        )
        history = ctx.metadata.get("attempt_history", [])
        if history:
            history[-1].update(
                {
                    "status": "failed",
                    "finished_at": datetime.now().isoformat(),
                    "failed_stage": stage,
                    "failure_reason": str(exc),
                }
            )
        self._save_metadata(ctx)

    def _save_metadata(self, ctx: PipelineContext) -> None:
        """Save pipeline metadata for reproducibility."""
        metadata = {
            "timestamp": datetime.now().isoformat(),
            "metaquest_version": _get_version(),
            "analysis_type": "fastq",
            "completed_stages": ctx.completed_stages,
            "effective_config": _jsonable(asdict(ctx.config)),
            "workflow": {
                "read_mode": ctx.read_mode,
                "taxonomy_only": ctx.skip_annotation,
                "skip_functional": ctx.skip_functional,
                "skip_abundance": ctx.skip_abundance,
                "low_memory": ctx.low_memory,
                "resume": ctx.resume,
            },
            "input_files": [
                {
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
                for path in ctx.input_files
            ],
            "tool_provenance": {
                "kraken2": _tool_version(["kraken2", "--version"]),
                "bracken": _tool_version(["bracken", "-v"]),
                "fastp": _tool_version(["fastp", "--version"]),
                "megahit": (
                    _tool_version(["megahit", "--version"]) if ctx.assembly else "not_run"
                ),
                "pyrodigal": _read_json(
                    ctx.output_dir / "gene_prediction" / "summary.json"
                ).get("tool_version", "not_run"),
                "eggnog_mapper": _read_json(
                    ctx.output_dir / "functional_annotation" / "summary.json"
                ).get("tool_version", "not_run"),
                "bbmap": _read_json(
                    ctx.output_dir / "gene_abundance" / "summary.json"
                ).get("tool_version", "not_run"),
                "kraken_database": {
                    "path": str(ctx.config.databases.kraken_db.resolve()),
                    "configured_release": ctx.config.databases.kraken_db_version,
                    "manifest": _read_json(
                        ctx.config.databases.kraken_db / "metaquest-db.json"
                    ),
                },
                "functional_database": (
                    {
                        "path": str(ctx.config.databases.functional_dir.resolve()),
                        "manifest": _read_json(
                            ctx.config.databases.functional_dir / "metaquest-db.json"
                        ),
                    }
                    if ctx.annotation and ctx.annotation.functional_annotations
                    else "not_used"
                ),
            },
        }
        if ctx.classification:
            metadata["tool_provenance"]["bracken_model_read_length"] = (
                ctx.classification.bracken_read_length
            )
        metadata.update(_jsonable(ctx.metadata))

        metadata_file = ctx.output_dir / "analysis_metadata.json"
        temporary = metadata_file.with_suffix(".json.tmp")
        with open(temporary, "w") as f:
            json.dump(metadata, f, indent=2)
            f.write("\n")
        temporary.replace(metadata_file)


def _get_version() -> str:
    try:
        from metaquest.config import __version__
        return __version__
    except ImportError:
        return "unknown"


def build_default_pipeline(
    config: MetaQuestConfig,
    skip_annotation: bool = False,
    skip_functional: bool = False,
    skip_abundance: bool = False,
) -> PipelineRunner:
    """
    Build the standard MetaQuest analysis pipeline.

    Stages:
      1. Taxonomic classification (Kraken2 + Bracken)
      2. Metagenomic assembly (MEGAHIT) [full workflow only]
      3. Gene prediction (Pyrodigal)
      4. Functional annotation (eggNOG-mapper) [unless skipped]
      5. Gene abundance (BBMap read mapping + TPM) [unless skipped]
      6. Stable reporting

    Custom pathogen detection, ML, HMM, ESM, island detection, and risk
    scoring are intentionally excluded until independently validated.
    """
    from metaquest.pipeline.stages.classification import run_classification_stage
    from metaquest.pipeline.stages.preprocessing import run_preprocessing_stage
    from metaquest.pipeline.stages.reporting import run_reporting_stage

    runner = PipelineRunner(config)
    runner.add_stage("Read Preprocessing", run_preprocessing_stage)
    runner.add_stage("Taxonomic Classification", run_classification_stage)

    if not skip_annotation:
        from metaquest.pipeline.stages.assembly import run_assembly_stage
        from metaquest.pipeline.stages.annotation import (
            run_annotation_stage,
            run_functional_annotation_stage,
        )

        runner.add_stage("Metagenomic Assembly", run_assembly_stage)
        runner.add_stage("Gene Prediction", run_annotation_stage)
        if not skip_functional:
            runner.add_stage("Functional Annotation", run_functional_annotation_stage)
        if not skip_abundance:
            from metaquest.pipeline.stages.abundance import run_abundance_stage

            runner.add_stage("Gene Abundance", run_abundance_stage)

    runner.add_stage("Reporting", run_reporting_stage)
    return runner
