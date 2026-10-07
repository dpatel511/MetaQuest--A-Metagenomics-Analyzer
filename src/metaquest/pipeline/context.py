"""
Pipeline Context — carries state between stages.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metaquest.settings import MetaQuestConfig


@dataclass
class AssemblyResult:
    contigs_fasta: Path
    total_contigs: int = 0
    total_bases: int = 0
    n50: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    l50: int = 0
    n90: int = 0
    l90: int = 0
    gc_percent: float = 0.0
    contigs_ge_1000: int = 0
    contigs_ge_5000: int = 0
    contigs_ge_10000: int = 0
    lengths: list[int] = field(default_factory=list)


@dataclass
class ClassificationResult:
    kraken_report: Path
    bracken_file: Path
    species_count: int = 0
    species_assigned_reads: int = 0
    total_reads: int = 0
    kraken_classified_reads: int = 0
    unclassified_reads: int = 0
    observed_read_length: int = 0
    bracken_read_length: int = 0


@dataclass
class AnnotationResult:
    gene_prediction_dir: Path
    protein_file: Path | None = None
    functional_annotations: Path | None = None
    functional_category_summary: Path | None = None
    gene_count: int = 0
    annotated_count: int = 0
    functional_reused: bool = False


@dataclass
class AbundanceResult:
    gene_abundance: Path
    functional_abundance: Path | None = None
    summary: dict[str, Any] = field(default_factory=dict)
    reused: bool = False


@dataclass
class PipelineContext:
    """Accumulates results as stages execute sequentially."""

    config: MetaQuestConfig
    input_files: list[Path]
    output_dir: Path
    read_mode: str = "paired"  # paired | single | interleaved
    skip_annotation: bool = False
    skip_functional: bool = False
    skip_abundance: bool = False
    low_memory: bool = False
    resume: bool = False
    analysis_input_files: list[Path] = field(default_factory=list)
    preprocessing: dict[str, Any] | None = None

    # Filled by stages
    classification: ClassificationResult | None = None
    assembly: AssemblyResult | None = None
    annotation: AnnotationResult | None = None
    abundance: AbundanceResult | None = None

    # Metadata for reproducibility
    metadata: dict[str, Any] = field(default_factory=dict)
    completed_stages: list[str] = field(default_factory=list)
