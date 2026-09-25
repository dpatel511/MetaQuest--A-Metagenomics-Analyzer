"""Gene abundance stage — read mapping to contigs and TPM normalization."""

from __future__ import annotations

import logging

from metaquest.exceptions import AnnotationError
from metaquest.pipeline.context import AbundanceResult, PipelineContext

logger = logging.getLogger(__name__)


def run_abundance_stage(ctx: PipelineContext) -> PipelineContext:
    """Quantify predicted genes and aggregate abundance by functional term."""
    from metaquest.core.abundance import run_gene_abundance
    from metaquest.io.output_formatter import get_formatter

    if ctx.annotation is None:
        raise AnnotationError("Gene prediction must complete before gene abundance")

    reads = ctx.analysis_input_files or ctx.input_files
    prediction_dir = ctx.annotation.gene_prediction_dir
    config = ctx.config.abundance
    logger.info("Mapping reads to contigs with BBMap for gene abundance")
    try:
        gene_table, functional_table, summary, reused = run_gene_abundance(
            reads,
            prediction_dir / "contigs.stable.fasta",
            prediction_dir / "genes.gff3",
            ctx.output_dir,
            functional_table=ctx.annotation.functional_annotations,
            read_mode=ctx.read_mode,
            threads=config.threads,
            min_identity=config.min_identity,
            ambiguous=config.ambiguous,
            term_attribution=config.term_attribution,
            java_memory=config.java_memory,
            reuse=ctx.resume,
        )
    except AnnotationError:
        raise
    except Exception as exc:
        raise AnnotationError(f"Gene abundance failed: {exc}", cause=exc) from exc

    ctx.abundance = AbundanceResult(
        gene_abundance=gene_table,
        functional_abundance=functional_table,
        summary=summary,
        reused=reused,
    )
    formatter = get_formatter()
    if reused:
        formatter.info("Reusing completed gene abundance mapping")
    if summary.get("mapping_rate", 0) < 0.5:
        formatter.warning(
            f"Only {summary.get('mapping_rate', 0):.1%} of {summary['count_unit']} mapped "
            "to the assembly; gene abundance describes the assembled fraction only"
        )
    logger.info(
        "Gene abundance complete: %.2f%% %s mapped, %.2f%% assigned to genes",
        summary.get("mapping_rate", 0) * 100,
        summary.get("count_unit", "reads"),
        summary.get("gene_assignment_rate", 0) * 100,
    )
    return ctx
