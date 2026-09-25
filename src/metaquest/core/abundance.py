"""Read-mapping gene abundance and functional abundance aggregation.

Quality-controlled reads are mapped back to the stable-ID contigs written by
gene prediction. Each primary alignment that passes the identity filter is
assigned to the predicted CDS containing its alignment midpoint. For paired
data, mates are combined so each fragment contributes one count, split equally
across the distinct genes its mates hit. Gene counts are normalized to
transcripts-per-million-style values (TPM) using CDS nucleotide length, then
summed per COG, KO, EC, and GO term from the eggNOG annotation table.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from ..exceptions import AnnotationError


GENE_ABUNDANCE_COLUMNS = (
    "gene_id",
    "contig_id",
    "start",
    "end",
    "strand",
    "partial",
    "length_bp",
    "count",
    "mean_depth",
    "rpk",
    "tpm",
    "annotation_status",
)
FUNCTIONAL_NAMESPACES = ("COG", "KO", "EC", "GO")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_REFERENCE_OPS = frozenset("MDN=X")


@dataclass(frozen=True)
class Gene:
    gene_id: str
    contig_id: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str
    partial: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class GeneIndex:
    """Per-contig interval lookup for predicted CDS features."""

    def __init__(self, genes: Iterable[Gene]):
        by_contig: dict[str, list[Gene]] = defaultdict(list)
        self.genes: list[Gene] = []
        for gene in genes:
            by_contig[gene.contig_id].append(gene)
            self.genes.append(gene)
        self._contigs: dict[str, tuple[list[int], list[Gene], int]] = {}
        for contig, items in by_contig.items():
            items.sort(key=lambda gene: (gene.start, gene.end))
            self._contigs[contig] = (
                [gene.start for gene in items],
                items,
                max(gene.length for gene in items),
            )

    def overlapping(self, contig: str, start: int, end: int) -> list[Gene]:
        """Return genes overlapping the 1-based inclusive interval."""
        entry = self._contigs.get(contig)
        if entry is None:
            return []
        starts, items, max_length = entry
        index = bisect_right(starts, end) - 1
        hits = []
        while index >= 0 and starts[index] > start - max_length:
            gene = items[index]
            if gene.end >= start:
                hits.append(gene)
            index -= 1
        return hits

    def containing(self, contig: str, position: int) -> list[Gene]:
        return self.overlapping(contig, position, position)


def read_gene_gff(gff_path: Path) -> list[Gene]:
    """Parse Pyrodigal CDS features from a GFF3 file."""
    genes = []
    with Path(gff_path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attributes = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            gene_id = attributes.get("ID")
            if not gene_id:
                raise AnnotationError(f"CDS feature without an ID in {gff_path}")
            genes.append(
                Gene(
                    gene_id=gene_id,
                    contig_id=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    partial=attributes.get("partial", ""),
                )
            )
    return genes


def _alignment_metrics(cigar: str, tags: list[str]) -> tuple[int, float | None]:
    """Return (reference span, identity) for a SAM alignment."""
    operations = [(int(length), op) for length, op in _CIGAR.findall(cigar)]
    span = sum(length for length, op in operations if op in _REFERENCE_OPS)
    totals = Counter()
    for length, op in operations:
        totals[op] += length
    columns = totals["="] + totals["X"] + totals["M"] + totals["I"] + totals["D"]
    if not columns:
        return span, None
    if totals["M"] == 0:
        return span, totals["="] / columns
    for tag in tags:
        if tag.startswith("NM:i:"):
            edits = int(tag[5:])
            return span, max(0.0, 1.0 - edits / columns)
    return span, None


@dataclass
class _Mate:
    genes: frozenset[str]
    mapped: bool


class AbundanceCounter:
    """Accumulate fragment counts and aligned-base depth from SAM records."""

    def __init__(self, index: GeneIndex, *, min_identity: float, paired: bool):
        self.index = index
        self.min_identity = min_identity
        self.paired = paired
        self.counts: Counter[str] = Counter()
        self.aligned_bases: Counter[str] = Counter()
        self.stats = Counter()
        self._pending: dict[str, _Mate] = {}

    def add_sam_lines(self, lines: Iterable[str]) -> None:
        for line in lines:
            if not line or line.startswith("@"):
                continue
            self._add_record(line.rstrip("\n").split("\t"))
        self.finish()

    def _add_record(self, fields: list[str]) -> None:
        if len(fields) < 11:
            raise AnnotationError("Malformed SAM record from read mapper")
        flag = int(fields[1])
        if flag & 0x900:  # secondary or supplementary
            return
        self.stats["reads"] += 1
        mate = _Mate(frozenset(), False)
        if not flag & 0x4 and fields[2] != "*":
            span, identity = _alignment_metrics(fields[5], fields[11:])
            if identity is None or identity >= self.min_identity:
                self.stats["reads_mapped"] += 1
                start = int(fields[3])
                end = start + max(span, 1) - 1
                for gene in self.index.overlapping(fields[2], start, end):
                    overlap = min(end, gene.end) - max(start, gene.start) + 1
                    self.aligned_bases[gene.gene_id] += overlap
                midpoint = (start + end) // 2
                mate = _Mate(
                    frozenset(gene.gene_id for gene in self.index.containing(fields[2], midpoint)),
                    True,
                )
            else:
                self.stats["reads_below_identity"] += 1
        if not self.paired or not flag & 0x1:
            self._record_fragment((mate,))
            return
        name = fields[0]
        previous = self._pending.pop(name, None)
        if previous is None:
            self._pending[name] = mate
        else:
            self._record_fragment((previous, mate))

    def _record_fragment(self, mates: tuple[_Mate, ...]) -> None:
        self.stats["fragments"] += 1
        if not any(mate.mapped for mate in mates):
            return
        self.stats["fragments_mapped"] += 1
        genes = sorted(set().union(*(mate.genes for mate in mates)))
        if not genes:
            self.stats["fragments_intergenic"] += 1
            return
        self.stats["fragments_in_genes"] += 1
        share = 1.0 / len(genes)
        for gene_id in genes:
            self.counts[gene_id] += share

    def finish(self) -> None:
        """Count mates whose partner never appeared as a lone fragment."""
        pending, self._pending = self._pending, {}
        for mate in pending.values():
            self._record_fragment((mate,))


def compute_tpm(genes: list[Gene], counts: Counter[str]) -> dict[str, tuple[float, float]]:
    """Return {gene_id: (reads per kilobase, TPM)}."""
    rates = {gene.gene_id: counts.get(gene.gene_id, 0.0) / (gene.length / 1000) for gene in genes}
    total = sum(rates.values())
    return {
        gene_id: (rate, rate / total * 1_000_000 if total else 0.0)
        for gene_id, rate in rates.items()
    }


def _terms(value: str) -> list[str]:
    if not value or value == "-":
        return []
    return [term.strip() for term in value.split(",") if term.strip() and term.strip() != "-"]


def _gene_terms(row: dict[str, str]) -> dict[str, list[str]]:
    cog = sorted({letter for letter in "".join(_terms(row.get("cog_categories", ""))) if letter.isalpha()})
    return {
        "COG": cog,
        "KO": sorted(set(_terms(row.get("kos", "")))),
        "EC": sorted(set(_terms(row.get("ecs", "")))),
        "GO": sorted(set(_terms(row.get("gos", "")))),
    }


def read_functional_table(path: Path) -> dict[str, dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return {row["query"]: row for row in csv.DictReader(handle, delimiter="\t")}


def aggregate_functional_abundance(
    gene_rows: list[dict],
    functional_rows: dict[str, dict[str, str]],
    output_path: Path,
    *,
    term_attribution: str = "full",
) -> dict:
    """Sum gene counts and TPM per functional term and write a long-form TSV."""
    if term_attribution not in ("full", "split"):
        raise AnnotationError("abundance.term_attribution must be 'full' or 'split'")
    totals: dict[str, dict[str, list[float]]] = {
        namespace: defaultdict(lambda: [0, 0.0, 0.0]) for namespace in FUNCTIONAL_NAMESPACES
    }
    covered_tpm = Counter()
    annotated_tpm = 0.0
    for row in gene_rows:
        functional = functional_rows.get(row["gene_id"])
        if functional is None:
            continue
        if functional.get("annotation_status") == "annotated":
            annotated_tpm += row["tpm"]
        for namespace, terms in _gene_terms(functional).items():
            if not terms:
                continue
            covered_tpm[namespace] += row["tpm"]
            weight = 1.0 if term_attribution == "full" else 1.0 / len(terms)
            for term in terms:
                entry = totals[namespace][term]
                entry[0] += 1
                entry[1] += row["count"] * weight
                entry[2] += row["tpm"] * weight

    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("namespace", "term", "gene_count", "count", "tpm"))
        for namespace in FUNCTIONAL_NAMESPACES:
            for term, (genes, count, tpm) in sorted(totals[namespace].items()):
                writer.writerow((namespace, term, genes, f"{count:.4f}", f"{tpm:.6f}"))
    return {
        "term_attribution": term_attribution,
        "tpm_in_annotated_genes": round(annotated_tpm, 6),
        "tpm_with_term": {ns: round(covered_tpm[ns], 6) for ns in FUNCTIONAL_NAMESPACES},
    }


def _bbmap_version() -> str:
    try:
        result = subprocess.run(
            ["bbmap.sh", "--version"], capture_output=True, text=True, check=False, timeout=60
        )
    except FileNotFoundError as exc:
        raise AnnotationError(
            "BBMap is required for gene abundance; install bbmap>=39 or use --skip-abundance"
        ) from exc
    except subprocess.TimeoutExpired:
        return "unknown"
    return _parse_bbmap_version(f"{result.stdout}\n{result.stderr}")


_BBMAP_VERSION = re.compile(r"\bversion\s+(\d+(?:\.\d+)+)", re.IGNORECASE)


def _parse_bbmap_version(text: str) -> str:
    """Extract the BBTools release from `bbmap.sh --version` or a BBMap log."""
    match = _BBMAP_VERSION.search(text)
    return match.group(1) if match else "unknown"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _bbmap_command(
    reads: list[Path],
    contigs: Path,
    *,
    interleaved: bool,
    threads: int,
    min_identity: float,
    ambiguous: str,
    java_memory: str | None,
) -> list[str]:
    command = ["bbmap.sh"]
    if java_memory:
        command.append(f"-Xmx{java_memory}")
    command += [
        f"ref={contigs}",
        f"in={reads[0]}",
        "out=stdout.sam",
        "nodisk=t",
        "sam=1.4",
        "secondary=f",
        f"ambiguous={ambiguous}",
        f"minid={min_identity}",
        f"threads={threads}",
        "overwrite=t",
    ]
    if len(reads) == 2:
        command.append(f"in2={reads[1]}")
    else:
        command.append(f"interleaved={'t' if interleaved else 'f'}")
    return command


def _stream_bbmap(command: list[str], log_path: Path, cwd: Path) -> Iterator[str]:
    with log_path.open("w", encoding="utf-8") as log:
        log.write("[MetaQuest] " + " ".join(command) + "\n")
        log.flush()
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=log, text=True, cwd=cwd, bufsize=1 << 20
        )
        try:
            assert process.stdout is not None
            yield from process.stdout
        finally:
            if process.stdout:
                process.stdout.close()
            returncode = process.wait()
    if returncode:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-20:])
        raise AnnotationError(
            f"BBMap failed with exit code {returncode}. See {log_path}. Last output:\n{tail}"
        )


def _write_gene_table(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=GENE_ABUNDANCE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "count": f"{row['count']:.4f}",
                    "mean_depth": f"{row['mean_depth']:.4f}",
                    "rpk": f"{row['rpk']:.6f}",
                    "tpm": f"{row['tpm']:.6f}",
                }
            )


def _read_gene_table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = []
        for row in csv.DictReader(handle, delimiter="\t"):
            for key in ("count", "mean_depth", "rpk", "tpm"):
                row[key] = float(row[key])
            rows.append(row)
    return rows


def run_gene_abundance(
    reads: list[Path],
    contigs_fasta: Path,
    genes_gff: Path,
    output_dir: Path,
    *,
    functional_table: Path | None = None,
    read_mode: str = "paired",
    threads: int = 8,
    min_identity: float = 0.95,
    ambiguous: str = "random",
    term_attribution: str = "full",
    java_memory: str | None = None,
    reuse: bool = False,
) -> tuple[Path, Path | None, dict, bool]:
    """Map reads to contigs and write gene and functional abundance tables."""
    reads = [Path(path).resolve() for path in reads]
    contigs_fasta = Path(contigs_fasta).resolve()
    genes_gff = Path(genes_gff).resolve()
    output_dir = (Path(output_dir) / "gene_abundance").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if not 0 < min_identity <= 1:
        raise AnnotationError("abundance.min_identity must be in (0, 1]")
    if ambiguous not in ("random", "best", "all", "toss"):
        raise AnnotationError("abundance.ambiguous must be random, best, all, or toss")
    for path in (*reads, contigs_fasta, genes_gff):
        if not path.is_file():
            raise AnnotationError(f"Gene abundance input not found: {path}")

    gene_table = output_dir / "gene_abundance.tsv"
    summary_path = output_dir / "summary.json"
    completion_path = output_dir / "completion.json"
    paired = len(reads) == 2 or read_mode == "interleaved"
    version = _bbmap_version()
    expected_state = {
        "contigs_sha256": _sha256(contigs_fasta),
        "genes_gff_sha256": _sha256(genes_gff),
        "reads": [{"path": str(path), "size_bytes": path.stat().st_size} for path in reads],
        "bbmap_version": version,
        "min_identity": min_identity,
        "ambiguous": ambiguous,
        "paired": paired,
    }

    reused = False
    if reuse and completion_path.is_file() and gene_table.is_file() and summary_path.is_file():
        existing = json.loads(completion_path.read_text(encoding="utf-8"))
        if all(existing.get(key) == value for key, value in expected_state.items()):
            rows = _read_gene_table(gene_table)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            reused = True

    if not reused:
        genes = read_gene_gff(genes_gff)
        counter = AbundanceCounter(GeneIndex(genes), min_identity=min_identity, paired=paired)
        if genes:
            command = _bbmap_command(
                reads,
                contigs_fasta,
                interleaved=read_mode == "interleaved" and len(reads) == 1,
                threads=threads,
                min_identity=min_identity,
                ambiguous=ambiguous,
                java_memory=java_memory,
            )
            counter.add_sam_lines(_stream_bbmap(command, output_dir / "bbmap.log", output_dir))
        # Cache keys keep the --version result; the log fills in provenance.
        reported_version = version
        if version == "unknown" and (output_dir / "bbmap.log").is_file():
            reported_version = _parse_bbmap_version(
                (output_dir / "bbmap.log").read_text(encoding="utf-8", errors="replace")
            )
        tpm = compute_tpm(genes, counter.counts)
        rows = [
            {
                "gene_id": gene.gene_id,
                "contig_id": gene.contig_id,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "partial": gene.partial,
                "length_bp": gene.length,
                "count": counter.counts.get(gene.gene_id, 0.0),
                "mean_depth": counter.aligned_bases.get(gene.gene_id, 0) / gene.length,
                "rpk": tpm[gene.gene_id][0],
                "tpm": tpm[gene.gene_id][1],
                "annotation_status": "not_run",
            }
            for gene in genes
        ]
        stats = counter.stats
        unit = "fragments" if paired else "reads"
        fragments = stats["fragments"]
        summary = {
            "tool": "BBMap",
            "tool_version": reported_version,
            "count_unit": unit,
            "assignment": "alignment midpoint within CDS; fragment split across distinct genes",
            "normalization": "TPM from CDS nucleotide length",
            "min_identity": min_identity,
            "ambiguous": ambiguous,
            "threads": threads,
            "total_genes": len(genes),
            "genes_with_counts": sum(1 for row in rows if row["count"] > 0),
            "input_reads": stats["reads"],
            "mapped_reads": stats["reads_mapped"],
            "reads_below_identity": stats["reads_below_identity"],
            f"input_{unit}": fragments,
            f"mapped_{unit}": stats["fragments_mapped"],
            f"{unit}_in_genes": stats["fragments_in_genes"],
            f"intergenic_{unit}": stats["fragments_intergenic"],
            "mapping_rate": stats["fragments_mapped"] / fragments if fragments else 0.0,
            "gene_assignment_rate": stats["fragments_in_genes"] / fragments if fragments else 0.0,
            "reused": False,
        }

    functional_path = None
    if functional_table is not None and Path(functional_table).is_file():
        functional_rows = read_functional_table(Path(functional_table))
        for row in rows:
            row["annotation_status"] = functional_rows.get(row["gene_id"], {}).get(
                "annotation_status", "unannotated"
            )
        functional_path = output_dir / "functional_abundance.tsv"
        summary["functional"] = aggregate_functional_abundance(
            rows, functional_rows, functional_path, term_attribution=term_attribution
        )
    else:
        summary.pop("functional", None)
    summary["reused"] = reused

    _write_gene_table(gene_table, rows)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    completion_path.write_text(
        json.dumps({**expected_state, "completed": True}, indent=2) + "\n", encoding="utf-8"
    )
    return gene_table, functional_path, summary, reused
