import csv
import json
import random
import shutil
from collections import Counter

import pytest

from metaquest.core import abundance
from metaquest.core.abundance import (
    AbundanceCounter,
    _parse_bbmap_version,
    Gene,
    GeneIndex,
    _alignment_metrics,
    aggregate_functional_abundance,
    compute_tpm,
    read_gene_gff,
    run_gene_abundance,
)
from metaquest.exceptions import AnnotationError
from metaquest.settings import load_config, validate_config


GENES = [
    Gene("c1_1", "c1", 101, 400, "+", "00"),
    Gene("c1_2", "c1", 501, 1100, "-", "00"),
    Gene("c2_1", "c2", 1, 300, "+", "10"),
]


def _sam(name, flag, contig, pos, cigar="100=", extra=""):
    fields = [name, str(flag), contig, str(pos), "60", cigar, "*", "0", "0", "*", "*"]
    return "\t".join(fields) + (f"\t{extra}" if extra else "")


def test_gene_index_finds_overlaps_and_midpoints():
    index = GeneIndex(GENES)
    assert [gene.gene_id for gene in index.containing("c1", 250)] == ["c1_1"]
    assert index.containing("c1", 450) == []
    assert sorted(g.gene_id for g in index.overlapping("c1", 350, 550)) == ["c1_1", "c1_2"]
    assert index.containing("missing", 10) == []


def test_alignment_identity_uses_extended_cigar_and_nm_fallback():
    assert _alignment_metrics("90=10X", []) == (100, 0.9)
    assert _alignment_metrics("50=2I48=", []) == (98, 98 / 100)
    assert _alignment_metrics("10S90M", ["NM:i:9"]) == (90, 0.9)
    assert _alignment_metrics("100M", []) == (100, None)


def test_single_end_counts_by_midpoint_and_reports_intergenic():
    counter = AbundanceCounter(GeneIndex(GENES), min_identity=0.95, paired=False)
    counter.add_sam_lines(
        [
            "@HD\tVN:1.4",
            _sam("r1", 0, "c1", 101),  # midpoint 150 -> c1_1
            _sam("r2", 16, "c1", 401),  # midpoint 450 -> intergenic
            _sam("r3", 0, "c1", 501, "80=20X"),  # identity 0.8 -> rejected
            _sam("r4", 4, "*", 0, "*"),  # unmapped
            _sam("r5", 256, "c1", 101),  # secondary -> ignored
        ]
    )
    assert counter.counts == Counter({"c1_1": 1.0})
    assert counter.stats["reads"] == 4
    assert counter.stats["fragments_mapped"] == 2
    assert counter.stats["fragments_intergenic"] == 1
    assert counter.stats["reads_below_identity"] == 1
    # r2 spans 401-500, overlapping no CDS base
    assert counter.aligned_bases == Counter({"c1_1": 100})


def test_paired_fragment_counts_once_and_splits_across_genes():
    counter = AbundanceCounter(GeneIndex(GENES), min_identity=0.95, paired=True)
    counter.add_sam_lines(
        [
            _sam("f1", 99, "c1", 151),  # both mates in c1_1 -> 1 count
            _sam("f1", 147, "c1", 251),
            _sam("f2", 99, "c1", 201),  # mates in c1_1 and c1_2 -> 0.5 each
            _sam("f2", 147, "c1", 701),
            _sam("f3", 73, "c2", 51),  # mate unmapped -> still one fragment
            _sam("f3", 133, "*", 0, "*"),
        ]
    )
    assert counter.counts == Counter({"c1_1": 1.5, "c1_2": 0.5, "c2_1": 1.0})
    assert counter.stats["fragments"] == 3
    assert counter.stats["fragments_in_genes"] == 3


def test_tpm_normalizes_by_length_and_sums_to_one_million():
    tpm = compute_tpm(GENES, Counter({"c1_1": 3, "c1_2": 6}))
    # c1_1: 3 / 0.3 kb = 10 RPK ; c1_2: 6 / 0.6 kb = 10 RPK -> equal TPM
    assert tpm["c1_1"][1] == pytest.approx(500_000)
    assert tpm["c1_2"][1] == pytest.approx(500_000)
    assert tpm["c2_1"] == (0.0, 0.0)
    assert sum(value for _, value in tpm.values()) == pytest.approx(1_000_000)
    assert compute_tpm(GENES, Counter())["c1_1"] == (0.0, 0.0)


def _functional_rows():
    return {
        "c1_1": {"query": "c1_1", "annotation_status": "annotated", "cog_categories": "KL",
                 "kos": "ko:K00001,ko:K00002", "ecs": "1.1.1.1", "gos": "-"},
        "c1_2": {"query": "c1_2", "annotation_status": "annotated", "cog_categories": "K",
                 "kos": "ko:K00001", "ecs": "", "gos": ""},
        "c2_1": {"query": "c2_1", "annotation_status": "unannotated", "cog_categories": "",
                 "kos": "", "ecs": "", "gos": ""},
    }


def _gene_rows():
    return [
        {"gene_id": "c1_1", "count": 2.0, "tpm": 600_000.0},
        {"gene_id": "c1_2", "count": 1.0, "tpm": 300_000.0},
        {"gene_id": "c2_1", "count": 1.0, "tpm": 100_000.0},
    ]


def _read_terms(path):
    with path.open() as handle:
        return {(r["namespace"], r["term"]): r for r in csv.DictReader(handle, delimiter="\t")}


def test_functional_abundance_full_attribution(tmp_path):
    output = tmp_path / "functional_abundance.tsv"
    info = aggregate_functional_abundance(_gene_rows(), _functional_rows(), output)
    terms = _read_terms(output)
    assert float(terms[("KO", "ko:K00001")]["tpm"]) == pytest.approx(900_000)
    assert float(terms[("KO", "ko:K00002")]["tpm"]) == pytest.approx(600_000)
    assert terms[("KO", "ko:K00001")]["gene_count"] == "2"
    assert float(terms[("COG", "K")]["tpm"]) == pytest.approx(900_000)
    assert float(terms[("COG", "L")]["tpm"]) == pytest.approx(600_000)
    assert info["tpm_in_annotated_genes"] == pytest.approx(900_000)
    assert info["tpm_with_term"]["GO"] == 0


def test_functional_abundance_split_attribution_preserves_mass(tmp_path):
    output = tmp_path / "functional_abundance.tsv"
    aggregate_functional_abundance(
        _gene_rows(), _functional_rows(), output, term_attribution="split"
    )
    terms = _read_terms(output)
    assert float(terms[("KO", "ko:K00001")]["tpm"]) == pytest.approx(600_000)
    ko_total = sum(float(r["tpm"]) for (ns, _), r in terms.items() if ns == "KO")
    assert ko_total == pytest.approx(900_000)
    with pytest.raises(AnnotationError):
        aggregate_functional_abundance([], {}, output, term_attribution="bogus")


def test_read_gene_gff_parses_pyrodigal_features(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "# Sequence Data: seqnum=1;seqlen=2000\n"
        "c1\tpyrodigal_v3.7.1\tCDS\t101\t400\t9.1\t+\t0\tID=c1_1;partial=00;conf=99\n",
        encoding="utf-8",
    )
    assert read_gene_gff(gff) == [Gene("c1_1", "c1", 101, 400, "+", "00")]


def test_abundance_config_defaults_validate():
    config = load_config()
    assert config.abundance.min_identity == 0.95
    assert config.abundance.term_attribution == "full"
    assert validate_config(config) == (True, [])


def test_run_gene_abundance_with_mocked_mapper_and_reuse(tmp_path, monkeypatch):
    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(">c1\n" + "A" * 1200 + "\n>c2\n" + "C" * 300 + "\n")
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "".join(
            f"{g.contig_id}\tpyrodigal\tCDS\t{g.start}\t{g.end}\t.\t{g.strand}\t0\t"
            f"ID={g.gene_id};partial={g.partial}\n"
            for g in GENES
        )
    )
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n")
    functional = tmp_path / "functional_annotations.tsv"
    with functional.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(_functional_rows()["c1_1"]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(_functional_rows().values())

    calls = []

    def fake_stream(command, log_path, cwd):
        calls.append(command)
        yield _sam("r1", 0, "c1", 101)
        yield _sam("r2", 0, "c1", 601)
        yield _sam("r3", 4, "*", 0, "*")

    monkeypatch.setattr(abundance, "_bbmap_version", lambda: "40.02")
    monkeypatch.setattr(abundance, "_stream_bbmap", fake_stream)

    table, functional_out, summary, reused = run_gene_abundance(
        [reads], contigs, gff, tmp_path, functional_table=functional, read_mode="single"
    )
    assert not reused and len(calls) == 1
    assert "interleaved=f" in calls[0] and "sam=1.4" in calls[0]
    assert summary["count_unit"] == "reads"
    assert summary["mapping_rate"] == pytest.approx(2 / 3)
    with table.open() as handle:
        rows = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert rows["c2_1"]["count"] == "0.0000"
    assert rows["c2_1"]["partial"] == "10"
    assert rows["c1_1"]["annotation_status"] == "annotated"
    # equal counts, c1_1 is half the length of c1_2 -> twice the TPM
    assert float(rows["c1_1"]["tpm"]) == pytest.approx(2 * float(rows["c1_2"]["tpm"]))
    assert functional_out.is_file()

    _, _, summary, reused = run_gene_abundance(
        [reads], contigs, gff, tmp_path, functional_table=functional,
        read_mode="single", reuse=True,
    )
    assert reused and len(calls) == 1
    assert summary["reused"] is True
    completion = json.loads((tmp_path / "gene_abundance" / "completion.json").read_text())
    assert completion["completed"] is True


@pytest.mark.skipif(shutil.which("bbmap.sh") is None, reason="BBMap not installed")
def test_bbmap_end_to_end_recovers_relative_abundance(tmp_path):
    rng = random.Random(7)
    high = "".join(rng.choice("ACGT") for _ in range(2000))
    low = "".join(rng.choice("ACGT") for _ in range(2000))
    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(f">high\n{high}\n>low\n{low}\n")
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "high\tpyrodigal\tCDS\t201\t1800\t.\t+\t0\tID=high_1;partial=00\n"
        "low\tpyrodigal\tCDS\t201\t1800\t.\t+\t0\tID=low_1;partial=00\n"
    )
    complement = str.maketrans("ACGT", "TGCA")
    r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
    with r1.open("w") as out1, r2.open("w") as out2:
        for index, (seq, copies) in enumerate([(high, 300), (low, 100)]):
            for copy in range(copies):
                start = rng.randint(200, 1500)
                forward = seq[start:start + 100]
                reverse = seq[start + 200:start + 300].translate(complement)[::-1]
                name = f"f{index}_{copy}"
                out1.write(f"@{name}/1\n{forward}\n+\n{'I' * 100}\n")
                out2.write(f"@{name}/2\n{reverse}\n+\n{'I' * 100}\n")

    table, _, summary, _ = run_gene_abundance(
        [r1, r2], contigs, gff, tmp_path, read_mode="paired", threads=2, java_memory="1g"
    )
    with table.open() as handle:
        rows = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert summary["count_unit"] == "fragments"
    assert summary["input_fragments"] == 400
    assert summary["mapping_rate"] > 0.99
    assert float(rows["high_1"]["count"]) == pytest.approx(300, abs=2)
    assert float(rows["low_1"]["count"]) == pytest.approx(100, abs=2)
    assert float(rows["high_1"]["tpm"]) / float(rows["low_1"]["tpm"]) == pytest.approx(3, rel=0.05)


def test_abundance_report_html_and_plots(tmp_path):
    from types import SimpleNamespace

    from metaquest.pipeline.context import AbundanceResult
    from metaquest.reporting.stable_reporter import generate_stable_reports

    directory = tmp_path / "gene_abundance"
    directory.mkdir()
    gene_table = directory / "gene_abundance.tsv"
    gene_table.write_text("gene_id\n")
    functional = directory / "functional_abundance.tsv"
    info = aggregate_functional_abundance(_gene_rows(), _functional_rows(), functional)
    summary = {
        "tool_version": "40.02", "count_unit": "fragments", "min_identity": 0.95,
        "ambiguous": "random", "input_fragments": 10, "mapped_fragments": 8,
        "fragments_in_genes": 6, "intergenic_fragments": 2, "mapping_rate": 0.8,
        "gene_assignment_rate": 0.6, "genes_with_counts": 3, "total_genes": 3,
        "normalization": "TPM from CDS nucleotide length", "functional": info,
    }
    ctx = SimpleNamespace(
        output_dir=tmp_path, completed_stages=[], config=load_config(db_dir=tmp_path / "db"),
        read_mode="paired", classification=None, assembly=None, annotation=None,
        preprocessing=None,
        abundance=AbundanceResult(gene_table, functional, summary, False),
    )
    generate_stable_reports(ctx)
    result = json.loads((tmp_path / "analysis_summary.json").read_text())
    assert result["abundance"]["mapping_rate"] == 0.8
    assert result["abundance"]["functional_abundance"] == "gene_abundance/functional_abundance.tsv"
    assert "Gene abundance" in (tmp_path / "report.html").read_text()
    assert "Mapping rate: 80.00%" in (tmp_path / "04_abundance_report.txt").read_text()
    assert (tmp_path / "plots" / "abundance_ko.svg").is_file()
    assert (tmp_path / "plots" / "plot_data" / "abundance_cog.tsv").is_file()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("WARNING: incubator\nBBTools version 40.02\n", "40.02"),
        ("BBMap version 39.06\n", "39.06"),
        ("java -ea -Xmx1g ...\nExecuting align2.BBMap [build=1]\nVersion 39.01\n", "39.01"),
        ("no version here", "unknown"),
    ],
)
def test_parse_bbmap_version_formats(text, expected):
    assert _parse_bbmap_version(text) == expected
