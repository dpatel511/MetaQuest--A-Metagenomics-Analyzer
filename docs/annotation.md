# Gene prediction and functional annotation

This document describes the maintained annotation boundary in
MetaQuest `2.0.0a1`.

## Current stage

~~~text
MEGAHIT contigs
  ↓
Pyrodigal metagenomic gene prediction
  ├── genes.faa
  ├── genes.fna
  ├── genes.gff3
  ├── contig_id_map.tsv
  └── summary.json
  ↓
eggNOG-mapper 2.1.15 + eggNOG 5.0.2
  ├── functional_annotations.tsv
  ├── functional_category_summary.tsv
  ├── metaquest.emapper.annotations
  ├── summary.json
  └── completion.json
~~~

## Pyrodigal

MetaQuest initializes `pyrodigal.GeneFinder(meta=True)` and processes one
contig at a time. This bounds memory use and applies Pyrodigal's metagenomic
models independently to assembled contigs.

Contigs shorter than `annotation.min_contig_length` are skipped. The default is
200 bp.

### Outputs

| File | Contents |
|---|---|
| `gene_prediction/genes.faa` | predicted protein sequences without terminal stop characters |
| `gene_prediction/genes.fna` | predicted coding nucleotide sequences |
| `gene_prediction/genes.gff3` | CDS coordinates and translation-table metadata |
| `gene_prediction/summary.json` | tool version, mode, contig counts, and gene count |

Stable contig identifiers derive from sequence hashes, making equivalent
assemblies independent of contig ordering. `contig_id_map.tsv` maps each
stable identifier back to the original assembler identifier and full sequence
checksum.

## eggNOG functional annotation

Every protein predicted by Pyrodigal is submitted to eggNOG-mapper using the
DIAMOND method and automatic taxonomic scope. MetaQuest left-joins mapper
results back to the complete gene catalog, so genes without an assignment
remain explicit `unannotated` rows rather than disappearing.

DIAMOND uses a `0.5` block size by default to keep the maintained workflow
practical on the supported 8-16 GB memory target. This value is recorded in
the functional summary and can be changed in YAML as
`annotation.diamond_block_size`.

| File | Contents |
|---|---|
| `functional_annotations.tsv` | one row per predicted gene, including unannotated genes |
| `functional_category_summary.tsv` | gene counts aggregated by COG, KO, EC, and GO term |
| `metaquest.emapper.annotations` | unmodified eggNOG-mapper annotation output |
| `eggnog_mapper.log` | mapper stdout and stderr for diagnostics |
| `summary.json` | tool/database versions, parameters, and annotation counts |
| `completion.json` | input checksum and parameters used for restart-safe reuse |

If the canonical protein-content checksum, mapper version, database release,
taxonomic scope, and E-value match a completed run, MetaQuest reuses it.
Record ordering does not invalidate the cache. Mapper temporary files remain
inside the functional output directory and are removed after success or
failure. Python output is unbuffered so `eggnog_mapper.log` can be monitored.

The files above describe gene presence. The gene abundance stage below adds
read-based quantification.

## Gene and functional abundance

~~~text
fastp-cleaned reads + contigs.stable.fasta + genes.gff3
  ↓
BBMap (nodisk, primary alignments, extended CIGAR)
  ↓
gene_abundance/
  ├── gene_abundance.tsv
  ├── functional_abundance.tsv
  ├── bbmap.log
  ├── summary.json
  └── completion.json
~~~

Reads are mapped to the stable-ID contigs produced by gene prediction, so
alignment coordinates match `genes.gff3` directly. SAM output is streamed into
Python and is not written to disk.

Counting rules:

- Only primary alignments are used; secondary and supplementary records are
  ignored. Multi-mapping reads are placed by BBMap's `ambiguous` policy
  (default `random`).
- Alignment identity is computed from the extended CIGAR (`=`/`X`/`I`/`D`) and
  must reach `abundance.min_identity` (default `0.95`).
- A read is assigned to the CDS containing its alignment midpoint.
- For paired data, mates are combined so each fragment counts once. When the
  two mates fall in different genes, the fragment is split equally between
  them. Counts are therefore reported in fragments for paired input and in
  reads for single-end input.
- Mean depth is aligned reference bases overlapping the CDS divided by CDS
  length.

Gene values are length-normalized as reads per kilobase (RPK) and converted to
TPM, so gene TPM sums to one million across all predicted genes.

| File | Contents |
|---|---|
| `gene_abundance.tsv` | one row per predicted gene: coordinates, partial flag, count, mean depth, RPK, TPM, annotation status |
| `functional_abundance.tsv` | per COG, KO, EC, and GO term: contributing genes, summed count, summed TPM |
| `bbmap.log` | BBMap command and statistics |
| `summary.json` | mapper version, parameters, mapping and gene-assignment rates |
| `completion.json` | contig, GFF3, read, and parameter state for restart-safe reuse |

With `abundance.term_attribution: full` (default), a gene annotated with
several terms contributes its full value to each, so a namespace can sum to
more than one million TPM. With `split`, the value is divided equally among
the gene's terms and the namespace total equals the TPM of genes carrying at
least one term. `summary.json` records the TPM share in eggNOG-annotated genes
and in genes carrying each namespace.

Use `--skip-abundance` to disable this stage. BBMap's Java heap can be capped
with `abundance.java_memory` (for example `"6g"`).

## Interpretation limits

- A sequence-similarity hit is not experimental confirmation of function.
- Presence of a gene is not proof of expression or phenotype.
- Assembly can merge, fragment, or omit low-abundance sequences.
- Closely related proteins can have different substrate specificity.
- Functional abundance depends on read mapping and normalization choices.
- Abundance covers only reads that map to assembled, predicted genes; the
  mapping and gene-assignment rates state how much of the sample that is.
- TPM is compositional. Cross-sample comparisons require compositional or
  count-based statistical methods rather than direct TPM ratios.
- Gene counts measure DNA copies in the sample, not expression.
