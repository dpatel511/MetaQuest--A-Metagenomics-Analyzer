# MetaQuest

MetaQuest is a research-use command-line pipeline for short-read metagenomic
FASTQ analysis. Version `2.0.0a1` is a stabilization release: the
maintained runtime is intentionally smaller than earlier prototypes and does
not make clinical or pathogen-risk claims.

## Current capabilities

- full FASTQ structural validation and sampled quality warnings
- fastp adapter and quality preprocessing with retained-read QC
- synchronized paired-read identifier and count validation
- Kraken2 taxonomic classification
- Bracken abundance re-estimation
- optional taxonomy-only execution
- MEGAHIT metagenomic assembly
- Pyrodigal metagenomic gene prediction
- eggNOG-mapper orthology-based functional annotation
- per-gene annotations and aggregated COG, KO, EC, and GO counts
- read-mapping gene abundance (BBMap, TPM) with TPM-weighted COG, KO, EC, and
  GO functional abundance
- descriptive text, JSON, offline HTML, and publication figure reporting
- explicit classified/unclassified denominators and reproducibility metadata
- versioned taxonomy and eggNOG database installation

MetaQuest does not currently perform read trimming, host-read removal, validated
AMR or virulence analysis, clinical diagnosis, or treatment recommendations.

## Pipeline

~~~text
FASTQ
  ├── validation
  ├── fastp preprocessing
  ├── Kraken2 → Bracken
  ├── MEGAHIT
  ├── Pyrodigal
  ├── eggNOG-mapper → COG / KO / EC / GO
  ├── BBMap read mapping → gene TPM → functional abundance
  └── descriptive reporting → HTML + figures + plotted-data TSV
~~~

The maintained pipeline currently runs directly in Python. Migration to
Snakemake and optional Bowtie2 host filtering remain planned.

## Quick start

Install the Python distribution from PyPI:

~~~bash
python -m pip install metaquest-bio
metaquest --version
~~~

The distribution is named `metaquest-bio`; the installed Python package and
command remain `metaquest`. The PyPI distribution does not bundle Kraken2,
Bracken, MEGAHIT, DIAMOND, or eggNOG-mapper. For the complete runtime, clone
the repository and create its Conda environment before installing the package:

~~~bash
conda env create -f environment/environment.yml
conda activate metaquest
python -m pip install .
~~~

Inspect and install the taxonomy database:

~~~bash
metaquest databases --list
metaquest databases --database taxonomy
metaquest check
~~~

By default, MetaQuest installs reference data in `./databases`. Use
`--db-dir` or `METAQUEST_DB_DIR` only when you intentionally want another
location.

Run taxonomic profiling:

~~~bash
metaquest run \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --output results/sample
~~~

Use `--skip-functional` to stop after Pyrodigal, `--skip-abundance` to skip
read mapping and gene abundance, or `--taxonomy-only` to skip assembly, gene
prediction, functional annotation, and abundance. The former
`--skip-annotation` spelling remains as a deprecated alias for
`--taxonomy-only`.

Each fresh run requires a new or empty output directory. Use `--resume` only
to reuse a matching MetaQuest run, or `--force` to move an existing output
directory to a timestamped backup before starting fresh.

## Commands

~~~text
metaquest check       # verify tools and databases
metaquest validate    # validate FASTQ input
metaquest run         # run the analysis pipeline
metaquest databases   # inspect or install reference data
metaquest init-config # create a configuration file
~~~

The previous `analyze` and `setup-db` names remain available as aliases for
backward compatibility.

Global output controls may appear before or after the command:

~~~bash
metaquest run --verbose --single reads.fastq.gz --output results/
metaquest databases --no-color --list
metaquest --quiet check
~~~

For systems that cannot load the Kraken2 database fully into RAM, enable
Kraken2 memory mapping with `--low-memory`:

~~~bash
metaquest run --low-memory --single reads.fastq.gz --output results/
~~~

Currently this flag only adds Kraken2's `--memory-mapping` option. It does not
change resource settings for Bracken, MEGAHIT, or gene prediction.

## Database storage

Reference databases are not stored in Git. From the repository root, install
them directly into the default `databases/` directory:

~~~bash
metaquest databases --database taxonomy
metaquest databases --database functional
metaquest check
~~~

The database path resolution order is:

1. `--db-dir`
2. `METAQUEST_DB_DIR`
3. `databases.base_dir` in YAML
4. `./databases`

## Documentation

- [Installation](docs/installation.md)
- [Database management](docs/databases.md)
- [Usage](docs/usage.md)
- [Gene prediction and functional annotation](docs/annotation.md)
- [Release process](docs/releasing.md)
- [Changelog](CHANGELOG.md)

## Development status

MetaQuest is alpha software intended for reproducible method development. Before
a stable release, the workflow requires end-to-end Snakemake execution,
functional and pathogen-associated method validation, resource benchmarks,
comparison-method validation, and publication datasets.

## License

MetaQuest is distributed under the GNU General Public License v3.0 or later.
Copyright (c) 2026 Dev Patel. External programs and reference databases are not
included in the Python distribution and remain subject to their own terms. See
[Third-party licenses](THIRD_PARTY_LICENSES.md) for the maintained dependency
and provenance inventory.
