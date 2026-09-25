# Third-party software and data

MetaQuest is licensed under GPL-3.0-or-later. This inventory documents the
maintained runtime dependencies for release engineering; it does not replace
the license text or terms supplied by each upstream project.

## Python distribution dependencies

These packages are installed as dependencies of the MetaQuest Python package.

| Package | Maintained constraint | Upstream license |
|---|---:|---|
| NumPy | `>=1.24` | BSD-3-Clause |
| pandas | `>=2.0` | BSD-3-Clause |
| Matplotlib | `>=3.8` | PSF-based |
| Seaborn | `>=0.13` | BSD-3-Clause |
| Biopython | `>=1.81` | Biopython License Agreement |
| Pyrodigal | `>=3.7.1,<4` | GPL-3.0-or-later |
| PyYAML | `>=6.0` | MIT |

Pyrodigal is imported directly by MetaQuest. Its GPL license is the reason the
MetaQuest distribution is provided under GPL-3.0-or-later instead of MIT.

## External executables

The following programs are installed separately, normally through Bioconda.
They are not copied into MetaQuest wheels or source distributions.

| Program | Maintained constraint | Upstream license | Project |
|---|---:|---|---|
| Kraken 2 | `>=2.1.3` | MIT | <https://github.com/DerrickWood/kraken2> |
| Bracken | `>=3.1` | GPL-3.0-or-later | <https://github.com/jenniferlu717/Bracken> |
| MEGAHIT | `>=1.2.9` | GPL-3.0 | <https://github.com/voutcn/megahit> |
| fastp | `>=0.23.4` | MIT | <https://github.com/OpenGene/fastp> |
| DIAMOND | `>=2.0.11,<2.1` | GPL-3.0 | <https://github.com/bbuchfink/diamond> |
| eggNOG-mapper | `2.1.15` | AGPL-3.0-only | <https://github.com/eggnogdb/eggnog-mapper> |
| BBTools/BBMap | `>=39` | BSD-3-Clause-LBNL | <https://sourceforge.net/projects/bbmap/> |

eggNOG-mapper is executed as a separate command. Version 2.1.15 is the final
v2 release compatible with the eggNOG v5.0.2 database. MetaQuest does not
vendor or modify eggNOG-mapper.

MetaQuest uses BBTools' `reformat.sh` when an interleaved FASTQ must be
split, and `bbmap.sh` to map reads to assembled contigs for gene abundance. The Bioconda package identifies its license as BSD-3-Clause-LBNL.

## Reference databases

Reference data is downloaded after installation and is never included in the
MetaQuest Python distribution.

| MetaQuest key | Pinned release | Upstream and provenance |
|---|---:|---|
| `taxonomy` | Standard-8, 2026-06-26 | Prebuilt Kraken 2 index from the Langmead Lab AWS index, derived from NCBI taxonomy and sequence resources |
| `functional` | eggNOG 5.0.2 | Files downloaded from the official eggNOG v5 server |

Software licenses do not automatically grant rights in third-party database
content. The downloader records URLs, versions, checksums, and provenance in
`metaquest-db.json`; users and redistributors remain responsible for reviewing
the current upstream data terms for their use case. MetaQuest redistributes
neither database.

## Release audit policy

- Update this inventory whenever a direct dependency, executable, or database
  source changes.
- Verify the exact upstream license for every pinned release.
- Keep external executables and database files out of wheels and source
  distributions.
- Treat an unknown or changed direct Python dependency license as a release
  blocker until it is reviewed.

This file records project due diligence and is not legal advice.
