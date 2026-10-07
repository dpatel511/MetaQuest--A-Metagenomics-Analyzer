"""
MetaQuest Settings
==================
YAML-based configuration system. Single source of truth for all pipeline parameters.

Resolution order for database path:
  1. METAQUEST_DB_DIR environment variable
  2. databases.base_dir in config YAML
  3. ./databases (relative to working directory)
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
import yaml

from metaquest.exceptions import ConfigError


# ---------------------------------------------------------------------------
# Config dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatabasesConfig:
    base_dir: Path
    kraken_db: Path
    functional_dir: Path
    kraken_db_version: str = "Standard-8 (2026-06-26)"


@dataclass(frozen=True)
class AssemblyConfig:
    assembler: str = "megahit"
    threads: int = 8


@dataclass(frozen=True)
class ClassificationConfig:
    threads: int = 8
    taxonomic_level: str = "S"
    min_hit_groups: int = 2
    bracken_threshold: int = 10


@dataclass(frozen=True)
class AnnotationConfig:
    tool: str = "pyrodigal"
    functional_tool: str = "eggnog-mapper"
    threads: int = 8
    evalue: float = 1e-6
    diamond_block_size: float = 0.5
    tax_scope: str = "auto"
    eggnog_version: str = "2.1.15"
    eggnog_database_release: str = "5.0.2"
    min_contig_length: int = 200


@dataclass(frozen=True)
class AbundanceConfig:
    mapper: str = "bbmap"
    threads: int = 8
    min_identity: float = 0.95
    ambiguous: str = "random"
    term_attribution: str = "full"
    java_memory: str | None = None


@dataclass(frozen=True)
class PreprocessingConfig:
    enabled: bool = True
    threads: int = 4
    qualified_quality_phred: int = 20
    length_required: int = 50


@dataclass(frozen=True)
class ReportingConfig:
    top_taxa: int = 20
    top_functional_terms: int = 20
    plot_formats: tuple[str, ...] = ("svg", "png")
    plot_dpi: int = 300
    color_palette: str = "colorblind"


@dataclass(frozen=True)
class MetaQuestConfig:
    databases: DatabasesConfig
    assembly: AssemblyConfig = AssemblyConfig()
    classification: ClassificationConfig = ClassificationConfig()
    annotation: AnnotationConfig = AnnotationConfig()
    abundance: AbundanceConfig = AbundanceConfig()
    preprocessing: PreprocessingConfig = PreprocessingConfig()
    reporting: ReportingConfig = ReportingConfig()


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG_PATH = Path(__file__).parent / "metaquest_default.yaml"

_cached_config: MetaQuestConfig | None = None


def _resolve_db_dir(raw: dict, override: Path | None = None) -> Path:
    """Resolve database directory from env > yaml > default."""
    if override is not None:
        return Path(override)
    env_val = os.environ.get("METAQUEST_DB_DIR")
    if env_val:
        return Path(env_val)
    yaml_val = raw.get("databases", {}).get("base_dir")
    if yaml_val:
        return Path(yaml_val)
    return Path("./databases")


def _build_databases_config(raw: dict, base_dir: Path) -> DatabasesConfig:
    db = raw.get("databases", {})
    kraken = Path(db.get("kraken_db") or str(base_dir / "taxonomy"))
    functional = Path(db.get("functional_dir") or str(base_dir / "functional"))
    version = db.get("kraken_db_version", "Standard-8 (2026-06-26)")
    return DatabasesConfig(
        base_dir=base_dir,
        kraken_db=kraken,
        functional_dir=functional,
        kraken_db_version=version,
    )


def _build_section(cls, raw: dict, key: str):
    """Build a frozen dataclass from a YAML section, ignoring unknown keys."""
    section = raw.get(key, {}) or {}
    valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
    filtered = {k: v for k, v in section.items() if k in valid_fields}
    if cls is ReportingConfig and "plot_formats" in filtered:
        filtered["plot_formats"] = tuple(filtered["plot_formats"])
    # Convert Path fields
    for fname, fobj in cls.__dataclass_fields__.items():
        if fname in filtered and filtered[fname] is not None:
            if fobj.type in ("Path", "Path | None"):
                filtered[fname] = Path(filtered[fname])
    return cls(**filtered)


def load_config(
    config_path: Path | None = None,
    *,
    db_dir: Path | None = None,
) -> MetaQuestConfig:
    """
    Load configuration from YAML file.

    Resolution:
      1. User-provided config_path
      2. Default config bundled with package

    Values in user config override defaults (shallow merge per section).
    """
    global _cached_config

    # Load defaults
    if not _DEFAULT_CONFIG_PATH.exists():
        raise ConfigError(f"Default config missing: {_DEFAULT_CONFIG_PATH}")

    with open(_DEFAULT_CONFIG_PATH) as f:
        defaults = yaml.safe_load(f) or {}

    # Merge user config on top
    raw = dict(defaults)
    if config_path:
        if not config_path.exists():
            raise ConfigError(f"Config file not found: {config_path}")
        with open(config_path) as f:
            user = yaml.safe_load(f) or {}
        for key, val in user.items():
            if isinstance(val, dict) and key in raw and isinstance(raw[key], dict):
                raw[key] = {**raw[key], **val}
            else:
                raw[key] = val

    resolved_db_dir = _resolve_db_dir(raw, db_dir)
    databases = _build_databases_config(raw, resolved_db_dir)

    config = MetaQuestConfig(
        databases=databases,
        assembly=_build_section(AssemblyConfig, raw, "assembly"),
        classification=_build_section(ClassificationConfig, raw, "classification"),
        annotation=_build_section(AnnotationConfig, raw, "annotation"),
        abundance=_build_section(AbundanceConfig, raw, "abundance"),
        preprocessing=_build_section(PreprocessingConfig, raw, "preprocessing"),
        reporting=_build_section(ReportingConfig, raw, "reporting"),
    )

    _cached_config = config
    return config


def get_config() -> MetaQuestConfig:
    """Return cached config or load defaults."""
    global _cached_config
    if _cached_config is None:
        _cached_config = load_config()
    return _cached_config


def reset_config() -> None:
    """Clear cached config (for testing)."""
    global _cached_config
    _cached_config = None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_config(config: MetaQuestConfig | None = None) -> tuple[bool, list[str]]:
    """Validate configuration, return (is_valid, errors)."""
    cfg = config or get_config()
    errors: list[str] = []

    if cfg.assembly.threads < 1:
        errors.append("assembly.threads must be >= 1")
    if cfg.annotation.threads < 1:
        errors.append("annotation.threads must be >= 1")
    if cfg.annotation.diamond_block_size <= 0:
        errors.append("annotation.diamond_block_size must be > 0")
    if cfg.annotation.tax_scope != "auto":
        errors.append("annotation.tax_scope currently supports only 'auto'")
    if cfg.abundance.mapper != "bbmap":
        errors.append("abundance.mapper currently supports only 'bbmap'")
    if cfg.abundance.threads < 1:
        errors.append("abundance.threads must be >= 1")
    if not 0 < cfg.abundance.min_identity <= 1:
        errors.append("abundance.min_identity must be in (0, 1]")
    if cfg.abundance.ambiguous not in ("random", "best", "all", "toss"):
        errors.append("abundance.ambiguous must be one of: random, best, all, toss")
    if cfg.abundance.term_attribution not in ("full", "split"):
        errors.append("abundance.term_attribution must be 'full' or 'split'")
    if cfg.preprocessing.qualified_quality_phred < 0:
        errors.append("preprocessing.qualified_quality_phred must be >= 0")
    if cfg.preprocessing.length_required < 1:
        errors.append("preprocessing.length_required must be >= 1")
    if cfg.preprocessing.threads < 1:
        errors.append("preprocessing.threads must be >= 1")
    if cfg.reporting.top_taxa < 1 or cfg.reporting.top_functional_terms < 1:
        errors.append("reporting top limits must be >= 1")
    if cfg.reporting.plot_dpi < 72:
        errors.append("reporting.plot_dpi must be >= 72")
    unsupported = set(cfg.reporting.plot_formats) - {"svg", "png", "pdf"}
    if unsupported:
        errors.append(f"unsupported reporting.plot_formats: {', '.join(sorted(unsupported))}")
    return len(errors) == 0, errors
