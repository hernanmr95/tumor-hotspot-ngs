"""
pipeline/hotspot_matcher.py

Cruza las variantes de una muestra (ya parseadas desde el VCF) contra el
panel curado de hotspots (data/hotspot_panel.csv), y clasifica cada
variante encontrada según criterios explícitos y trazables:

  - MATCH: el gen y el cambio de proteína (HGVS_P) coinciden exactamente
    con una entrada del panel, y la VAF supera el umbral mínimo.
  - LOW_VAF_MATCH: coincide con el panel, pero la VAF está por debajo del
    umbral configurado — se reporta aparte en vez de descartarse en
    silencio, porque una VAF baja puede ser real (subclonal) o un
    artefacto, y esa decisión no debe tomarla el software solo.
  - NOT_IN_PANEL: la variante no coincide con ninguna entrada del panel
    (puede ser benigna, o un hotspot real no incluido en este panel de
    demostración — ver docs/scope_and_limitations.md).

El matching es por (gen, cambio de proteína), NO por coordenada genómica,
porque distintos anotadores pueden reportar la misma variante con
transcritos de referencia distintos y coordenadas ligeramente distintas;
el cambio de proteína en el transcrito canónico es más estable para
este propósito de demostración.
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from pipeline.vcf_parser import Variant

DEFAULT_MIN_VAF = 0.05


@dataclass
class MatchResult:
    variant: Variant
    status: str  # "MATCH", "LOW_VAF_MATCH", "NOT_IN_PANEL"
    panel_entry: Dict[str, str] = None


def load_hotspot_panel(panel_path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """Carga el panel CSV como dict {(gene, protein_change): fila_completa}."""
    panel = {}
    with open(panel_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["gene"].strip(), row["protein_change"].strip())
            panel[key] = row
    return panel


def match_variant(variant: Variant, panel: Dict[Tuple[str, str], Dict[str, str]],
                   min_vaf: float = DEFAULT_MIN_VAF) -> MatchResult:
    if not variant.gene or not variant.hgvs_p:
        return MatchResult(variant=variant, status="NOT_IN_PANEL")

    key = (variant.gene, variant.hgvs_p)
    panel_entry = panel.get(key)

    if panel_entry is None:
        return MatchResult(variant=variant, status="NOT_IN_PANEL")

    if variant.vaf is not None and variant.vaf < min_vaf:
        return MatchResult(variant=variant, status="LOW_VAF_MATCH", panel_entry=panel_entry)

    return MatchResult(variant=variant, status="MATCH", panel_entry=panel_entry)


def match_all_variants(variants: List[Variant], panel: Dict[Tuple[str, str], Dict[str, str]],
                        min_vaf: float = DEFAULT_MIN_VAF) -> List[MatchResult]:
    return [match_variant(v, panel, min_vaf=min_vaf) for v in variants]
