"""
pipeline/report.py

Genera un informe en Markdown a partir de los MatchResult del pipeline.
El formato imita (deliberadamente, a nivel de estructura) el tipo de
informe que produciría un pipeline de tumor profiling: resumen ejecutivo,
tabla de hallazgos, y una sección de variantes fuera del panel para
trazabilidad completa (nada se descarta en silencio).
"""

from datetime import datetime, timezone
from pathlib import Path
from typing import List

from pipeline.hotspot_matcher import MatchResult

DISCLAIMER = (
    "> **Aviso**: este informe es generado por una herramienta de demostración "
    "con fines educativos y de portfolio. NO ha sido validada analíticamente "
    "ni cuenta con certificación regulatoria (CE-IVD / FDA / CLIA), y **no debe "
    "usarse para ninguna decisión clínica real**. Ver `docs/scope_and_limitations.md`."
)


def _format_match_row(result: MatchResult) -> str:
    v, p = result.variant, result.panel_entry
    vaf_pct = f"{v.vaf * 100:.1f}%" if v.vaf is not None else "N/D"
    return (
        f"| {v.gene} | {v.hgvs_p} | {vaf_pct} | {v.depth or 'N/D'} "
        f"| {p['associated_cancers']} | {p['clinical_significance']} |"
    )


def _format_background_row(result: MatchResult) -> str:
    v = result.variant
    vaf_pct = f"{v.vaf * 100:.1f}%" if v.vaf is not None else "N/D"
    return f"| {v.chrom}:{v.pos} | {v.gene or 'N/D'} | {v.hgvs_p or 'N/D'} | {vaf_pct} |"


def build_report(results: List[MatchResult], sample_name: str, min_vaf: float) -> str:
    matches = [r for r in results if r.status == "MATCH"]
    low_vaf = [r for r in results if r.status == "LOW_VAF_MATCH"]
    background = [r for r in results if r.status == "NOT_IN_PANEL"]

    lines = []
    lines.append(f"# Informe de marcadores tumorales — {sample_name}")
    lines.append("")
    lines.append(f"Generado: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}  ")
    lines.append(f"Variantes analizadas: {len(results)}  ")
    lines.append(f"Umbral mínimo de VAF: {min_vaf * 100:.0f}%")
    lines.append("")
    lines.append(DISCLAIMER)
    lines.append("")

    lines.append("## Resumen")
    lines.append("")
    lines.append(f"- **{len(matches)}** hotspot(s) del panel detectado(s) con VAF suficiente.")
    lines.append(f"- **{len(low_vaf)}** coincidencia(s) con el panel pero con VAF por debajo del umbral "
                  f"(revisión manual recomendada).")
    lines.append(f"- **{len(background)}** variante(s) fuera del panel (posible ruido, variantes benignas, "
                  f"o hotspots reales no incluidos en este panel de demostración).")
    lines.append("")

    if matches:
        lines.append("## Hallazgos: hotspots detectados")
        lines.append("")
        lines.append("| Gen | Cambio proteico | VAF | Profundidad | Cánceres asociados | Significancia |")
        lines.append("|---|---|---|---|---|---|")
        for r in matches:
            lines.append(_format_match_row(r))
        lines.append("")

    if low_vaf:
        lines.append("## Coincidencias con VAF baja (revisar manualmente)")
        lines.append("")
        lines.append("| Gen | Cambio proteico | VAF | Profundidad | Cánceres asociados | Significancia |")
        lines.append("|---|---|---|---|---|---|")
        for r in low_vaf:
            lines.append(_format_match_row(r))
        lines.append("")

    if background:
        lines.append("## Variantes fuera del panel (trazabilidad completa)")
        lines.append("")
        lines.append("| Posición | Gen | Cambio proteico | VAF |")
        lines.append("|---|---|---|---|")
        for r in background:
            lines.append(_format_background_row(r))
        lines.append("")

    return "\n".join(lines)


def write_report(report_text: str, out_path: Path):
    out_path.write_text(report_text, encoding="utf-8")
