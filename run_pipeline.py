#!/usr/bin/env python3
"""
run_pipeline.py

Punto de entrada del pipeline de anotación de marcadores tumorales.

Flujo:
  1. Carga el panel curado de hotspots (CSV).
  2. Parsea el VCF de la muestra (ya anotado con GENE/HGVS_P en INFO).
  3. Cruza cada variante contra el panel (por gen + cambio de proteína).
  4. Genera un informe en Markdown con hallazgos, casos de VAF baja para
     revisión manual, y trazabilidad completa de variantes fuera de panel.

Uso:
    python run_pipeline.py \
        --vcf data/sample_tumor.vcf \
        --panel data/hotspot_panel.csv \
        --sample-name "TUMOR_SAMPLE" \
        --min-vaf 0.05 \
        --outdir results/
"""

import argparse
from pathlib import Path

from pipeline.vcf_parser import parse_vcf
from pipeline.hotspot_matcher import load_hotspot_panel, match_all_variants, DEFAULT_MIN_VAF
from pipeline.report import build_report, write_report


def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline de anotación de marcadores tumorales (demo)")
    parser.add_argument("--vcf", default="data/sample_tumor.vcf",
                         help="VCF de la muestra, ya anotado con GENE/HGVS_P")
    parser.add_argument("--panel", default="data/hotspot_panel.csv",
                         help="CSV con el panel curado de hotspots")
    parser.add_argument("--sample-name", default="TUMOR_SAMPLE",
                         help="Nombre de la muestra para el encabezado del informe")
    parser.add_argument("--min-vaf", type=float, default=DEFAULT_MIN_VAF,
                         help="VAF mínima para considerar un match como MATCH (no LOW_VAF_MATCH)")
    parser.add_argument("--outdir", default="results",
                         help="Directorio de salida para el informe")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(">> Cargando panel de hotspots...")
    panel = load_hotspot_panel(Path(args.panel))
    print(f"   {len(panel)} entradas en el panel")

    print(">> Parseando VCF de la muestra...")
    variants = list(parse_vcf(Path(args.vcf)))
    print(f"   {len(variants)} variantes leídas")

    print(">> Cruzando variantes contra el panel...")
    results = match_all_variants(variants, panel, min_vaf=args.min_vaf)

    n_match = sum(1 for r in results if r.status == "MATCH")
    n_low_vaf = sum(1 for r in results if r.status == "LOW_VAF_MATCH")
    n_background = sum(1 for r in results if r.status == "NOT_IN_PANEL")
    print(f"   {n_match} match(es), {n_low_vaf} con VAF baja, {n_background} fuera de panel")

    print(">> Generando informe...")
    report_text = build_report(results, sample_name=args.sample_name, min_vaf=args.min_vaf)
    out_path = outdir / f"{args.sample_name}_report.md"
    write_report(report_text, out_path)

    print(f"\nListo. Informe en: {out_path}")


if __name__ == "__main__":
    main()
