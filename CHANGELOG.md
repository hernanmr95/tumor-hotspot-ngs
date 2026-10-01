# Changelog

## [Unreleased]

### Planeado
- Soporte para VCFs multi-muestra (tumor/normal pareado).
- Normalización de variantes (left-alignment) antes del matching por coordenada.
- Carga del panel de datos desde una fuente externa.

## [1.0.0] - 2026-09-11

### Añadido
- Parser de VCF mínimo basado en stdlib (`pipeline/vcf_parser.py`).
- Matcher de variantes contra panel de hotspots por (gen, cambio de proteína).
- Clasificación explícita MATCH / LOW_VAF_MATCH / NOT_IN_PANEL.
- Generación de informe en Markdown con trazabilidad completa.
- Datos sintéticos de ejemplo (panel + VCF de muestra).
- Tests unitarios de parser y matcher.
