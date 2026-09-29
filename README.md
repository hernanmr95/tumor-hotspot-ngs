# tumor-hotspot-ngs

Pipeline en Python puro que cruza variantes de una muestra tumoral
(VCF ya anotado) contra un panel curado de mutaciones "hotspot" bien
establecidas en oncología molecular, y genera un informe trazable con
Markdown.

## ⚠️ Antes de nada: qué es y qué no es esto

Ver [`docs/scope_and_limitations.md`](docs/scope_and_limitations.md) —
en resumen: esto es una demostración técnica, no un dispositivo de
diagnóstico clínico. Los hotspots incluidos (BRAF V600E, EGFR L858R,
KRAS G12D, etc.) son marcadores reales y bien documentados en la
literatura, pero para uso clínico real hace falta validación analítica
formal y certificación regulatoria.

## Qué hace:
1. Carga un panel curado de hotspots (`data/hotspot_panel.csv`): gen,
   cambio de proteína, cánceres asociados, significancia clínica.
2. Parsea un VCF de la muestra que ya viene anotado (campo INFO con
   `GENE` y `HGVS_P`) — este pipeline asume que el llamado de variantes y
   la anotación funcional (tipo GATK Mutect2 + VEP/SnpEff) ya se hicieron
   antes; no los reimplementa.
3. Cruza cada variante contra el panel por (gen, cambio de proteína).
4. Clasifica cada variante en `MATCH`, `LOW_VAF_MATCH` (coincide pero con
   VAF por debajo del umbral — para revisión manual) o `NOT_IN_PANEL`.
5. Genera un informe en Markdown con trazabilidad completa: no solo los
   hallazgos positivos, también las variantes fuera de panel.

## Estructura:
.
├── run_pipeline.py                 # Orquestador
├── pipeline/
│   ├── vcf_parser.py               # Parser de VCF (stdlib, sin dependencias pesadas)
│   ├── hotspot_matcher.py          # Lógica de cruce variante ↔ panel
│   ├── report.py                   # Generación del informe Markdown
│   └── generate_synthetic_sample.py  # Genera el VCF de ejemplo
├── data/
│   ├── hotspot_panel.csv           # Panel curado (conocimiento público establecido)
│   └── sample_tumor.vcf            # VCF sintético de ejemplo
├── docs/
│   └── scope_and_limitations.md    # Qué es y qué NO es esta herramienta
├── tests/
│   ├── test_vcf_parser.py
│   └── test_hotspot_matcher.py
└── .github/workflows/ci.yml

## Uso:
```bash
pip install -r requirements.txt   # solo pytest; el pipeline usa stdlib

# (opcional) regenerar la muestra sintética
python pipeline/generate_synthetic_sample.py

# correr el pipeline
python run_pipeline.py \
  --vcf data/sample_tumor.vcf \
  --panel data/hotspot_panel.csv \
  --sample-name "DEMO_001" \
  --min-vaf 0.05 \
  --outdir results/

cat results/DEMO_001_report.md
```

## Datos de ejemplo:
`data/sample_tumor.vcf` es sintético: incluye 5 hotspots conocidos con
VAF de tumor realista, 1 hotspot con VAF muy baja (para ejercitar el
caso `LOW_VAF_MATCH`), variantes de "fondo" en los mismos genes que
deliberadamente NO son el hotspot catalogado (para verificar que no hay
falsos positivos por coincidir solo el gen), y variantes en genes fuera
del panel. Las coordenadas genómicas son inventadas.
