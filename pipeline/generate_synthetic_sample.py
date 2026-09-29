"""
pipeline/generate_synthetic_sample.py

Genera un VCF sintético (data/sample_tumor.vcf) que simula la salida de
un pipeline de llamada + anotación de variantes (p. ej. GATK Mutect2 +
VEP/SnpEff) sobre una muestra tumoral. El pipeline de este repo empieza
DESPUÉS de ese paso: asume un VCF ya anotado con GENE y HGVS_P en el
campo INFO, tal como lo dejaría VEP/SnpEff en un flujo real.

Incluye:
  - Algunas variantes que coinciden exactamente con el panel de hotspots
    (data/hotspot_panel.csv), con VAFs realistas de tumor.
  - Variantes "de fondo" en los mismos genes pero que NO son hotspots
    conocidos (para comprobar que el matcher no genera falsos positivos
    por coincidir solo el gen).
  - Una variante con VAF muy baja (posible artefacto/subclonal) para
    ejercitar el filtrado por profundidad/VAF mínima.

Las posiciones genómicas son inventadas (no corresponden a coordenadas
reales de referencia) porque el matching de este pipeline se hace por
gen + cambio de proteína (HGVS_P), no por coordenada — ver README para
la justificación de este diseño.
"""

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

HEADER = """##fileformat=VCFv4.2
##source=tumor-hotspot-ngs-synthetic-generator
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=HGVS_P,Number=1,Type=String,Description="Protein change (HGVS p. notation)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR_SAMPLE
"""

# Cada tupla: (chrom, pos, ref, alt, gene, hgvs_p, ref_depth, alt_depth)
VARIANTS = [
    # --- Hotspots conocidos, VAF razonable de tumor ---
    ("chr7",  140753336, "A", "T", "BRAF",   "p.Val600Glu",   45, 55),
    ("chr7",   55191822, "T", "G", "EGFR",   "p.Leu858Arg",   60, 40),
    ("chr12",  25245350, "C", "T", "KRAS",   "p.Gly12Asp",    50, 50),
    ("chr3",  178952085, "A", "G", "PIK3CA", "p.His1047Arg",  70, 30),
    ("chr17",   7577120, "C", "T", "TP53",   "p.Arg175His",   65, 35),

    # --- Variante hotspot con VAF muy baja (posible subclonal/artefacto) ---
    ("chr1",  114716126, "C", "T", "NRAS",   "p.Gln61Arg",    97, 3),

    # --- Variantes "de fondo": mismo gen, pero NO son el hotspot catalogado ---
    ("chr7",  140753301, "G", "A", "BRAF",   "p.Asp594Asn",   80, 20),
    ("chr12",  25245347, "G", "A", "KRAS",   "p.Ala11Thr",    75, 25),
    ("chr17",   7577558, "G", "A", "TP53",   "p.Pro72Leu",    90, 10),

    # --- Variantes en genes fuera del panel, ruido normal de la muestra ---
    ("chr2",   29223528, "C", "G", "ALK",    "p.Phe1174Leu",  85, 15),
    ("chr9",   21971120, "C", "T", "CDKN2A", "p.Arg80*",      88, 12),
]


def make_vcf_line(chrom, pos, ref, alt, gene, hgvs_p, ref_depth, alt_depth):
    dp = ref_depth + alt_depth
    info = f"GENE={gene};HGVS_P={hgvs_p}"
    genotype = "0/1" if alt_depth < dp else "1/1"
    sample_field = f"{genotype}:{ref_depth},{alt_depth}:{dp}"
    return "\t".join([
        chrom, str(pos), ".", ref, alt, "60", "PASS",
        info, "GT:AD:DP", sample_field,
    ])


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_path = DATA_DIR / "sample_tumor.vcf"

    with open(out_path, "w") as f:
        f.write(HEADER)
        for variant in VARIANTS:
            f.write(make_vcf_line(*variant) + "\n")

    print(f"VCF sintético escrito en: {out_path} ({len(VARIANTS)} variantes)")


if __name__ == "__main__":
    main()
