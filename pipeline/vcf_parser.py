"""
pipeline/vcf_parser.py

Parser mínimo de VCF (v4.x) basado en la librería estándar, suficiente
para este pipeline: extrae CHROM/POS/REF/ALT, los campos de INFO
(GENE, HGVS_P) y la profundidad/VAF de la muestra tumoral desde FORMAT.

En un pipeline de producción se usaría `pysam` o `cyvcf2` (mucho más
robustos frente a la variabilidad real de VCFs de distintos callers),
pero para mantener el repo sin dependencias binarias pesadas, y dado
que el formato de entrada aquí es un VCF ya normalizado y anotado,
un parser simple y bien testeado es suficiente y más transparente.
"""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str]
    hgvs_p: Optional[str]
    ref_depth: Optional[int]
    alt_depth: Optional[int]

    @property
    def depth(self) -> Optional[int]:
        if self.ref_depth is None or self.alt_depth is None:
            return None
        return self.ref_depth + self.alt_depth

    @property
    def vaf(self) -> Optional[float]:
        """Variant Allele Fraction: fracción de reads que soportan el alelo alternativo."""
        if not self.depth:
            return None
        return round(self.alt_depth / self.depth, 4)


def _open_maybe_gzip(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_info(info_field: str) -> dict:
    info = {}
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def _parse_sample_depths(format_field: str, sample_field: str) -> tuple:
    """Extrae (ref_depth, alt_depth) del campo AD si está presente en FORMAT."""
    keys = format_field.split(":")
    values = sample_field.split(":")
    fields = dict(zip(keys, values))

    ad = fields.get("AD")
    if not ad or "," not in ad:
        return None, None

    parts = ad.split(",")
    try:
        ref_depth, alt_depth = int(parts[0]), int(parts[1])
        return ref_depth, alt_depth
    except (ValueError, IndexError):
        return None, None


def parse_vcf(path: Path) -> Iterator[Variant]:
    """Itera los objetos Variant de un VCF (opcionalmente gzip)."""
    with _open_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue  # línea de cabecera de columnas, ya conocemos el orden estándar

            fields = line.split("\t")
            if len(fields) < 8:
                continue  # línea malformada, se ignora de forma defensiva

            chrom, pos, _id, ref, alt, _qual, _filter, info_field = fields[:8]
            info = _parse_info(info_field)

            ref_depth = alt_depth = None
            if len(fields) >= 10:
                format_field, sample_field = fields[8], fields[9]
                ref_depth, alt_depth = _parse_sample_depths(format_field, sample_field)

            yield Variant(
                chrom=chrom,
                pos=int(pos),
                ref=ref,
                alt=alt,
                gene=info.get("GENE"),
                hgvs_p=info.get("HGVS_P"),
                ref_depth=ref_depth,
                alt_depth=alt_depth,
            )
