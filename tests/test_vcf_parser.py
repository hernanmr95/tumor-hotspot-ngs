import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.vcf_parser import parse_vcf, Variant

SAMPLE_VCF = """##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=HGVS_P,Number=1,Type=String,Description="Protein change">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR_SAMPLE
chr7\t140753336\t.\tA\tT\t60\tPASS\tGENE=BRAF;HGVS_P=p.Val600Glu\tGT:AD:DP\t0/1:45,55:100
chr1\t123456\t.\tC\tG\t60\tPASS\tGENE=FOO;HGVS_P=p.Xxx1Yyy\tGT:AD:DP\t0/1:90,10:100
"""


def _write_temp_vcf(content: str) -> Path:
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    f.write(content)
    f.close()
    return Path(f.name)


def test_parse_vcf_extracts_basic_fields():
    path = _write_temp_vcf(SAMPLE_VCF)
    variants = list(parse_vcf(path))

    assert len(variants) == 2
    v = variants[0]
    assert isinstance(v, Variant)
    assert v.chrom == "chr7"
    assert v.pos == 140753336
    assert v.ref == "A"
    assert v.alt == "T"
    assert v.gene == "BRAF"
    assert v.hgvs_p == "p.Val600Glu"


def test_parse_vcf_computes_depth_and_vaf():
    path = _write_temp_vcf(SAMPLE_VCF)
    variants = list(parse_vcf(path))

    v = variants[0]
    assert v.ref_depth == 45
    assert v.alt_depth == 55
    assert v.depth == 100
    assert v.vaf == 0.55


def test_parse_vcf_handles_low_vaf_variant():
    path = _write_temp_vcf(SAMPLE_VCF)
    variants = list(parse_vcf(path))

    v = variants[1]
    assert v.vaf == 0.1


def test_parse_vcf_skips_header_lines():
    content = "##fileformat=VCFv4.2\n##some=other\n#CHROM" + SAMPLE_VCF.split("#CHROM")[1]
    path = _write_temp_vcf(content)
    variants = list(parse_vcf(path))
    assert len(variants) == 2
