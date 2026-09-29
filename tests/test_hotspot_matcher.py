import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.vcf_parser import Variant
from pipeline.hotspot_matcher import load_hotspot_panel, match_variant, match_all_variants

SAMPLE_PANEL_CSV = """gene,protein_change,variant_type,associated_cancers,clinical_significance,notes
BRAF,p.Val600Glu,missense,Melanoma,Pathogenic,hotspot clasico
KRAS,p.Gly12Asp,missense,"Cancer colorrectal",Pathogenic,hotspot frecuente
"""


def _write_temp_panel() -> Path:
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="")
    f.write(SAMPLE_PANEL_CSV)
    f.close()
    return Path(f.name)


def _variant(gene, hgvs_p, ref_depth=50, alt_depth=50):
    return Variant(
        chrom="chr1", pos=1, ref="A", alt="T",
        gene=gene, hgvs_p=hgvs_p,
        ref_depth=ref_depth, alt_depth=alt_depth,
    )


def test_load_hotspot_panel_indexes_by_gene_and_protein_change():
    panel = load_hotspot_panel(_write_temp_panel())
    assert ("BRAF", "p.Val600Glu") in panel
    assert ("KRAS", "p.Gly12Asp") in panel
    assert len(panel) == 2


def test_match_variant_returns_match_for_known_hotspot_with_sufficient_vaf():
    panel = load_hotspot_panel(_write_temp_panel())
    v = _variant("BRAF", "p.Val600Glu", ref_depth=50, alt_depth=50)  # VAF 0.5

    result = match_variant(v, panel, min_vaf=0.05)

    assert result.status == "MATCH"
    assert result.panel_entry["associated_cancers"] == "Melanoma"


def test_match_variant_returns_low_vaf_match_below_threshold():
    panel = load_hotspot_panel(_write_temp_panel())
    v = _variant("BRAF", "p.Val600Glu", ref_depth=97, alt_depth=3)  # VAF 0.03

    result = match_variant(v, panel, min_vaf=0.05)

    assert result.status == "LOW_VAF_MATCH"


def test_match_variant_returns_not_in_panel_for_unknown_variant():
    panel = load_hotspot_panel(_write_temp_panel())
    v = _variant("TP53", "p.Arg175His")  # no está en este panel de test

    result = match_variant(v, panel, min_vaf=0.05)

    assert result.status == "NOT_IN_PANEL"
    assert result.panel_entry is None


def test_match_variant_does_not_match_on_gene_alone():
    """Un cambio de proteína distinto en el mismo gen NO debe dar match."""
    panel = load_hotspot_panel(_write_temp_panel())
    v = _variant("BRAF", "p.Asp594Asn")  # mismo gen, otro cambio

    result = match_variant(v, panel, min_vaf=0.05)

    assert result.status == "NOT_IN_PANEL"


def test_match_all_variants_processes_a_list():
    panel = load_hotspot_panel(_write_temp_panel())
    variants = [
        _variant("BRAF", "p.Val600Glu"),
        _variant("UNKNOWN", "p.Foo1Bar"),
    ]
    results = match_all_variants(variants, panel, min_vaf=0.05)

    assert len(results) == 2
    assert results[0].status == "MATCH"
    assert results[1].status == "NOT_IN_PANEL"
