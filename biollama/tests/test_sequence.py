import sys
import biollama.tools.sequence as seq
from biollama.core.annotation import LlamaEnsembl


def test_sequence():
    print(seq.Sequence())


def test_rsid():
    llama = LlamaEnsembl()
    rsids = """rs1517114
    rs4646
    rs55886062
    rs3918290
    rs67376798
    rs75017182
    rs115232898
    rs1801158
    rs11615
    rs1800566
    rs7779029
    rs151264360
    rs4148323
    rs8175347"""
    rsidstr = [r.lstrip() for r in rsids.split('\n')]
    df = llama.annotate_variants(rsidstr, extra_cols=['MAF', 'ambiguity'])
    print(df)


def test_cds_convert():
    llama = LlamaEnsembl()
    result = llama.get_cds_region('NM_015506.2', 'c.1A>G')
    result2 = llama.get_cds_region('NM_015506.2', 'c.445_446del')
    print(result)
    print(result2)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_cds_convert()
        test_sequence()
        test_rsid()
