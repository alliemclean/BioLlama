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


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_sequence()
        test_rsid()
