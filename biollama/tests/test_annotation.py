import sys
from biollama.core.annotation import CosmicLlama, LlamaEnsembl
import pandas as pd


def test_cosmic_id():
    cl = CosmicLlama()
    print(cl.query('COSM476'))


def test_annotation():
    llama = LlamaEnsembl()
    df = pd.DataFrame({'CHROM': ['chr1', 'chr17', '7', '7'], 'START': [153330766, 41244573, 116412043, 116412044], 'END': [153330797, 41245953, 116412043, 116412044],
                       'gene_exp': ['S100A9', 'BRCA1', 'MET', 'MET'], 'exon_exp': ['2', '10', '14', ''],
                       'strand_exp': ['+', '-', '+', '+']})
    ndf = llama.annotate_dataframe(df)
    mdf = pd.concat([df, ndf], axis=1)
    print(mdf)
    assert(mdf[mdf['genes'] != mdf['gene_exp']].shape[0] == 0)
    assert(mdf[mdf['exons'] != mdf['exon_exp']].shape[0] == 0)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_cosmic_id()
        test_annotation()
