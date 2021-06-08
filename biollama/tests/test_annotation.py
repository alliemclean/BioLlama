import sys
from biollama.core.annotation import CosmicLlama, LlamaEnsembl, UCSCapi
from biollama.tools.get_gene_transcripts import annotate
import pandas as pd


def test_cosmic_id():
    cl = CosmicLlama()
    print(cl.query('COSV51769364'))


def test_cosmic_id_v3():
    cl = CosmicLlama(version="v3")
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


def test_ucsc_annotation():
    # 0-based as opposed to ensembl which is 1-based query
    llama = UCSCapi()
    df = pd.DataFrame({'chrom': ['chr1', 'chr1', 'chr1', 'chr17', '7', '7'], 'start': [2139276, 3585000, 153330766, 41244573, 116412042, 116412044], 'end': [2139437, 3585163, 153330797, 41245953, 116412043, 116412044],
                       'gene_exp': ['FAAP20', 'TP73', 'S100A9', 'BRCA1', 'MET', 'MET'], 'exon_exp': ['I2', 'I1', '2', '10', '14', 'I14'],
                       'strand_exp': ['-', '+', '+', '-', '+', '+']})
    ndf = llama.annotate_dataframe(df)
    mdf = pd.concat([df, ndf], axis=1)
    print(mdf)
    assert(mdf[mdf['gene'] != mdf['gene_exp']].shape[0] == 0)
    assert(mdf[mdf['exons'] != mdf['exon_exp']].shape[0] == 0)


def test_annotation_hg38():
    llama = LlamaEnsembl(genome='hg38')
    df = pd.DataFrame({'CHROM': ['chr1', 'chr17', '7', '7'], 'START': [153358330, 43092618, 116771890, 116773071], 'END': [153358350, 43092648, 116771920, 116773072],
                       'gene_exp': ['S100A9', 'BRCA1', 'MET', 'MET'], 'exon_exp': ['2', '10', '14', ''],
                       'strand_exp': ['+', '-', '+', '+']})
    ndf = llama.annotate_dataframe(df)
    mdf = pd.concat([df, ndf], axis=1)
    print(mdf)
    assert(mdf[mdf['genes'] != mdf['gene_exp']].shape[0] == 0)
    assert(mdf[mdf['exons'] != mdf['exon_exp']].shape[0] == 0)


def test_ucsc_annotation_hg38():
    # 0-based as opposed to ensembl which is 1-based query
    llama = UCSCapi(genome='hg38')
    df = pd.DataFrame({'chrom': ['chr1', 'chr1', 'chr1', 'chr17', '7', '7'], 'start': [2204473, 3665788, 153358367, 43092374, 116771988, 116771990], 'end': [2204493, 3665808, 153358397, 43092394, 116771989, 116771991],
                       'gene_exp': ['FAAP20', 'TP73', 'S100A9', 'BRCA1', 'MET', 'MET'], 'exon_exp': ['I2', 'I1', '2', '10', '14', 'I14'],
                       'strand_exp': ['-', '+', '+', '-', '+', '+']})
    ndf = llama.annotate_dataframe(df)
    mdf = pd.concat([df, ndf], axis=1)
    print(mdf)
    assert(mdf[mdf['gene'] != mdf['gene_exp']].shape[0] == 0)
    assert(mdf[mdf['exons'] != mdf['exon_exp']].shape[0] == 0)


def test_ucsc():
    ucsc = UCSCapi()
    res = ucsc.query("chr20:39788239-39788373")
    transcript = res.longest()
    try:
        assert(transcript['transcript'] == 'NM_002660.3')
    except AssertionError as e:
        print("{} is not equal to {}".format(transcript['transcript'], "NM_002660.3"))
        raise e


def test_gene_transcripts():
    llama = LlamaEnsembl()
    ucsc = UCSCapi()
    gene = 'BRCA1'
    chrom, start, end = llama.get_gene_pos(gene)
    res = ucsc.query("chr{}:{}-{}".format(chrom, start, end))
    transcript = res.longest(gene)
    print(transcript['transcript'])


def test_table_annotation():
    df = pd.DataFrame({'chrom': ['17', 'chr20'], 'start': [41196312, 39765265],
                       'end': [41277500, 39767773], 'gene': ['BRCA1', 'PLCG1']})
    ndf1 = annotate(df, ['gene'])
    ndf2 = annotate(df, ['chrom', 'start', 'end'])
    assert(ndf1.shape[0] == 2)
    assert(ndf2.shape[0] == 3)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_cosmic_id()
        test_cosmic_id_v3()
        test_ucsc_annotation()
        test_annotation()
        test_ucsc_annotation_hg38()
        test_annotation_hg38()
        test_ucsc()
        test_gene_transcripts()
        test_table_annotation()
