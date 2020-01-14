

__author__ = "Allison MacLeay"

import argparse
import pandas as pd
from biollama.core.annotation import UCSCapi, LlamaEnsembl


def annotate(df, cols):
    """ get longest transcript for gene column or regions """
    def get_regions(row, cols):
        prefix = '' if row[cols[0]].startswith('chr') else 'chr'
        return "{}{}:{}-{}".format(prefix, row[cols[0]], row[cols[1]], row[cols[2]])

    llama = LlamaEnsembl()
    ucsc = UCSCapi()
    dd = {'query': [], 'transcripts': []}
    if len(cols) == 1:
        genes = df[cols[0]].values
        for gene in genes:
            chrom, start, end = llama.get_gene_pos(gene)
            res = ucsc.query("chr{}:{}-{}".format(chrom, start, end))
            dd['query'].append(gene)
            dd['transcripts'].append(res.longest(gene)['transcript'])
            mergecol = cols[0]
    elif len(cols) == 3:
        dd['genes'] = []
        df['region'] = df.apply(get_regions, axis=1, args=[cols])
        regions = df['region'].values
        for region in regions:
            res = ucsc.query(region)
            for gene in res.genes():
                dd['genes'].append(gene)
                dd['query'].append(region)
                dd['transcripts'].append(res.longest(gene)['transcript'])
                mergecol = 'region'
    ndf = df.merge(pd.DataFrame(dd), left_on=mergecol, right_on='query', how='outer')
    return ndf[~ndf.duplicated()]


def main(input, output, columns):
    df = pd.read_csv(input, sep='\t')
    cols = columns.split(',')
    ndf = annotate(df, cols)
    ndf.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input table file')
    parser.add_argument('--output', help='output file')
    parser.add_argument('--columns', default='gene', help='column to annotate.  Default is "gene".  If there are 3 '
                                                          'columns it assumes they are chrom,start,end in that order '
                                                          'and will treat the column names as such.')
    args = parser.parse_args()
    main(args.input, args.output, args.columns)