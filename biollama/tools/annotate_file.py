"""
http://genome.ucsc.edu/cgi-bin/hgBlat

"""
import pandas as pd

from argparse import ArgumentParser

from biollama.core.annotation import LlamaEnsembl


HG_BLAT_COLS = ['ACTIONS', 'QUERY', 'SCORE', 'START', 'END', 'QSIZE', 'IDENTITY', 'CHROM', 'STRAND',
                              'START', 'END', 'SPAN']


class SeqDoc(object):
    """ a sequence document to be annotated """
    def __init__(self, df=None):
        self.df = df

    def read_file(self):
        """ Not defined """
        raise NotImplementedError

    def write(self, fpath):
        self.df.to_csv(fpath, sep='\t', index=False)

    def annotate(self):
        """
        Annotate data frame
        :param df: pandas dataframe
        :return: merged dataframe
        """
        llama = LlamaEnsembl()
        new_df = llama.annotate_dataframe(self.df)
        return self.df.join(new_df)

    @staticmethod
    def detect_type(fpath):
        with open(fpath, 'r') as fhead:
            lines = [fhead.readline()]
        if all([word in lines[0] for word in HG_BLAT_COLS]):
            return HgBlatDoc
        cols = lines[0].strip().split('\t')
        try:
            if isinstance(int(cols[1]), int) and isinstance(int(cols[2]), int) and cols[5] in ['-', '+']:
                return BedDoc
        except ValueError:
            pass
        return None


class HgBlatDoc(SeqDoc):
    """ copied and pasted data from Blat output at UCSC
        http://genome.ucsc.edu/cgi-bin/hgBlat
    """
    def __init__(self, pfile):
        super(HgBlatDoc, self).__init__()
        self.df = self.read_file(pfile)

    def read_file(self, pfile):
        df = pd.read_csv(pfile, sep="\s+", names=self.columns, skiprows=2)
        df.index = range(df.shape[0])
        return df

    @property
    def columns(self):
        """ change first start and end columns to match_start, match_end"""
        for i, col in enumerate(HG_BLAT_COLS):
            if col == 'START':
                ind_start = i
            if col == 'END':
                ind_end = i
                break
        names = HG_BLAT_COLS[:]
        names[ind_start] = 'MATCH_START'
        names[ind_end] = 'MATCH_END'
        return names


class BedDoc(SeqDoc):
    """ bed file
        6 columns with no header
    """
    def __init__(self, pfile):
        super(BedDoc, self).__init__()
        self.df = self.read_file(pfile)

    def read_file(self, pfile):
        df = pd.read_csv(pfile, sep="\t", names=self.columns, header=None)
        return df

    @property
    def columns(self):
        """ change first start and end columns to match_start, match_end"""
        return ['CHROM', 'START', 'END', 'NAME', 'INFO', 'STRAND']


def main(args=None):
    """The main method"""
    parser = ArgumentParser(description="Generates gene and exon locations from copied blat info text at "
                                        "http://genome.ucsc.edu/cgi-bin/hgBlat")
    parser.add_argument('--file', help="input file")
    parser.add_argument('--out', help="output file")
    args = parser.parse_args()
    doctype = SeqDoc.detect_type(args.file)
    doc = doctype(args.file)
    new_doc = SeqDoc(doc.annotate())
    new_doc.write(args.out)



if __name__ == "__main__":
    main()