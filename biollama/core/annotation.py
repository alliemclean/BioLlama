import numpy as np
import pandas as pd
from pyensembl import EnsemblRelease


class LlamaEnsembl(object):
    """ Ensembl tools """
    def __init__(self, version=75):
        self.version = version
        self.db = EnsemblRelease(version)

    def load_ensembl_ref(self, rid=None):
        """ Download, load, and index ensembl data """
        self.db.download(self.version)
        self.db.index()
        if rid is not None:
            return self.db.transcript_by_id(rid)
        else:
            return None

    def get_exon_numbers(self, gene):
        """ This creates exon areas from the biggest transcript """
        dct = {'start': [], 'id': [], 'stop': [], 'transcript': []}
        gene_id = self.db.gene_ids_of_gene_name(gene)[0]
        transcripts = self.db.transcript_ids_of_gene_id(gene_id)
        longest = 0
        e = None
        for trans in transcripts:
            tsc = self.db.exon_ids_of_transcript_id(trans)
            tsize = len(tsc)
            if tsize > longest:
                longest = tsize
                e = tsc
                longest_transcript = trans
        for exid in e:
            exon = self.db.exon_by_id(exid)
            dct['start'].append(exon.start)
            dct['stop'].append(exon.end)
            dct['id'].append(exid)
            dct['transcript'].append(longest_transcript)
        df = pd.DataFrame(dct)
        df['number'] = df.index + 1
        return df

    def get_genes(self, chrom, start, stop):
        return [gobj.gene_name for gobj in self.db.genes_at_locus(chrom, start, stop)]

    def parse_ref_exons(self, chrom, start, stop):
        """ Return fasta reference with only the sequences needed"""
        ens_db = self.db
        try:
            exons = ens_db.exons_at_locus(chrom, start, stop)
        except ValueError as e:
            # Load pyensembl db
            raise e
        if not len(exons):
            return '', ''
        exon_numbers = self.get_exon_numbers(exons[0].gene_name)
        transcript = exon_numbers['transcript'].values[0]
        trx_exons = []
        for ex in exons:
            nrow = exon_numbers[exon_numbers['id'] == ex.exon_id]
            if nrow.shape[0] > 0:
                trx_exons.extend(nrow['number'].values)
        return transcript, ','.join([str(number) for number in trx_exons])

    def annotate_dataframe(self, df, chrom_col='CHROM', start_col='START', end_col='END'):
        genes = []
        exons = []
        transcripts = []
        for i, row in df.iterrows():
            genes_row = self.get_genes(row[chrom_col], row[start_col], row[end_col])
            genes.append(','.join(genes_row))
            if len(genes_row) < 2:
                trans_row, exons_row = self.parse_ref_exons(row[chrom_col], row[start_col], row[end_col])
            else:
                trans_row = ''
                exons_row = ''
            exons.append(exons_row)
            transcripts.append(trans_row)
        new_df = pd.DataFrame({'genes': genes, 'exons': exons, 'transcript': transcripts}, index=df.index)
        return new_df