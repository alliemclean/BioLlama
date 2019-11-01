"""
https://pypi.org/project/pyensembl/
https://grch37.rest.ensembl.org

"""

import numpy as np
import pandas as pd
from pyensembl import EnsemblRelease
import requests
import json
import sys


class LlamaEnsembl(object):
    """ Ensembl tools """
    def __init__(self, version=75):
        self.version = version
        self.db = EnsemblRelease(version)
        self.rest_url = "http://grch37.rest.ensembl.org"

    def rest_call(self, ext, data=None):
        if data:
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            r = requests.post(self.rest_url + ext, headers=headers, data=data)
        else:
            headers = {"Content-Type": "application/json"}
            r = requests.get(self.rest_url + ext, headers=headers)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        # print(repr(decoded))
        return decoded

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

    # Rest client calls
    def get_rsids(self, rsids):
        ext = "/variation/homo_sapiens"
        data = {"ids": rsids}
        return self.rest_call(ext, json.dumps(data))

    def get_cds_region(self, transcript, position):
        """ get location of variant to """
        ext = "/variation/human/{}:{}?".format(transcript, position)
        try:
            mappings = self.rest_call(ext)['mappings'][0]
        except requests.exceptions.HTTPError:
            return '', '', ''
        return mappings['seq_region_name'], mappings['start'], mappings['end']

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

    # Annotate DataFrames
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

    def annotate_variants(self, rsid_array, extra_cols=[]):
        """ Get chom:start-end for a list of variants """
        result = {'chrom': [], 'start': [], 'end': [], 'rsid': [], 'allele': [], 'vartype': [], 'consequence': []}
        for extra in extra_cols:
            result[extra] = []
        response = self.get_rsids(rsid_array)
        for var in rsid_array:
            mapping = response[var]['mappings'][0]
            result['chrom'].append(mapping['seq_region_name'])
            result['start'].append(mapping['start'])
            result['end'].append(mapping['end'])
            result['rsid'].append(var)
            result['allele'].append(mapping['allele_string'])
            result['vartype'].append(response[var]['var_class'])
            result['consequence'].append(response[var]['most_severe_consequence'])
            for extra in extra_cols:
                result[extra].append(response[var][extra])
        return pd.DataFrame(result)

    def annotate_cds_regions(self, df, tx_col='NM', cds_col='MutationName'):
        chroms = []
        starts = []
        ends = []
        for _, row in df.iterrows():
            location = self.get_cds_region(row[tx_col], row[cds_col])
            chroms.append(location[0])
            starts.append(location[1])
            ends.append(location[2])
        df['chrom'] = chroms
        df['start'] = starts
        df['end'] = ends
        return df


class CosmicResult(object):
    def __init__(self, decoded, cosmid):
        self.raw = decoded
        self.internal_id = decoded[0]
        record_names = decoded[1]
        records = decoded[3]
        self.records = {name: record for name, record in zip(record_names, records)}
        self.id = cosmid
        self.gene = self.records[cosmid][1]
        self.hgvs = self.records[cosmid][2]
        self.hgvs_p = self.records[cosmid][3]

    def __str__(self):
        return "id:{} gene:{} hgvs:{} hgvs_p:{}".format(self.id, self.gene, self.hgvs, self.hgvs_p)


class CosmicLlama(object):
    """ query clintable for cosmic ids """
    def __init__(self):
        self.url = "https://clinicaltables.nlm.nih.gov/api/cosmic/v3/search"

    def query(self, cosmid):
        qstring = "?terms={}".format(cosmid)
        res = requests.get("{}{}".format(self.url, qstring), headers={"Content-Type": "application/json"})
        if not res.ok:
            res.raise_for_status()
            sys.exit()

        decoded = res.json()
        # print(repr(decoded))
        return CosmicResult(decoded, cosmid)