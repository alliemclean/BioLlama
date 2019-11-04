import pandas as pd

def get_pos_flds(region):
    """ return chrom, start, end from genomic region chr:start-end"""
    chrom, span = region.split(':')
    start, end = span.split('-')
    return chrom, start, end

def exon_df_from_ref(exon_starts, exon_ends, cds_start=0, cds_end=None, strand='+'):
    """
    :param exon_starts:
    :param exon_ends:
    :param cds_start:
    :param cds_end:
    :param strand:
    :return: dataframe
    if cds_start and end aren't supplied cds_status is unset.  if it is
      1 = entirely in coding region
      0 = not in coding region
      .5 = partially in coding region
    """
    starts = [int(exon) for exon in exon_starts.split(',') if exon != '']
    ends = [int(exon) for exon in exon_ends.split(',') if exon != '']
    exon_count = len(ends)
    cds_length = 0
    dd = {k: [] for k in ['exon_id', 'cds_pos', 'start', 'end', 'cds_status']}
    cds_unset = False
    if cds_end is None:
        cds_end = ends[-1] + 1
        cds_unset = True
    for i, (s, e) in enumerate(zip(starts, ends)):
        if s < cds_start:
            if e < cds_start:
                cds_status = 0
            else:
                cds_length += e - cds_start
                cds_status = .5
        elif e > cds_end:
            if s > cds_end:
                cds_status = 0
            else:
                cds_length += cds_end - s
                cds_status = .5
        else:
            cds_length += e - s
            cds_status = 1
        if strand == '+':
            exon = i + 1
        else:
            exon = exon_count - i
        dd['exon_id'].append(exon)
        dd['cds_pos'].append(cds_length)
        dd['start'].append(s)
        dd['end'].append(e)
        dd['cds_status'].append(cds_status)
    df = pd.DataFrame(dd)
    df['strand'] = strand
    df['cds_length'] = cds_length
    if cds_unset:
        df['cds_status'] = 'No CDS info'
    return df