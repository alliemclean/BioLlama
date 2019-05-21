import os
import subprocess as sp
import time

from plotly.offline import plot
import plotly.graph_objs as go

import biollama.tools.sequence as seq

def search(sra_path, sequence, score=100, max=5, sample='SRR390728'):
    """ run search algorithm and benchmark time """
    start = time.time()
    cmd = "{}/sra-search \"{}\" -S {} -m {} -a SmithWaterman {}".format(sra_path, sequence, score, max, sample)
    print(cmd)
    print(sp.check_output(cmd, shell=True))
    end = time.time()
    return end - start


def time_search(sra_path, limit):
    times_length = time_search_length(sra_path, limit)
    times_max = time_search_max(sra_path, limit)
    print(times_length)
    print(times_max)
    data = [go.Scatter(x=list(times_length.keys()), y=list(times_length.values()))]
    plot({"data": data})

def time_search_max(sra_path, limit):
    times = {}
    sequence = seq.Sequence(length=5)
    for i in range(5, 5 + limit):
        times[i] = search(sra_path, sequence, max=i)
    return times


def time_search_length(sra_path, limit):
    times = {}
    for i in range(limit):
        length = i * 2 + 1
        sequence = seq.Sequence(length)
        times[length] = search(sra_path, sequence, max=10)
    return times


if __name__ == '__main__':
    """These are just to test the SRA toolkit
    https://trace.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc&f=std
    """
    sra_path = os.getenv('SRA')  # Path to SRA toolkit bin directory
    time_search(sra_path, 5)


