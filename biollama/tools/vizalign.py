from argparse import ArgumentParser
from itertools import product
import numpy as np

from plotly.offline import plot
import plotly.graph_objs as go

def size_sort(seq1, seq2):
    """ return longer sequence before shorter one """
    if len(seq2) > len(seq1):
        tmp = seq1
        seq1 = seq2
        seq2 = tmp
    return seq1.upper(), seq2.upper()


def get_matrix(seqa, seqb):
    cmp = np.zeros((len(seqa), len(seqb)))
    for i, j in product(range(len(seqa)), range(len(seqb))):
        a = seqa[i]
        b = seqb[j]
        cmp[i][j] = 1 if a == b else 0
    return cmp

def heatmap(marr, seqa, seqb):
    data = [go.Heatmap(z=marr.T, x=[l for l in seqa], y=[m for m in seqb])]
    plot({"data": data,
          "layout": go.Layout(title="Sequence Similarity", xaxis={"scaleratio": 1}, yaxis={"scaleanchor": "x"})})
    pass

def compare(seq1, seq2):
    seqa, seqb = size_sort(seq1, seq2)
    cmpm = get_matrix(seqa, seqb)
    heatmap(cmpm, seqa, seqb)


def main(args):
    parser = ArgumentParser()
    parser.add_argument('seq1')
    parser.add_argument('seq2')
    args = parser.parse_args()
    compare(args.seq1, args.seq2)


if __name__ == "__main__":
    main()