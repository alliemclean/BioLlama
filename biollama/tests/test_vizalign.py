
from pkg_resources import resource_filename
from biollama.tools.vizalign import compare
import sys

def get_data():
    data = {}
    key = None
    seqfile = resource_filename('biollama.tests.data', 'sequences')
    with open(seqfile, 'r') as fileh:
        for line in fileh:
            if line.startswith('#'):
                key = line.strip().lstrip('#')
                data[key] = []
            else:
                if key is not None:
                    data[key].append(line.strip())
    return data

def get_seq(sarr):
    seqs = [line for line in sarr if not line.startswith('>')]
    return ''.join(seqs)

def test_vizalign():
    seqs_d = get_data()
    keys = list(seqs_d.keys())
    seqs = get_seq(seqs_d[keys[0]]), get_seq(seqs_d[keys[1]])
    compare(seqs[0], seqs[1])
    #compare('aaaaaatgcgcgcgc', 'aatgc')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_vizalign()
