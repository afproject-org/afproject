import argparse
import json
import sys
import itertools
import numpy as np

from libs import utils
from libs import treeutils

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', dest="input_file", required=True,
                    type=argparse.FileType('r'), metavar="FILE",
                    help="Input file")
parser.add_argument('-f', '--format', dest="format", choices=['tsv', 'phylip'],
                    help="File format", default='tsv')
args = parser.parse_args()


fh = open('datasets/genetree/swisstree/dataset.json')
DATASET = json.load(fh)
fh.close()


# Read and valite input file.
if args.format == 'tsv':
    err, dpair = utils.validate_tsv(args.input_file)
    if err:
        print(err)
        sys.exit(1)
else:
    err, matrix = utils.validate_phy(args.input_file)
    if err:
        print(err)
        sys.exit(1)
    dpair = utils.matrix_to_dpair(matrix)



REFERENCES = ['ST001', 'ST002', 'ST003', 'ST004', 'ST005', 'ST007', 'ST008',
'ST009', 'ST010', 'ST011', 'ST012']
nrfs = []
for r in REFERENCES:
    fh = open('datasets/genetree/swisstree/dataset/{}/ids.json'.format(r))
    IDS = json.load(fh)
    fh.close()
    fh = open('datasets/genetree/swisstree/dataset/{}/tree.ids.newick'.format(r))
    reftree = fh.read().strip()
    fh.close()


    ids = IDS["seqids"]
    l = []
    for id1, id2 in itertools.combinations(ids, 2):
        pair = tuple(sorted([id1, id2]))
        if pair not in dpair:
            err = "Input file does not have all sequence ids."
            print(err)
            sys.exit(1)
        l.append('{}\t{}\t{}'.format(pair[0], pair[1], dpair[pair]))
    d = utils.read_clean_tsv_data_dpair(l)
    matrix = utils.dpair_to_matrix(d, ids)
    phystr = matrix.format()
    usr_tree_str = treeutils.run_fneigbor(phystr)
    res = treeutils.ete_compare(
        usr_tree_str, reftree
    )
    qt, rt, nrf, rf, maxrf, src_br, ref_br, treekod = res
    print('{}\tnRF\t{}'.format(r, nrf))
    nrfs.append(nrf)

print('MEAN nRF: {}'.format(np.mean(nrfs)))
print('STDEV nRF: {}'.format(np.std(nrfs)))