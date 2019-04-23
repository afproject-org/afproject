import argparse
import json
import sys

from libs import crm
from libs import utils, treeutils

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', dest="input_file", required=True,
                    type=argparse.FileType('r'), metavar="FILE",
                    help="Input file")
parser.add_argument('-f', '--format', dest="format", choices=['tsv', 'phylip', 'newick'],
                    help="File format", default='tsv')
parser.add_argument('-r', '--reference', dest="reference", 
                    choices=['ecoli', 'plants', 'ecoli_shigella', 'yersinia', 'sim_hgt', 'fish_mito'],
                    help="File format", default='tsv')
args = parser.parse_args()


# Read reference
fh = open('datasets/genome/{}/dataset.json'.format(args.reference))
j = json.load(fh)
fh.close()


# Read and valite input file.
if args.format == 'newick':
    err, tree = utils.validate_newick(args.input_file)
    if err:
        print(err)
        sys.exit(1)
    err, usr_tree_str = utils.validate_tree(tree, j['seqids'])
    if err:
        print(err)
        sys.exit(1)
else:
    if args.format == 'tsv':
        err, dpair = utils.validate_tsv(args.input_file)
        if err:
            print(err)
            sys.exit(1)
        err, clean_data = utils.validate_dpair(dpair, j['seqids'])
        if err:
            print(err)
            sys.exit(1)
        matrix = utils.dpair_to_matrix(dpair, j['seqids'])       
        phystr = matrix.format()
    else:
        err, matrix = utils.validate_phy(args.input_file)
        if err:
            print(err)
            sys.exit(1)
        err, phystr = utils.validate_matrix(matrix, j['seqids'])
        if err:
            print(err)
            sys.exit(1)
    usr_tree_str = treeutils.run_fneigbor(phystr)


# Read tree
fh = open('datasets/genome/{}/tree.newick'.format(args.reference))
ref_tree_str = fh.read().strip()
fh.close()


# Compare user and reference tree
res = treeutils.ete_compare(usr_tree_str, ref_tree_str)
qt, rt, nrf, rf, maxrf, src_br, ref_br, treekod = res
qd, nqd = treeutils.quartet_dist(usr_tree_str, ref_tree_str)
print('nRF:     {:.2f}'.format(nrf))
print('RF:      {}'.format(rf))
print('max RF:  {}'.format(maxrf))
print('nQD:     {:.4f}'.format(nqd))





