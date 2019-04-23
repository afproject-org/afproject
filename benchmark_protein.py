import argparse
import json
import sys

import numpy as np

from libs import utils
from libs import protein

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', dest="input_file", required=True,
                    type=argparse.FileType('r'), metavar="FILE",
                    help="Input file")
parser.add_argument('-f', '--format', dest="format", choices=['tsv', 'phylip'],
                    help="File format", default='tsv')
parser.add_argument('-r', '--reference', dest="reference", 
                    choices=['low-ident', 'high-ident'],
                    help="File format", default='tsv')
args = parser.parse_args()


# Read reference
fh = open('datasets/protein/{}/dataset.json'.format(args.reference))
DATASET = json.load(fh)
fh.close()
fh = open('datasets/protein/{}/ids.json'.format(args.reference))
IDS = json.load(fh)
fh.close()
REF_SEQ_IDS = {}
for sid, sstr in zip(IDS["seqids"], IDS["seqstr"]):
    REF_SEQ_IDS[sid] = tuple(sstr.strip().split('.'))
strlev_queryset = [l["name"] for l in DATASET["levels"]]

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



err, clean_data = protein.validate_dpair(dpair, IDS['seqids'])
if err:
    print(err)
    sys.exit(1)
dpair = utils.read_clean_tsv_data_dpair(clean_data.split('\n'))


benchmark = protein.benchmark(
            dpair, 
            REF_SEQ_IDS,
            strlev_queryset,
            1000
)

l = []
for strlevel, roc_auc, tpr_str, fpr_str in benchmark:
    print("AUC\t{}\t{}".format(strlevel, roc_auc))
    l.append(roc_auc)

print("Mean AUC: {}".format(np.mean(l)))
print("Stdev AUC: {}".format(np.std(l)))
