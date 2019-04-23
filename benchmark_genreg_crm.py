import argparse
import json
import sys

from libs import crm
from libs import utils

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', dest="input_file", required=True,
                    type=argparse.FileType('r'), metavar="FILE",
                    help="Input file")
parser.add_argument('-f', '--format', dest="format", choices=['tsv', 'phylip'],
                    help="File format", default='tsv')
args = parser.parse_args()


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


# Read reference
fh = open('datasets/genreg/crm/ids.json')
ids_json = json.load(fh)
fh.close()
D = crm.read_as_dict(ids_json)

# See if input file is complete
err, clean_data = crm.validate_dpair(dpair, D)
if err:
    print(err)
    sys.exit(1)
d = utils.read_clean_tsv_data_dpair(clean_data.split('\n'))

#Benchmark
ps = []
ns = []
print('Tissue\tk\tn\tpercent')
for tissue, k, n, percent, line in crm.benchmark(d, D):
    print('{}\t{}\t{}\t{:.1f}'.format(tissue, k, n, percent))
    ps.append(percent)
    ns.append(n)

# Summary stats
average, weighted_average, std = crm.stats(ps, ns)
print()
print('Weighted average:\t{:.2f}'.format(weighted_average))
print('Standard deviat.:\t{:.2f}'.format(std))
print('Average:         \t{:.2f}'.format(average))