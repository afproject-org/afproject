# afproject
Codebase that benchmarks alignment-free sequence comparison methods. The code is used in publicly available webservice AFproject (http://afproject.org). 

This repository contains Python scripts to benchmark an alignment-free methods offline.

## Requirements
The following software tools are required:

1. [fneighbor](http://emboss.sourceforge.net/apps/cvs/embassy/phylipnew/fneighbor.html) from the EMBOSS package
2. [tqdist](http://users-cs.au.dk/cstorm/software/tqdist/)

The following Python packages are required:

1. [numpy](http://www.numpy.org/) for general calculations
2. [scikit-learn](https://scikit-learn.org/stable/) for calculating ROC curves and AUC.
3. [ete3](http://etetoolkit.org/) for calculating Robinson-Foulds distance

## Usage

To calculate accuracy with `fish_mito` dataset from the genome-based phylogeny category for input file in `TSV` format, just run:

```bash
python benchmark_genome.py --input example_input/genome/std/assembled-fish_mito_random.tsv --format tsv --reference fish_mito
```

You will get the information about normalized Robinson-Foulds (nRF) distance and normalized Quartet Distance (nQD).

```bash
nRF:     1.00
RF:      44.0
max RF:  44.0
nQD:     0.6897
```

You can provide your input files in TSV / PHYLIP / Newick format.

```bash
python benchmark_genome.py --input example_input/genome/std/assembled-fish_mito_random.phy --format phylip --reference fish_mito
```

```bash
python benchmark_genome.py --input example_input/genome/std/assembled-fish_mito_random.newick --format phylip --reference fish_mito
```

To calculate performance in the CRM dataset for your input PHYLIP file, just run:

```bash
python benchmark_genreg_crm.py --input example_input/genreg/crm_random.phy --format phylip
```

You will get information about your method's performance across all 7 tissues:

```
Tissue	k	n	percent
FB	151	300	50.3
FP	127	253	50.2
FT	17	36	47.2
HM	149	300	49.7
HL	17	36	47.2
HH	72	136	52.9
FE	68	136	50.0

Weighted average:	50.21
Standard deviat.:	1.83
Average:         	49.65
```

To calculate performance in the protein classification category, jus run:

```bash
python benchmark_protein.py -i example_input/protein/low-ident/low-ident_random.tsv -f tsv -r low-ident
```

You will get information about your method's performance across 4 structural levels:

```
AUC	class	0.5005448508889172
AUC	fold	0.4998414124679509
AUC	superfamily	0.5042921080665331
AUC	family	0.5021884324506269
Mean AUC: 0.5017167009685071
Stdev AUC: 0.0017135629461110453
```


To calculate performance in the gene tree category, jus run:

```bash
python benchmark_genetree.py -i example_input/genetree/swisstree_random.tsv -f tsv
```

You will get information about your method's performance across 4 structural levels:

```
ST001	nRF	0.97
ST002	nRF	1.0
ST003	nRF	1.0
ST004	nRF	1.0
ST005	nRF	0.96
ST007	nRF	1.0
ST008	nRF	1.0
ST009	nRF	1.0
ST010	nRF	1.0
ST011	nRF	0.99
ST012	nRF	1.0
MEAN nRF: 0.993
STDEV nRF: 0.0124
```

