## Installation
**TARPMIR REQUIRES PYTHON 3.8**

```sh
$ git clone https://github.com/nntavantzis/ptixiaki.git
$ cd ptixiaki
$ python3 -m venv venv
$ source .venv/bin/activate
$ pip install -r requirements.txt
```

## Usage
*To predict disrupted miRNAs:*

```sh
$ source .venv/bin/activate
$ python3 main.py -i input_file

INPUT FILE NEEDS TO BE A VEP OUTPUT FILE (.tsv)
```


*For a single TarPmiR prediction:*
```sh
$ source .venv/bin/activate
$ python3 runTarPmiR.py -mrna mrna_file [-mirna mirna_file] [-pcut probability_cut] [-num number_of_processes]
```

The main program outputs 2 files, the modified VEP input (default `output/output_disrupted.tsv`) and the clean output that contains non-duplicate gene-miRNA pairs (default `output/clean_output_disrupted.tsv`).