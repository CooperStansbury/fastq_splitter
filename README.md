# fastq_splitter
A tool to split pore-c reads into pairwise reads based on a restriction fragment. 

## Usage:

First, create the environment:

```
conda env create
conda activate splitq
```

Then run the command: 

```
python splitq/ -i <input fastq> -o <output directory> -n <number of sourrounding base pairs> -overwrite <overwrite destination>
```

Input arguments are:

```
usage: [-h] [-i I] [-n [N]] [-o [O]] [-overwrite [OVERWRITE]]

A Python3 commandline tool to virtually digest sequences based on NlaIII sites.

optional arguments:
  -h, --help            show this help message and exit
  -i I                  Input merged fastq file.
  -n [N]                Number of base pairs surrounding a cut site.
  -o [O]                The path to a directory for all output files.
  -overwrite [OVERWRITE]
                        If True, overwrite the directory at `-o`
```