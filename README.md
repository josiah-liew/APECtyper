![version](https://img.shields.io/badge/version-1.0.0-blue)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-3.0.en.html)
[![made-with-bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)

# APECtyper

*In silico* command line tool for typing Avian Pathogenic *Escherichia coli*
 
## Description

`APECtyper` is a Bash shell script that wraps several well-establised tools ([ECTyper](), [mlst], [BLAST+], and [R]) into a single pipeline to classify  

This tool is the *in silico* version of the *in vitro* multiplex PCR assays developed by [Johnston et al. *in prep*]().

## Requirements

* [ECTyper](https://github.com/phac-nml/ecoli_serotyping)
* [mlst](https://github.com/tseemann/mlst)
* [NCBI BLAST+ blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/?report=reader&%2F%3Freport=reader)
    * This should be installed as an mlst dependency 
* [R](https://cran.r-project.org) version XXXXXXX or higher 

### Installing dependencies via [Conda](https://bioconda.github.io/user/install.html)

```
conda install -c bioconda -c r ectyper mlst r
```

## Installation

This will install the latest version directly from GitHub:

```
git clone https://github.com/JohnsonSingerLab/APECtyper.git
```

Change permissions to make `APECtyper.sh` executable:
```
cd APECtyper
chmod +x APECtyper.sh
```

### Check Installation

Ensure the desired APECtyper version is installed:

```
./APECtyper.sh -v
```

## Usage

### Input Requirements

* assembly in FASTA format or directory containing multiple assemblies in FASTA format; assemblies can be in multiple contigs
* directory to output the results into

### Command Line Options

% ./APECtyper.sh
Usage: APECtyper.sh [OPTIONS] -i [FASTA or DIR] -o [DIR]
	-h		print this message
	-v		print the version
	-c		check SeqKit is in path
	-i		fasta contigs file or directory containing multiple files
	-o		output directory
  -r    prints citation
```

### Example Usage

Single FASTA file:
```
% ./APECtyper.sh -i data/assemblies/BS448.fasta -o example_output
```

Directory containing multiple FASTA files:
```
% ./APECtyper.sh -i data/assemblies -o example_output
```

### Output

Within the user-defined output directory, there will be XXXXXXXXX items:



## Citation

```
./APECtyper.sh -r
```
If you use APECtyper in your work, please cite:  




