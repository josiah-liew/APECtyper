![version](https://img.shields.io/badge/version-1.0.0-blue)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-3.0.en.html)
[![made-with-bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)

# APECtyper

*In silico* command line tool for typing Avian Pathogenic *Escherichia coli* (APEC)
 
## Description

`APECtyper` is a Bash shell script that wraps several well-establised tools into a single pipeline in order to classify Avian Pathogenic *Escherichia coli* (APEC) based on the revised APEC pathotyping scheme developed by [Johnson et al. *in prep*](). 


([ECTyper](), [mlst], [BLAST+], and [R]) 

## Dependencies

* [ECTyper](https://github.com/phac-nml/ecoli_serotyping)
* [mlst](https://github.com/tseemann/mlst)
* [NCBI BLAST+ blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/?report=reader&%2F%3Freport=reader)
    * This should already be installed as an mlst dependency. 
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

```
% ./APECtyper.sh
Usage: APECtyper.sh [OPTIONS] -f [FASTA or DIR] -o [DIR]
	-h		print this usage message
	-v		print the version
	-r		print citation
	-f		FASTA contig file or directory containing multiple FASTA files
	-o		output directory
        -i              minimum blast % identity [default: 90]
	-c              minimum blast % coverage [default: 90]
	-t              number of threads to use [default: 1]
	-s              combine reports from multiple samples into single TSV file
```

### Example Usage

Single FASTA file:
```
% ./APECtyper.sh -f data/assemblies/sample1.fasta -o sample1_output
```

Directory containing multiple FASTA files:
```
% ./APECtyper.sh -f data/assemblies -o multisample_output
```

## Output

Within the user-defined output directory, there will be three directories – `serotype`, `mlst`, and `blast` – containing the individual output files from ECTyper, mlst, and blastn, respectively. There will also be two TSV files per sample: 

* `pathotype_results_SAMPLE.tsv`: A tab-separated file summarizing the outputs from ECtyper and mlst, as well as the pathotype classification.
* `blast_results_SAMPLE.tsv`: A tab-separated summary of the blastn results against the custom APEC virulence gene database. Note that the percent identity and percent coverage thresholds can be set by the user with the `-i` and `-c` flags.  

If the `-s` flag is included when running `APECtyper.sh`, there will be two additional files that are compilations of all blast reports and all pathotype reports: `blast_results_summary.tsv` and `pathotype_results_summary.tsv`.

The `pathotype_results_SAMPLE.tsv` file has the following columns:  

Column | Example | Description
-------|---------|------------
Sample | `sample1` | Name of input assembly file with file extension removed
Species | `Escherichia coli` | Species identified by ECTyper
Serotype | `O25:H4` | Serotype identified by ECTyper
SerotypeQC | `PASS (REPORTABLE)` | QC message produced by ECTyper (See [ECTyper page](https://github.com/phac-nml/ecoli_serotyping#quality-control-qc-module) for a description of all possible QC codes) 
ST | `131` | Sequence type identified by mlst
APEC_plasmid | `Present` | Whether or not APEC plasmid markers *hlyF* and *ompT* were found
Pathotype | `High Risk APEC` | Pathotype classification (possible values include `High Risk APEC`, `APEC`, `High Risk non-APEC`, `non-APEC`, or `Not E. coli`)

The `blast_results_SAMPLE.tsv` file has the following columns:

Column | Example | Description
-------|---------|------------
Sample | `sample1` | Name of input assembly file with file extension removed














## Citations

```
./APECtyper.sh -r
```

If you use APECtyper in your work, please cite:  






Please also cite:

* ECtyper
* mlst










