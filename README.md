# HAM
Measuring Annotated Homology
While annotated protein sequences are widely used in machine learning applications, pre-processing these sequences regarding homology is mainly limited to clustering complete sequences based on global alignment without considering their annotations. Here, I am introducing new tools that identify all possible local homologies between annotated sequences within the same or across two datasets and then resolve these homologies.

#### Please reference the following preprint:

Malhis N. Pre-processing annotated homologous regions in protein sequences concerning machine-learning applications" *bioRxiv* (2024). [doi.org/10.1101/2024.10.25.620288] (https://doi.org/10.1101/2024.10.25.620288).  

## Minimum Hardware Requirements

OS: Linux (tested on Ubuntu).

RAM: 8 GB minimum, 16 GB recommended.

CPU: Multicore with 4+ cores recommended.

## To install:

```bash
# Clone the HAM software:	
git clone https://github.com/NawarMalhis/HAM.git
# Clone the "annotated fasta format" library:	
git clone https://github.com/NawarMalhis/AFF.git
# Change directory:	
cd HAM
# Update the path to the AFF (annotated fasta format) folder in ham_param.py
aff_path = '/xxx/xxx/AFF/'
# Create a conda ham environment:
conda env create -f ham.yml
```

## To run:

```bash
# Activate the ham_env environment:
conda activate ham_env
```

## Input data:
Input data should be in an annotated fasta format such that each sequence is annotated with a single line of annotation; for each sequence, we need a fasta header line with a unique accession, a sequence line, and an annotation line. Annotations can include '1', '0', and '-'. Where '-'is used as a mask that is neither '1' nor '0'. The annotated fasta format file can start with lines marked with '#' as comments. Example of a sequence in an annotated fasta format with a single annotation line:
```bash
>Q86FP9
MKHFAILILAVVASAVVMAYPERDSAKEGNQEQERALHVKVQKRTDGDADYDEYEEDGTTPTPDPTAPTAKPRLRGNKP
000000000000000000000000000000000---------0000000000000000000000000001111100000
```
All input files (datasets) need to be in the same data directory. Our input files TR2008.af and TS2008.af are in the data directory' data/'. 
The data directory can be located anywhere.

## Optional parameters for ham.py and hac.py:

    • Identity cut off, default 80%:	-ico 80	
    • Minimum aligned size, default 10:	-msz	10
    • Number of threads, default 8:	-num_threads 8

## Example:
First, we run ham.py for each of our training/testing datasets to identify shared homologous regions.
```bash
(ham_env) ~/Tools/HAM$ python3 ham.py -in1 TS2008.af -in2 TR2008.af -p ~/data/
```
A results directory is created inside our data directory, and three files are added:
    • ham-details-TS2008-TR2008.tsv: includes a list of the one-to-one residue homology between the two input files.
    • TR2008-homology-to-TS2008.af: this is the same TR2008.af input file with three extra annotation lines added, H0, shows the ‘0’ annotations of homologous residues in TS2008.af. 'H1' shows the '1' annotations of homologous residues, and 'H-'shows the '-'annotations of homologous residues. Example line:
```bash
>PDB:1a3b_I
ITYTDCTESGQDLCLCEGSDVCGKGNKCILGSNGEENQCVTGEGTPKPQSHNDGDFEEIPEEYLQ
00000000000000000000000000000000000000000000000000000111111111110
.................000000000000....................................
........................................................111111...
..........................--------...............................
```
    • TS2008-homology-to-TR2008.af: just like with TR2008-homology-to-TS2008.af. This file identifies regions in TS2008.af that are homologous to those in TR2008.af.






