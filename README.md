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
# Create a conda ham environment:
conda env create -f ham.yml
```

## To run:
Activate the ham_env environment:
```bash
conda activate ham_env
```
