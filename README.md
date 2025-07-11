# Homology-based Annotation Masking (HAM) and Homology Annotation Conflict (HAC)

Accurate preprocessing of annotated protein sequences with regard to homologies is essential for maintaining the integrity of machine-learning applications. This study presents two new tools—HAM and HAC—developed to address key challenges, including data leakage and annotation conflicts. HAM detects local homologous regions across datasets and applies annotation masking to prevent contamination between training and test sets. HAC identifies local homologies within a single dataset, exposing and resolving conflicting annotations. Applied to three benchmark datasets, these tools uncover significant overlaps and annotation inconsistencies, underscoring the need to integrate them into standard preprocessing pipelines to improve the reliability of models predicting intrinsically disordered protein binding sites.

**Please reference the following preprint:**

Malhis N. Pre-processing homologous regions in annotated protein sequences concerning machine-learning applications" *bioRxiv* (2024). https://doi.org/10.1101/2024.10.25.620288   

### Minimum Hardware Requirements

OS: Linux (tested on Ubuntu).

RAM: 8 GB minimum, 16 GB recommended.

CPU: Multicore with 4+ cores recommended.

### To install:

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

### To run:

```bash
# Activate the ham_env environment:
conda activate ham_env
```

### Input data:
Input data should be in a simple annotated FASTA format, where each sequence is annotated with a single line of annotation. For each sequence, we require a FASTA header line with a unique accession, a sequence line, and an annotation line. Annotations can include '1', '0', and '-'. Where '-' is used as a mask that is neither '1' nor '0'. The simple annotated fasta format file can start with lines marked with '#' as comments. Example of a sequence in a simple annotated fasta format:
```bash
>Q86FP9
MKHFAILILAVVASAVVMAYPERDSAKEGNQEQERALHVKVQKRTDGDADYDEYEEDGTTPTPDPTAPTAKPRLRGNKP
000000000000000000000000000000000---------0000000000000000000000000001111100000
```
All input files (datasets) need to be in the same data directory. Our input files TR2008.af and TS2008.af are in the data directory' data/'. 
The data directory can be located anywhere.

### Optional parameters for ham.py and hac.py:

    • Identity cut off, default 80%:	-ico 80	
    • Minimum aligned size, default 10:	-msz	10
    • Number of threads, default 8:	-num_threads 8

## HAM: Homology-based Annotation Masking
This includes two tools, ham.py and ham_mask_homologous.py. The first, ham.py, identifies homologous regions between the two input files. The second, ham_mask_homologous.py, masks those homologous regions identified by ham.py.

### HAM Example:
First, we run ham.py for each of our training/testing datasets to identify shared homologous regions.
```bash
(ham_env) ~/Tools/HAM$ python3 ham.py -in1 TS2008.saf -in2 TR2008.saf -p ~/data/
```
A results directory is created inside our data directory, and three files are added:
1. ham-details-TS2008-TR2008.tsv: includes a list of the one-to-one residue homology between the two input files.
2. TR2008-homology-to-TS2008.af: this is the same TR2008.saf input file with three extra annotation lines added, H0, shows the ‘0’ annotations of homologous residues in TS2008.saf. 'H1' shows the '1' annotations of homologous residues, and 'H-'shows the '-' annotations of homologous residues. Example line:
```bash
>PDB:1a3b_I
ITYTDCTESGQDLCLCEGSDVCGKGNKCILGSNGEENQCVTGEGTPKPQSHNDGDFEEIPEEYLQ
00000000000000000000000000000000000000000000000000000111111111110
.................000000000000....................................
........................................................111111...
..........................--------...............................
```
3. TS2008-homology-to-TR2008.af: just like with TR2008-homology-to-TS2008.af. This file identifies regions in TS2008.af that are homologous to those in TR2008.saf.

Note that we can run ham multiple times against different files:
```bash
(ham_env) ~/Tools/HAM$ python3 ham.py -in1 TS2008.af -in2 xxxx.saf -p ~/data/
(ham_env) ~/Tools/HAM$ python3 ham.py -in1 TS2008.af -in2 yyyy.saf -p ~/data/
```
Then, we run ham_mask_homologous.py to mask homologous residues identified by all previous runs:
```bash
(ham_env) ~/Tools/HAM$ python3 ham_mask_homologous.py -in TS2008.saf -p data/
```
ham_mask_homologous reads TS2008.saf from the data directory, then searches for all files in the results directory that start with "TS2008-homology-to-", and processes them one at a time by masking TS2008.saf residues that are homologous. In this example, ham_mask_homologous will consider TS2008 Homology to TR2008.af, TS2008-homology-to-xxxx.af, and TS2008-homology-to-yyyy.af.
The final masked dataset is saved in the data directory as ham-masked-TS2008.af

## HAC: Homology Annotation Conflict
This includes two tools, hac.py and hac_resolve_conflict.py. The first, hac.py, identifies homologous regions with conflicting annotations between the sequences of the input file. The second, hac_resolve_conflict.py, enables us to resolve these conflicting annotations identified by hac.py by reannotating them with either '1', '0', or '-' (masking).

### HAC Example:
Our input file, TR2008.saf, is located in the data directory 'data/'. First, we run hac.py for each input file:

```bash
(ham_env) ~/Tools/HAM$ python3 hac.py -in TR2008.saf -p data
```
This generates two files: 
    • hac-details-TR2008.tsv: includes a list of the one-to-one residue homology between the sequences of the input file.
    • TR2008-annotation-conflict.af: the same TR2008.saf input file with two extra annotation lines added, H0 shows the '0' annotations of the homologous residues in TR2008.saf, and 'H1' shows the '1' annotations. Example sequence with annotations in the TR2008-annotation-conflict.af file:
```bash
>PDB:1a3b_I
ITYTDCTESGQDLCLCEGSDVCGKGNKCILGSNGEENQCVTGEGTPKPQSHNDGDFEEIPEEYLQ
00000000000000000000000000000000000000000000000000000111111111110
000000000000000000000000000000000000000000000000000000......00000
..................................................111111111111111
```
To resolve the conflict, we can choose one of three priorities: '01', '10', and '-'.
```bash
(ham_env) ~/Tools/HAM$ python3 hac_resolve_conflict.py -in TR2008.saf -p data/ -pr '01'
```

'01' converts conflicting annotations of '0' and '1' into '1'. Thus, the updated annotation to the above sequence is:
```bash
00000000000000000000000000000000000000000000000000111111111111111
```
The final dataset with resolved annotations is saved in the data directory as TR2008-resolved-01.af

'10' converts conflicting annotations of '0' and '1' into '0'. Thus, the updated annotation to the above sequence is:
```bash
00000000000000000000000000000000000000000000000000000011111100000
```
The final dataset with resolved annotations is saved in the data directory as TR2008-resolved-10.af

'-' converts conflicting annotations of '0' and '1' into '-'. Thus, the updated annotation to the above sequence is:
```bash
00000000000000000000000000000000000000000000000000----111111-----
```
The final dataset with resolved annotations is saved in the data directory as TR2008-resolved-masked.af

