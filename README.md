# RAT-DNA-screening

Code for data processing, software methods, and pipelines related to the manuscript 
"Random adversarial threshold search enables automated DNA screening"

This code has been validated on Mac OS Ventura 13.3.1 with 8-core Apple silicon (M1 chip).

## Installation

Installation time is about 5 minutes.

First, [install Python](https://www.python.org/about/gettingstarted/)
[version 3.7](https://www.python.org/downloads/release/python-373/) or higher.

Install virtualenv, create a new Python virtual environment (version 3.7.3 has been
validated) and activate it:

```
pipx install virtualenv # or, `pip install virtualenv`. restart shell afterward.
git clone https://github.com/dgretton/RAT-DNA-screening.git
cd RAT-DNA-screening
virtualenv venv
```

Then, install dependencies and run.

```
pip install --upgrade pip
pip install -r requirements.txt
python mhv_parallelized.py
```

Run time is about 2 minutes. Sample databases generated using the Metropolis-Hastings
(MH) algorithm appear in resources/. Due to security concerns, the exact function
proportional to the probability density we utilized for MH is replaced with randomness
in funtrp_blosum_predict.py, but the methodology is identical.

To run experimental-data-dependent code, download and place data according to instructions
below. Then, for example:

```
python rat_heatmap.py
```

Extended Data Figures 6 and 7 appear.

## Data

Metropolis-Hastings sampler does not require data download. For other functions:

"2021-03-29 compiled_library_counts.csv" should be placed in 
"resources/NGS/2021_02_08_multiple_libraries/2021-03-29 compiled_library_counts.csv"

"PQSVEC_variant_counts.csv" should be placed in
"resources/NGS/2020_12_10_synthesis_screening_lib_12/PQSVEC_variant_counts.csv"

"screening_databases.zip" should be decompressed into "resources/", creating
"resources/screening_databases/"

## resources/NGS/2021_02_08_multiple_libraries

Window PQSVECRPFVFGAGKPYEF was studied first as a pilot. The data processing pipeline
was as follows.

The following fastq-format files containing Next-Generation Sequencing (NGS) reads for
window PQSVECRPFVFGAGKPYEF were saved in the directory raw_data.

- Sample point 1: 201202Esv_D20-6637_1_sequence.fastq. Pre-electroporation

- Sample point 2: 201202Esv_D20-6638_1_sequence.fastq. Post-electroporation, pre-extrusion
of initial phage

- Sample point 3: 201202Esv_D20-6639_1_sequence.fastq. Post-extrusion, pre-infection

- Sample point 4: 201202Esv_D20-6640_1_sequence.fastq. Post-infection (one full round of
selection)

- Sample point 5: 201202Esv_D20-6641_1_sequence.fastq. Two full rounds of selection (not
collected for any windows other than this pilot, not used in this study)

- Sample point 6: 201202Esv_D20-6642_1_sequence.fastq. Three full rounds of selection (not
collected for any windows other than this pilot, not used in this study)

These files may be found in the archive all_raw_m13_ngs.zip (7.5GB, decompresses to 40GB)
along with raw NGS data for the 8 other windows.

The summary of variant counts across all sample points may be found in
PQSVEC_variant_counts.csv.

**process_fastq.py**

Converts the .fastq files into two different csv files.

-> *_raw_sequences.csv contains a list of all the sequences (same as fastq) with no
processing

-> *_sequencing_counts.csv contains a list of all the UNIQUE sequences and how many times
they occurred

**count_lib_members_in_data.py**

Takes the #_sequencing_counts.csv, and the csv file containing all of the library members,
and computes the number of times each library member was observed. It outputs this to
*_counted.csv

**compile_counted.py**

Process each of the libraries separately just because they're slow; this file compiles
all the results into a single csv that lists the number of counts of each library member
in each sample

**graphs.py**

Processes the results of Brian_all_counted.csv and produces analysis of library coverage.
