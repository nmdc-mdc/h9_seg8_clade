
# H9 Influenza Sequence Clade and Genotype Classification

## Overview

This software is designed for classifying H9 influenza sequences into specific clades or genotypes, providing valuable insights for H9 influenza research and epidemiological analysis.

## Environment Setup Guide

This guide provides step-by-step instructions to set up essential bioinformatics tools using Anaconda.

### Step 1: Download and install Anaconda

1. Visit the [Anaconda](https://www.anaconda.com/download) to find the latest version suitable for your operating system.
2. Download and install Anaconda by following the instructions on the website.

Download Anaconda
For Linux users, you can download the latest version of the Anaconda installer script via the command line using the wget command. Please note that the version number in the command below (e.g., latest-version) may need to be adjusted according to the latest released version.
  ```bash
   wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh 
  ```
Install Anaconda3
  ```bash
   bash Anaconda3-2024.10-1-Linux-x86_64.sh  
  ```
Initialize Anaconda
  ```bash
  export PATH="/Path/anaconda3/bin:$PATH"
  source ~/.bashrc  
  ```
### Step 2: Install Bioinformatics Tools

#### Option 1: Using setup.sh

Execute the setup.sh script to set up the environment.

  ```bash[README.md](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2FGit%2Finflu_seg8%2FREADME.md)
  bash setup.sh 
  ```
#### Option 2: Using environment.yml
If the setup.sh script does not work as expected, you can set up the environment using the environment.yml file.
  ```bash
  conda env create -f environment.yml 
  ```
#### Option 3: Manual installation guide for tools
If the setup.sh or environment.yml files are not working for you, you can manually install all the required tools. Follow the steps below to set up your environment and install the necessary tools.

1. **Install Python 3.11**

    ```bash
   conda create --name h9-env python=3.11
   ```

2. **Install tools from conda**
   ```bash
   conda activate h9-env
   conda install pandas
   conda install -c conda-forge biopython=1.84
   conda install bioconda::mafft
   ```
3. **Install blast**

BLAST version 2.16 is required. Since conda does not support direct installation of version 2.16, please download the Linux version of BLAST 2.16 from the [NCBI FTP server](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/).

1. Extract the package:
   ```bash
   tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
   ```

2. Add Blast to h9-env environment
   ```bash
   mv "Path_to/ncbi-blast-2.16.0+/bin/"* "Path_to/anaconda3/envs/h9-env/bin/"
   export H9_BLAST_PATH="Path_to/anaconda3/envs/h9-env/bin/blastn"
   ```
3. Verify Installation
    ```bash
   blastn -version
   ```
## Script Execution

### Step 1: Activate the h9-env Environment

  ```bash
  conda activate h9-env 
  ```

### Step 2: Run the Script

#### Mode 1: Segment Clade Classification Only

If you want to classify by segment clade only, use the following command:

```bash
python "h9_clade_seg8.py" -i "path_to/input_file" -n none -o "path_to/output_directory"
```

**Example:**

```bash
python "h9_clade_seg8.py" -i "example/example.fasta" -n none -o "output"
```

#### Mode 2: Segment Clade and Genotype Classification

If you want to classify by both segment clade and genotype, use this command:

```bash
python "h9_clade_seg8.py" -i "path_to/input_file(s)" -n "isolate_name(s)" -o "path_to/output_directory"
```

**Example:**

```bash
python "h9_clade_seg8.py" -i "example/isolate_1.fasta,example/isolate_2.fasta" -n "isolate_1,isolate_2" -o "output"
```

In this example, the input consists of multiple isolates separated by commas. The `-n` option specifies the names of each isolate in the same order as the input files.

---

## Output

After running the classification script, navigate to the `path_to/output_directory` directory to find the results.

### Mode 1: Segment Clade Classification Only

For classification by segment clade only, the output will include:
- A directory named `none`.
- Two JSON files:
  - `cladeResultJson.json` — contains clade classification results.
  - `no_pass_number.json` — lists any sequences that did not pass the classification criteria.

Inside the `none` directory, you will find:
- Three subdirectories:
  - `feature` — containing extracted features variation of the segment sequence from our database.
  - `meta` — containing metadata for the sequences that close to input sequences from our data.
  - `tree` — holding the phylogenetic tree data for visualizing clade relationships.
- A `cladeResultJson.json` file specific to this classification.

### Mode 2: Segment Clade and Genotype Classification

For classification by both segment clade and genotype, the output will include:
- One or more directories, each named after the corresponding isolate(s) provided in the input.
- Two JSON files:
  - `cladeResultJson.json` — containing both clade and genotype classification results for all input sequences.
  - `no_pass_number.json` — listing sequences that did not pass the classification criteria.

Within each isolate-named directory, you will find:
- Three subdirectories:
  - `feature` — containing extracted features variation of the segment sequence from our database.
  - `meta` — containing metadata for the sequences that close to input sequences from our data.
  - `tree` — storing phylogenetic tree data for clade and genotype visualization.
- A `cladeResultJson.json` file specific to that isolate’s classification results. 

---
