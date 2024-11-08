
# H9 Influenza Sequence Clade and Genotype Classification

## Overview

This software is designed for classifying H9 influenza sequences into specific clades or genotypes, providing valuable insights for H9 influenza research and epidemiological analysis.

## Environment Setup Guide

This guide provides step-by-step instructions to set up essential bioinformatics tools using Miniconda.

### Step 1: Download Miniconda

1. Visit the [Miniconda Documentation](https://docs.anaconda.com/free/miniconda/) to find the latest version suitable for your operating system.
2. Download and install Miniconda by following the instructions on the website.

### Step 2: Install Bioinformatics Tools

Using `conda`, you can easily install bioinformatics tools from the Bioconda channel. After each installation, create a symbolic link to make the tool accessible from any directory.

1. **Install BLAST**
   ```bash
   conda install bioconda::blast
   ln -s /root/miniconda3/bin/blastn /bin/
   ```

2. **Install IQ-TREE**
   ```bash
   conda install bioconda::iqtree
   ln -s /root/miniconda3/bin/iqtree /bin/
   ```

3. **Install MAFFT**
   ```bash
   conda install bioconda::mafft
   ln -s /root/miniconda3/bin/mafft /bin/
   ```

### Step 3: Install Python Modules

For data analysis and bioinformatics scripting, install the following Python modules:

1. **Pandas**
   ```bash
   conda install pandas
   ```

2. **Biopython**
   ```bash
   conda install biopython
   ```

3. **Bio**
   ```bash
   conda install Bio
   ```

### Verification

After installation, verify that each tool is correctly set up by running its command:

```bash
blastn -h
iqtree -h
mafft -h
```

If the commands run successfully and display help options, your setup is complete.

---

## Script Run

### Step 1: Set Up Paths

1. Open the Python script located at `h9_seg8_clade/script/h9_clade_seg8.py`.
2. Go to the `if __name__ == "__main__":` section near the end of the script.
3. Set the paths as follows:

   ```python
   work_dir = "E:"                   # Replace with the path to the h9_seg8_clade directory
   python = "/usr/bin/python3.9"      # Path to Python executable
   blastn = "/usr/bin/blastn"         # Path to BLASTN executable
   mafft = "/usr/bin/mafft"           # Path to MAFFT executable
   iqtree = "/usr/bin/iqtree"         # Path to IQ-TREE executable
   ```

### Step 2: Run the Script

#### Mode 1: Segment Clade Classification Only

If you want to classify by segment clade only, use the following command:

```bash
python "path_to/h9_seg8_clade/script/h9_clade_seg8.py" -i "path_to/input_file" -n none
```

**Example:**

```bash
python "E:/h9_seg8_clade/script/h9_clade_seg8.py" -i example.fasta -n none
```

#### Mode 2: Segment Clade and Genotype Classification

If you want to classify by both segment clade and genotype, use this command:

```bash
python "path_to/h9_seg8_clade/script/h9_clade_seg8.py" -i "path_to/input_file(s)" -n "isolate_name(s)"
```

**Example:**

```bash
python "E:/h9_seg8_clade/script/h9_clade_seg8.py" -i "isolate_1.fasta,isolate_2.fasta" -n "isolate_1,isolate_2"
```

In this example, the input consists of multiple isolates separated by commas. The `-n` option specifies the names of each isolate in the same order as the input files.

---

## Output

After running the classification script, navigate to the `path_to/h9_seg8_clade/output/` directory to find the results.

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
