# How to use the pipeline

1. Clone this repository
2. Make sure you have the following tools available
    - fastqc (available [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - cutadapt (available [here](https://cutadapt.readthedocs.io/en/stable/index.html))
    - R
    - count_constructs (available [here](https://github.com/Luk13Mad/count_constructs))
    - python
3. Change into the scripts folder and run ``./setup.sh``. This creates other needed folders.
4. Load your raw data into the raw_sequencing_data folder.
5. Use snakemake to generate final count tables and corresponding quality control graphics.
6. Use jupyter notebook to run analysis.

# Example

Below is an example of the pipeline in action, step by step.

- Step 1. up to and including step 3. from above have been completed.

- Assume you have the following sequencing data in the raw_sequencing_data folder:  
```
/.../raw_sequencing_data/X204SC20113232-Z01-F010/raw_data/
├── PA_TU8988_A1
│   ├── PA_TU8988_A1_EKDL210008849-1a-D701-N504_HTK7WDSX2_L4_1.fq.gz
│   └── PA_TU8988_A1_EKDL210008849-1a-D701-N504_HTK7WDSX2_L4_2.fq.gz
├── PA_TU8988_B2
│   ├── PA_TU8988_B2_EKDL210008849-1a-D702-N505_HTK7WDSX2_L4_1.fq.gz
│   └── PA_TU8988_B2_EKDL210008849-1a-D702-N505_HTK7WDSX2_L4_2.fq.gz
├── PA_TU8988_C3
│   ├── PA_TU8988_C3_EKDL210008849-1a-DY0088-N506_HTK7WDSX2_L4_1.fq.gz
│   └── PA_TU8988_C3_EKDL210008849-1a-DY0088-N506_HTK7WDSX2_L4_2.fq.gz
├── PA_TU8988_D4
│   ├── PA_TU8988_D4_EKDL210008849-1a-D704-N504_HTK7WDSX2_L4_1.fq.gz
│   └── PA_TU8988_D4_EKDL210008849-1a-D704-N504_HTK7WDSX2_L4_2.fq.gz
├── PA_TU8988_E5
│   ├── PA_TU8988_E5_EKDL210008849-1a-D701-N505_HTK7WDSX2_L4_1.fq.gz
│   └── PA_TU8988_E5_EKDL210008849-1a-D701-N505_HTK7WDSX2_L4_2.fq.gz
└── PA_TU8988_F6
    ├── PA_TU8988_F6_EKDL210008849-1a-D702-N506_HTK7WDSX2_L4_1.fq.gz
    └── PA_TU8988_F6_EKDL210008849-1a-D702-N506_HTK7WDSX2_L4_2.fq.gz

```
In this example you conducted a CRISPR screen with three timepoints TP1, TP2 and TP3 each of which has two replicates. So we have TP1 (PA_TU8988_A1 and PA_TU8988_B2), TP2 (PA_TU8988_C3 and PA_TU8988_D4) and TP3 (PA_TU8988_E5 and PA_TU8988_F6).

- To obtain the raw count table and associated quality control graphic for TP1 replicate 1 (that is PA_TU8988_A1) you have to run the following snakemake command ``snakemake --profile LSF processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/graphics/summary_graphics_PA_TU8988_A1.pdf`` You have to set up the profile yourself so snakemake work with your cluster. This command will prompt snakemake to execute the following sequence of rules:  
```
Building DAG of jobs...
Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
count_constructs             1              1              1
cutadapt_raw_data            1              1              1
fastqc_raw_data              1              1              1
fastqc_trimmed_data          1              1              1
make_summary_graphics        1              1              1
total                        5              1              1

rule fastqc_raw_data:
    output: raw_sequencing_data/FastQC/X204SC20113232-Z01-F010/PA_TU8988_A1/FastQC_raw_data.done
    wildcards: run=X204SC20113232-Z01-F010, type=PA_TU8988_A1

rule cutadapt_raw_data:
    input: raw_sequencing_data/FastQC/X204SC20113232-Z01-F010/PA_TU8988_A1/FastQC_raw_data.done
    output: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_1_trimmed.fq.gz, processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_2_trimmed.fq.gz
    wildcards: run=X204SC20113232-Z01-F010, type=PA_TU8988_A1

rule fastqc_trimmed_data:
    input: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_1_trimmed.fq.gz, processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_2_trimmed.fq.gz
    output: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/FastQC/trimmed_data.done
    wildcards: run=X204SC20113232-Z01-F010, type=PA_TU8988_A1

rule count_constructs:
    input: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/FastQC/trimmed_data.done, processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_1_trimmed.fq.gz, processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/PA_TU8988_A1_L1_2_trimmed.fq.gz
    output: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/count_results_PA_TU8988_A1.tsv
    wildcards: run=X204SC20113232-Z01-F010, type=PA_TU8988_A1

rule make_summary_graphics:
    input: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/count_results_PA_TU8988_A1.tsv
    output: processed_sequencing_data/X204SC20113232-Z01-F010/PA_TU8988_A1/graphics/summary_graphics_PA_TU8988_A1.pdf
    wildcards: run=X204SC20113232-Z01-F010, type=PA_TU8988_A1
```
You obtain the other count tables in the same way. With this the snakemake part is done and we proceed inside the jupyter notebook.

- We have to collect the count tables generated by snakemake into one file. For this we use the notebook section **Collect the count results**. In our example we would set the variables to ``subfolders_with_countresults = ["X204SC20113232-Z01-F010"]
celllines = ["PATU"]`` This will aggregate the count results from the folder X204SC20113232-Z01-F010 using the script aggregate_counts.py. You have the option to normalize counts to counts per million (CPM). 

- Now that we have collected the count files (and optinally normalized them as CPM), we have to merge them into one file. We do this with the script merge_counts.py.

- The next few cells in the notebook are optional. The columns in the merged count tables are named after the files they originated from, if you want to rename them for example to something shorter, now is the time.

- Now we have to generate TOML files with some options for later processing. One TOML file per cell line, placed right next to the merged_counts file. In this example we need only one TOML file.  
```
cellline = "PATU"

[prepareLFC]
earlyTP = "PATUTP1"
comparisons = [[1,2],[1,3],[2,3]]

[GEMINI]
sample_anno_rows = [14,15,18,19] #comma separated list of row numbers needed from sample_annotation
counts_cols = ["PATUTP1","PATUTP3"] #comma separated list of column names needed from counts
earlyTP = ["PATUTP1.REP1","PATUTP1.REP2"] # comma separated list of column names with ETP from counts
ETP_col_index = [1,6] #column indices with early timepoint columns in count
LTP_col_index = [5,2] #column indices with late timepoint columns in count
mean_x = 1.0 #mean for x variable for GEMINI
sd_x = 1.0 #sd for x variable for GEMINI
mean_y = 0.0 #mean for y variable for GEMINI
sd_y = 0.5 #sd for y variable for GEMINI
mean_xx = 1.0
sd_xx = 0.5
mean_s = 0.0 #mean for s variable for GEMINI
sd_s = 0.5 #mean for r variable for GEMINI
pc_genes = ["AQR","CDK9","PABPN1"]
```


Tested with fastqc 0.11.5, cutadapt 1.18, R 4.1.0 (with reshape2, dplyr, ggrepel, ggplot2, gemini, configr), count_constructs 1.0 and python 3.10.11 (with snakemake_7.25.3, pandas, numpy, click_8.1.3, statsmodels, toml, scipy).
