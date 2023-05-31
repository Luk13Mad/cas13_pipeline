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
/.../raw_sequencing_data/X204SC20113232-Z01-F010/raw_data
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

Tested with fastqc 0.11.5, cutadapt 1.18, R 4.1.0 (with reshape2, dplyr, ggrepel, ggplot2, gemini, configr), count_constructs 1.0, python 3.10.11 (with snakemake_7.25.3, pandas, numpy, click_8.1.3, statsmodels, toml, scipy)
