# How to use the pipeline

1. Clone this repository
2. Make sure you have the following tools available
    - fastqc (available [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - cutadapt (available [here](https://cutadapt.readthedocs.io/en/stable/index.html))
    - R
    - count_constructs (available [here](https://github.com/Luk13Mad/count_constructs))
    - python
3. Change into the scripts folder and run ``./setup.sh``. This creates other needed folders.
4. Use snakemake to generate final count tables and corresponding quality control graphics.
5. Use jupyter notebook to run analysis.

# Example

Tested with fastqc 0.11.5, cutadapt 1.18, R 4.1.0 (with reshape2, dplyr, ggrepel, ggplot2, gemini, configr), count_constructs 1.0, python 3.10.11 (with snakemake_7.25.3, pandas, numpy, click_8.1.3, statsmodels, toml, scipy)