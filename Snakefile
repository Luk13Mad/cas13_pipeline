import glob
import os

WORKING_DIR = workflow.basedir

def get_input_raw_data_fastqc(wildcards):
	files=[x for x in os.listdir("raw_sequencing_data/{}/raw_data/{}/".format(wildcards.run,wildcards.type)) if x.endswith(".fq.gz")]
	return [os.path.join("raw_sequencing_data/{}/raw_data/{}/".format(wildcards.run,wildcards.type),x) for x in files]	

rule fastqc_raw_data:
	input:
		fastqfile=get_input_raw_data_fastqc
	output:
		"raw_sequencing_data/FastQC/{run}/{type}/FastQC_raw_data.done"
	shell:
		'''
		mkdir -p raw_sequencing_data/FastQC/{wildcards.run}/{wildcards.type}
		touch raw_sequencing_data/FastQC/{wildcards.run}/{wildcards.type}/FastQC_raw_data.done
		fastqc --noextract -o raw_sequencing_data/FastQC/{wildcards.run}/{wildcards.type} -t 4 {input.fastqfile}
		'''

def get_fwd_cutadapt(wildcards):
	return glob.glob("raw_sequencing_data/{}/raw_data/{}/*_1.fq.gz".format(wildcards.run,wildcards.type))

def get_rev_cutadapt(wildcards):
	return glob.glob("raw_sequencing_data/{}/raw_data/{}/*_2.fq.gz".format(wildcards.run,wildcards.type))

rule cutadapt_raw_data:
	input:
		"raw_sequencing_data/FastQC/{run}/{type}/FastQC_raw_data.done",
		forward_reads=get_fwd_cutadapt,
		reverse_reads=get_rev_cutadapt
	output:
		"processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
		"processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz"
	params:
		wd = WORKING_DIR
	shell:
		'''
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}
		cutadapt --minimum-length=100 -a CTGTCTCTTATA -A CTGTCTCTTATA -o {params.wd}/processed_sequencing_data/{wildcards.run}/{wildcards.type}/{wildcards.type}_L1_1_trimmed.fq.gz -p {params.wd}/processed_sequencing_data/{wildcards.run}/{wildcards.type}/{wildcards.type}_L1_2_trimmed.fq.gz {params.wd}/{input.forward_reads} {params.wd}/{input.reverse_reads}
		'''


rule fastqc_trimmed_data:
	input:
		forward_reads="processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
                reverse_reads="processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz"
	output:
		"processed_sequencing_data/{run}/{type}/FastQC/trimmed_data.done"
	params:
		wd = WORKING_DIR
	shell:
		'''
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC
		touch processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC/trimmed_data.done
		fastqc --noextract -o {params.wd}/processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC -t 4 {params.wd}/{input.forward_reads} {params.wd}/{input.reverse_reads}
		'''

rule count_constructs:
	input:
		"processed_sequencing_data/{run}/{type}/FastQC/trimmed_data.done",
		forward_reads="processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
		reverse_reads="processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz",
	output:
		"processed_sequencing_data/{run}/{type}/count_results_{type}.tsv"
	params:
		wd = WORKING_DIR
	shell:
		'''
		count_constructs 2D --resfile count_results_{wildcards.type}.tsv --samplefolder {params.wd}/processed_sequencing_data/{wildcards.run}/{wildcards.type}/ -f {wildcards.type}_L1_1_trimmed.fq.gz -r {wildcards.type}_L1_2_trimmed.fq.gz -e {wildcards.type}_constructs.csv -s1 58 -e1 81 -s2 50 -e2 90 -w 20 --batchsize 100000
		'''




rule make_summary_graphics:
	input:
		"processed_sequencing_data/{run}/{type}/count_results_{type}.tsv"
	output:
		"processed_sequencing_data/{run}/{type}/graphics/summary_graphics_{type}.pdf"
	shell:
		'''
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}/graphics
		sleep 20
		Rscript scripts/summary_graphics.R {input} {wildcards.run} {wildcards.type}
		'''
