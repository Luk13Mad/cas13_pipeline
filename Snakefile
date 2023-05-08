import glob


def get_input_raw_data_fastqc(wildcards):
	files=[x for x in os.listdir("raw_sequencing_data/{}/raw_data/{}/".format(wildcards.run,wildcards.type)) if x.endswith(".fq.gz")]
	return [os.path.join("raw_sequencing_data/{}/raw_data/{}/".format(wildcards.run,wildcards.type),x) for x in files]	

rule fastqc_raw_data:
	input:
		fastqfile=get_input_raw_data_fastqc
	output:
		#directory("raw_sequencing_data/FastQC/{run}/{type}")
		"raw_sequencing_data/FastQC/{run}/{type}/FastQC_raw_data.done"
	shell:
		'''
		module load FastQC/0.11.5-Java-1.8.0_112
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
		#directory("processed_sequencing_data/{run,[^/]*}/{type,[^/]*}")
		"processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
		"processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz"
	shell:
		'''
		module load cutadapt/1.9.1-foss-2017a
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}
		cutadapt --minimum-length=100 -a CTGTCTCTTATA -A CTGTCTCTTATA -o /b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/{wildcards.run}/{wildcards.type}/{wildcards.type}_L1_1_trimmed.fq.gz -p /b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/{wildcards.run}/{wildcards.type}/{wildcards.type}_L1_2_trimmed.fq.gz /b06x-isi/b062/a-c/Braun_CRISPR_2D/{input.forward_reads} /b06x-isi/b062/a-c/Braun_CRISPR_2D/{input.reverse_reads}
		'''


rule fastqc_trimmed_data:
	input:
		forward_reads="processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
                reverse_reads="processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz"
	output:
		#directory("processed_sequencing_data/{run}/{type}/FastQC")
		"processed_sequencing_data/{run}/{type}/FastQC/trimmed_data.done"
	shell:
		'''
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC
		touch processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC/trimmed_data.done
		module load FastQC/0.11.5-Java-1.8.0_112
		fastqc --noextract -o /b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/{wildcards.run}/{wildcards.type}/FastQC -t 4 /b06x-isi/b062/a-c/Braun_CRISPR_2D/{input.forward_reads} /b06x-isi/b062/a-c/Braun_CRISPR_2D/{input.reverse_reads} 
		'''

rule count_constructs:
	input:
		"processed_sequencing_data/{run}/{type}/FastQC/trimmed_data.done",
		forward_reads="processed_sequencing_data/{run}/{type}/{type}_L1_1_trimmed.fq.gz",
		reverse_reads="processed_sequencing_data/{run}/{type}/{type}_L1_2_trimmed.fq.gz",
		expected="processed_sequencing_data/{run}/{type}/{type}_constructs.csv"
	output:
		"processed_sequencing_data/{run}/{type}/count_results_{type}.tsv"
	shell:
		'''
		module load anaconda3/2021.05
		source activate cas13_pipeline 
		python /b06x-isi/b062/a-c/Braun_CRISPR_2D/scripts/count_constructs/main.py 2D --resfile count_results_{wildcards.type}.tsv --samplefolder  /b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/{wildcards.run}/{wildcards.type}/ -f {wildcards.type}_L1_1_trimmed.fq.gz -r {wildcards.type}_L1_2_trimmed.fq.gz -e {wildcards.type}_constructs.csv -s1 58 -e1 81 -s2 50 -e2 90 -w 20 --batchsize 100000
		'''




rule make_summary_graphics:
	input:
		"processed_sequencing_data/{run}/{type}/count_results_{type}.tsv"
	output:
		"processed_sequencing_data/{run}/{type}/graphics/summary_graphics_{type}.pdf"
	shell:
		'''
		module load R-bundle/20210323-foss-2020a-R-4.0.4
		mkdir -p processed_sequencing_data/{wildcards.run}/{wildcards.type}/graphics
		sleep 20
		Rscript /b06x-isi/b062/a-c/Braun_CRISPR_2D/scripts/summary_graphics.R {input} {wildcards.run} {wildcards.type}
		'''
