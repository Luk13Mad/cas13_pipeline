import click
import os
import pandas as pd
import glob

#sum over different sort of counts and remove not needed columns



def clean_one_file(filename,outputdir,normalize):
	df = pd.read_csv(filename,sep="\t")
	filename = filename.split(".")[0]
	filename = filename.split("/")[-1]
	filename = "_".join(filename.split("_")[2:])
	df["sum"] = df.loc[:,["count_0MM_0MM","count_0MM_1MM","count_1MM_0MM","count_1MM_1MM"]].sum(axis=1)
	df = df.drop(columns=["gene1", "gene2", "gene1_sgrnaid", "gene2_sgrnaid","count_0MM_0MM","count_0MM_1MM","count_1MM_0MM","count_1MM_1MM"])
	df.index = df[['seq1', 'seq2']].agg(';'.join, axis=1)
	df = df.drop(columns=["seq1","seq2"])
	df.columns=[filename]
	read_sum = df[filename].sum()
	if normalize == True: #normalize to counts per million
		df[filename] = df[filename].map(lambda x: (x/read_sum)*1000000)
		df.to_csv(outputdir+filename+"_normalized_cleaned.csv",sep="\t")
	else:
		df.to_csv(outputdir+filename+"_raw_cleaned.csv",sep="\t")

def collect_all_files(directory,outputdir,normalize):
	count_files = glob.glob("processed_sequencing_data/"+directory+"/*/count_results_*")
	for f in count_files:
		clean_one_file(f,outputdir,normalize)
	
@click.command()
@click.option("--normalize",is_flag = True, help = "Normalize as CPM.")
@click.option("--directory",required = True,help = "Subfolder with count results.")
@click.option("--outputdir",type = click.Path(exists = True),required = True,help = "Path where results will be written to.")
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
def main(workingdir,normalize,outputdir,directory):
	os.chdir(workingdir)
	collect_all_files(directory,outputdir,normalize)

if __name__ == "__main__":
	main()
