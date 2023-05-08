import pandas as pd
import os
import click
import glob
from functools import reduce


def collect_all_files(counts_dir,normalize):
	if normalize:
		return glob.glob(counts_dir + "*normalized_cleaned.csv")
	else:
		return glob.glob(counts_dir + "*raw_cleaned.csv")


def merge_all_files(filelist,cellline,counts_dir):
	dataframes = [pd.read_csv(f,sep="\t",index_col=0) for f in filelist]
	dfmerged = reduce(lambda  left,right: pd.merge(left,right,left_index=True,right_index=True), dataframes)
	
	if all(["normalized" in x for x in filelist]) == True:
		status = "norm"
	else:
		status = "raw"
	
	dfmerged.to_csv(counts_dir+"merged_counts_{}_{}.csv".format(cellline,status),sep="\t")


def remove_files(filelist):
	for f in filelist:
		try:
			os.remove(f)
		except Exception as e:
			print(e)

@click.command()
@click.option("--cellline",required = True,help = "Name of cell line, will be used in name of merged file.")
@click.option("--counts_dir",type = click.Path(exists = True),required = True,help = "Path where count results are found.")
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
@click.option("--normalize",is_flag = True,required = True,help = "Flag if merge should be run on  normalized or raw files.")
def main(cellline,counts_dir,workingdir,normalize):
	os.chdir(workingdir)
	files_to_merge = collect_all_files(counts_dir,normalize)
	merge_all_files(files_to_merge,cellline,counts_dir)
	remove_files(files_to_merge)

if __name__ == "__main__":
        main()

