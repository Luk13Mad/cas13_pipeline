import pandas as pd
import numpy as np
import os
import sys
import click
import toml

def MAD(data: np.ndarray) -> float:
    """
    Median absolute deviation.
    MAD = median(|X_i - median(X)|)
    """
    return np.median(np.abs(data - np.median(data)))

def print_LFC_to_file(LFC_col,df,output_dir):
	'''
	print specified LFC to file
	'''
	_,cellline,_comparison = LFC_col.split("_")
	filename = "raw_" + LFC_col.replace("LFC","LFCZ") + ".tsv"
	out = df.loc[:,[LFC_col,'spot1', 'spot2','sgrna1', 'sgrna2']]
	out.columns = ["LFCZ","gene1","gene2","guide1","guide2"]
	out.to_csv(output_dir + filename,sep="\t")

def calc_LFC(comparison:list,cellline,df):
    '''
    calc LFC between two timepoints for specified cellline
    '''
    ETP,LTP = comparison
    colname_ETP = cellline + "TP" + str(ETP)
    colname_LTP = cellline + "TP" + str(LTP)
    return (df.loc[:,colname_LTP] - df.loc[:,colname_ETP])

def load_data(path_count_file,path_anno_file,path_exclusion_file,path_sample_anno):
	counts_df = pd.read_csv(path_count_file,sep="\t",index_col=0)
	TTTT_gate = pd.read_csv(path_exclusion_file,sep="\t",index_col=0)
	sample_anno_df = pd.read_csv(path_sample_anno,sep="\t")
	guide_anno = pd.read_csv(path_anno_file,sep="\t",index_col=0)
	guide_anno["merged_genes"] = guide_anno["spot1"] + ";" + guide_anno["spot2"]
	guide_anno["sgrna1"] = [x[0] for x in guide_anno.index.str.split(";")]
	guide_anno["sgrna2"] = [x[1] for x in guide_anno.index.str.split(";")]
	
	counts_df = counts_df.loc[~counts_df.index.isin(TTTT_gate["combined_guides"]),:] #remove TTTT_gate

	return counts_df,guide_anno,sample_anno_df

def data_massaging(count_df,sample_anno_df,cellline,earlyTP):
	'''calcs median over replicates and removes rows with too little counts in earlyTP'''	
	sample_anno_df = sample_anno_df.loc[sample_anno_df.samplename.str.contains(cellline),:]

	#calc mean over replicates
	for tpname,subdf in sample_anno_df.groupby("samplename"):
		count_df[tpname] = count_df.loc[:,subdf["colname"].tolist()].mean(axis=1)


	return count_df



@click.command()
@click.option("--output_dir",type = click.Path(exists = True),required = True,help = "Path where results will be written to.")
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
@click.option("--path_count_file",type = click.Path(exists=True),required = True, help = "Path to file with counts.")
@click.option("--path_config_file",type = click.Path(exists=True),required = True, help = "Path to file with config.")
@click.option("--path_anno_file",type = click.Path(exists=True),required = True, help = "Path to file with guide annotation.")
@click.option("--path_exclusion_file",type = click.Path(exists=True),required = True, help = "Path to file with exclusion file.")
@click.option("--path_sample_anno",type = click.Path(exists=True),required = True, help = "Path to file with sample annotation.")
@click.option("--constant",type = click.INT,default = 32, help = "Pseudocount to add before calculating LFC.")
def main(output_dir,workingdir,path_count_file,path_config_file,path_anno_file,path_exclusion_file,path_sample_anno,constant):
	os.chdir(workingdir)
	config = toml.load(path_config_file)
	counts_df,guide_anno_df,sample_anno_df = load_data(path_count_file,path_anno_file,path_exclusion_file,path_sample_anno)

	#data massaging, averaging over replicates
	counts_df = data_massaging(counts_df,sample_anno_df,config["cellline"],config["prepareLFC"]["earlyTP"])

	#guide_anno_df cointains same rows as counts_df
	guide_anno_df = guide_anno_df.loc[counts_df.index,:]

	#transforming counts wiht log and pseudocount 
	counts_df = counts_df.apply(lambda x: np.log2(x + constant),axis=0)

	#merge annotation with count data
	counts_df = counts_df.merge(guide_anno_df,left_index=True,right_index=True)

	#comparisons is list of lists e.g. [[1,2],[1,3]]	
	for comp in config["prepareLFC"]["comparisons"]:
		comp_str = "".join([str(x) for x in comp])
		counts_df["LFC_{}_{}".format(config["cellline"],comp_str)] = calc_LFC(comp,"{}".format(config["cellline"]),counts_df)

	LFC_columns = [x for x in counts_df.columns if "LFC" in x]

	#replace INF with 0 in LFC columns
	counts_df.loc[:,LFC_columns] = counts_df.loc[:,LFC_columns].replace([np.inf, -np.inf], 0)

	#mask with all NT_NT pairs
	NTNTmask = (counts_df.spot1.str.contains("control") & counts_df.spot2.str.contains("control"))
	NTNTmean = counts_df.loc[NTNTmask,LFC_columns].mean(axis=0)

	#subtract NTNT mean from all guide pairs
	counts_df.loc[:,LFC_columns] = counts_df.loc[:,LFC_columns] - NTNTmean


	for n in [x for x in counts_df.columns if "LFC" in x]:
		print_LFC_to_file(n,counts_df,output_dir)

if __name__ == "__main__":
        main()

