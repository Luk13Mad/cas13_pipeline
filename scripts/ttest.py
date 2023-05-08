import click
import pandas as pd
import numpy as np
import pickle
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import norm,ttest_ind


@click.command()
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
@click.option("--comparison",required = True, help = "Comparison of timepoints.")
@click.option("--cellline",required = True, help = "Name of cellline.")
@click.option("--exclusionlist",type = click.Path(exists = True),help = "Path to fiel with guide combinations to exclude.")
def main(cellline : str,comparison : str,workingdir : str, exclusionlist : str = None) -> None:
	dicts_file = "all_" + cellline + "_" + comparison
	
	if exclusionlist is not None:
		ex_list = pd.read_csv(exclusionlist,sep="\t",index_col=0)

	with open(workingdir + '/dLFC_analysis/all_{0}_{1}/expected_vs_observed_dict_{0}_{1}.pickle'.format(cellline,comparison), 'rb') as handle:
		expected_vs_observed = pickle.load(handle)


	df = pd.DataFrame.from_dict(expected_vs_observed,orient="index")
	anno = pd.read_csv(workingdir + "metainfo/guide_anno.csv",sep="\t",index_col=0)
	df["dLFCZ_sum"] = df["observed"] - df["exp_sum"]
	df["guide1"] = [x.split(";")[0] for x in df.index]
	df["guide2"] = [x.split(";")[1] for x in df.index]
	df = df.merge(anno,left_index=True,right_index=True,how="left")
	df["sortby"] = df.loc[:,["spot1","spot2"]].apply(lambda x: ";".join(sorted(x)),axis=1)
	df["sortby_guides"] = df.loc[:,["guide1","guide2"]].apply(lambda x: ";".join(sorted(x)),axis=1)
	if exclusionlist is not None:
		df = df.loc[~df.index.isin(ex_list.index),:]

	res_dict = dict()

	for i,row in df.iterrows():

		if row.sortby in res_dict:
			continue
		if "control" in row.spot1 or "control" in row.spot2:
			continue

		mask_ctrl1 = (df.spot1.str.contains("control")) & (df.spot2.isin([row.spot1,row.spot2]))
		mask_ctrl2 = (df.spot1.isin([row.spot1,row.spot2])) & (df.spot2.str.contains("control"))

		mask_ctrl = mask_ctrl1 ^ mask_ctrl2 
		ctrl_all = df.loc[mask_ctrl,"dLFCZ_sum"].values

		mask_val = df.sortby.str.contains(row.spot1) & df.sortby.str.contains(row.spot2)
		val_all = df.loc[mask_val,"dLFCZ_sum"].values

		#pvalue
		if len(ctrl_all) > 2 and len(val_all) > 2:
			pval = ttest_ind(ctrl_all,val_all,axis=None,random_state = 13,equal_var = False).pvalue
			res_dict.update({row.sortby : {"mean_dLFC":np.mean(val_all),"std_dLFC":np.std(val_all),"pvalue":pval}})


	res_df = pd.DataFrame.from_dict(res_dict,orient = "index")
	mask = res_df.pvalue.isnull()
	res_df = res_df.loc[~mask,:]
	res_df["FDR"] = fdrcorrection(res_df.pvalue.values)[1]


	res_df.to_csv(workingdir + '/dLFC_analysis/all_{0}_{1}/ttest_{0}_{1}.tsv'.format(cellline,comparison),sep="\t")

    
if __name__=="__main__":
    main()
