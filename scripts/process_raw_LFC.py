import pandas as pd
import numpy as np
import sys,os
import math
import pickle
import click



def get_single_effect(iguide,idf,icontrol,level):
	'''
	iguide : current guide name
	idf : dataframe containing all needed LFC data
	icontrol : List with control guides
	'''
	if level == "guide":
		mask = (idf.gene1.isin(icontrol) & idf.guide2.isin([iguide])) | (idf.gene2.isin(icontrol) & idf.guide1.isin([iguide]))
		return idf.loc[mask,"LFCZ"].mean()
	elif level == "gene":
		mask = (idf.gene1.isin(icontrol) & idf.gene2.isin([iguide])) | (idf.gene2.isin(icontrol) & idf.gene1.isin([iguide]))
		return idf.loc[mask,"LFCZ"].mean()
    
def get_expected_and_observed(iguidepair,idf,icontrol,level):
	'''
	iguidepair : guidepair for which the expected vs the observed effect are to be measured
	method : method to use to generate the expected effect
	idf : dataframe containing all needed LFC data
	'''
	if level == "guide":
		partner1,partner2 = iguidepair.split(";")
		mask = (idf.guide1.isin([partner1]) & idf.guide2.isin([partner2])) | (idf.guide1.isin([partner2]) & idf.guide2.isin([partner1]) )
		observed = idf.loc[mask,"LFCZ"].mean()

		expected_sum,expected_product,expected_min,expected_max,expected_log = calc_expected(partner1,partner2,idf,icontrol,level)
		return iguidepair,{"observed":observed,"exp_sum":expected_sum,"exp_product":expected_product,"exp_min":expected_min,"exp_max":expected_max,"exp_log":expected_log}
	elif level == "gene":
		partner1,partner2 = iguidepair.split(";")
		mask = (idf.gene1.isin([partner1]) & idf.gene2.isin([partner2])) | (idf.gene1.isin([partner2]) & idf.gene2.isin([partner1]) )
		observed = idf.loc[mask,"LFCZ"].mean()
		expected_sum,expected_product,expected_min,expected_max,expected_log = calc_expected(partner1,partner2,idf,icontrol,level)
		return iguidepair,{"observed":observed,"exp_sum":expected_sum,"exp_product":expected_product,"exp_min":expected_min,"exp_max":expected_max,"exp_log":expected_log}




def calc_expected(partner1,partner2,idf,icontrol,level):
    single_fx_partner1 = get_single_effect(partner1,idf,icontrol,level)
    single_fx_partner2 = get_single_effect(partner2,idf,icontrol,level)

    expected_sum = single_fx_partner1 + single_fx_partner2
    expected_product = single_fx_partner1 * single_fx_partner2
    expected_min = min(single_fx_partner1,single_fx_partner2)
    expected_max = max(single_fx_partner1,single_fx_partner2)
    expected_log = math.log2((2**single_fx_partner1+1)*(2**single_fx_partner2+1)+1)

    return expected_sum,expected_product,expected_min,expected_max,expected_log


@click.command()
@click.option("--output_dir",type = click.Path(exists = True),required = True,help = "Path where results will be written to.")
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
@click.option("--cellline",required = True,help = "Name of cell line.")
@click.option("--comparison",required = True, help = "Specify comparison of timepoints, e.g. 13 for TP1vsTP3.")
@click.option("--lfc_file",type=click.Path(exists=True),required = True, help = "Path to file containing LFC values.")
def main(output_dir,workingdir,cellline,lfc_file,comparison):
	os.chdir(workingdir)
	dicts_file = "all_" + cellline + "_" + comparison + "/"

	df = pd.read_csv(lfc_file,sep="\t",low_memory=False,index_col=0)
	df["sortby"] = df.loc[:,["gene1","gene2"]].apply(lambda x: ";".join(sorted(x)),axis=1)


	#all control 
	all_controls = np.unique(np.concatenate((df.gene1.loc[df.gene1.str.contains("control")].unique(),df.gene2.loc[df.gene2.str.contains("control")].unique())))

	out_dict = {}

	#guidepair is combination of sgRNA
	for guidepair in df.index:
		pair,observed_expected_dict = get_expected_and_observed(guidepair,df,all_controls,"guide")
		if pair in out_dict:
			print("#Was already in ")
		out_dict.update({pair:observed_expected_dict})
	
	if not os.path.exists(output_dir + dicts_file):
		os.mkdir(output_dir + dicts_file)

	with open(output_dir + dicts_file + 'expected_vs_observed_dict_{}_{}.pickle'.format(cellline,comparison), 'wb') as handle:
		pickle.dump(out_dict, handle)


	out_dict = {}
	#do same on gene level
	for guidepair in df.sortby.unique():
		pair,observed_expected_dict = get_expected_and_observed(guidepair,df,all_controls,"gene")
		if pair in out_dict:
			print("#Was already in ")
		out_dict.update({pair:observed_expected_dict})

	with open(output_dir + dicts_file + 'expected_vs_observed_dict_{}_{}_genelevel.pickle'.format(cellline,comparison), 'wb') as handle:
		pickle.dump(out_dict, handle)

if __name__ == "__main__":
	main()
