import click
import pandas as pd
import numpy as np
import pickle
from statsmodels.stats.multitest import fdrcorrection

def geo_mean(x : np.ndarray) -> np.ndarray:
    '''calc geometric mean for array x'''
    return np.exp(np.log(x).mean())

def get_unique_sgrna(gene_name : str,df : pd.DataFrame) -> np.ndarray:
    '''Return unique sgrnas associated with gene_name'''
    mask1 = df.spot1 == gene_name
    mask2 = df.spot2 == gene_name
    return np.unique(np.concatenate((df.loc[mask1,"guide1"].values, df.loc[mask2,"guide2"].values), axis=None))

def get_ordered_genesymbols(sgrna : str,gene_name : str,sortby : str,df : pd.DataFrame) -> list[str]:
    '''Returns gene symbols ordered by dLFC from lowest to highest associated with sgrna'''
    mask = (df.guide1.isin([sgrna]) | df.guide2.isin([sgrna]))
    df = df.loc[mask,:]
    out = []
    for i,row in df.sort_values(sortby,ascending=True).iterrows():
        if row["spot1"] == gene_name:
            out.append(row["spot2"])
        else:
            out.append(row["spot1"])
    return out

def normalize_ranks(ranks : list[int],max_rank : int,min_rank : int =1) -> list[float]:
    '''normalize the ranks to the interval [0,1]'''
    if len(ranks)!=0:
        out = (np.array(ranks) - min_rank) / (max_rank+min_rank)
        return list(out)
    else:
        return []

def get_averaged_rank(ordered_dict,current_test_gene,all_test_genes) -> dict:
    '''Returns average rank of partner gene with test gene'''
    
    #all genes from all_test_genes except the current gene
    all_test_genes_small = [x for x in all_test_genes if x!=current_test_gene]
    
    out_dict = dict()
    for gene in all_test_genes_small:
        out_value = []
        for sgrna in ordered_dict.keys():
            max_rank = len(ordered_dict[sgrna])
            out_value = out_value + normalize_ranks([i+1 for i,e in enumerate(ordered_dict[sgrna]) if e == gene],max_rank)
        out_key = current_test_gene + ";" + gene
        if len(out_value)<3:#at least 3 ranks have to be there to compute average
            out_dict.update({out_key:np.nan})
        else:
            out_dict.update({out_key:(np.mean(out_value),np.std(out_value))})
    return out_dict


def calc_pvalue(test_rank,null_dist,null_dist_amount,kind = "leftside"):
    '''Calcs pvalue based on permutation simulations. Vectorized using numpy'''
    #simulations_as_or_more_extreme = np.sum(np.array(null_dist) >= test_rank)
    if kind == "leftside":
        simulations_as_or_more_extreme = np.sum(np.less_equal(np.array(null_dist), test_rank[:,np.newaxis]),axis=1)
        return simulations_as_or_more_extreme/null_dist_amount
    elif kind == "rightside":
        simulations_as_or_more_extreme = np.sum(np.greater_equal(np.array(null_dist), test_rank[:,np.newaxis]),axis=1)
        return simulations_as_or_more_extreme/null_dist_amount
    else:
        raise ValueError

def get_shuffled_df(original_df,new_index):
    '''Shuffles label of df to be used to generate null distribution'''
    rank_and_old_index = original_df.loc[:,["dLFCZ_sum"]].copy()
    rank_and_old_index.reset_index(drop=True,inplace=True)
    
    shuffled_spot1_sgrna1 = original_df.loc[new_index,["spot1","guide1"]].copy()
    shuffled_spot1_sgrna1.reset_index(drop=True,inplace=True)
    
    shuffled_spot2_sgrna2 = original_df.loc[new_index,["spot2","guide2"]].copy()
    shuffled_spot2_sgrna2.reset_index(drop=True,inplace=True)
    return pd.concat([rank_and_old_index,shuffled_spot1_sgrna1,shuffled_spot2_sgrna2],axis=1)


def run_dLFCZ_rank_generation(test_genes_all : np.ndarray,test_df : pd.DataFrame) -> pd.DataFrame:
    all_results = []

    for t_gene in test_genes_all:
        test_sgrnas = get_unique_sgrna(t_gene,test_df)
        test_sgrnas_dict = dict()
        for t_sgrnas in test_sgrnas:
            ordered_symbols = get_ordered_genesymbols(t_sgrnas,t_gene,"dLFCZ_sum",test_df)
            test_sgrnas_dict.update({t_sgrnas:ordered_symbols})
        result = get_averaged_rank(test_sgrnas_dict,t_gene,test_genes_all)
        all_results.append(result)
    
    #flatten list of dicts to one dict
    all_results_dict = {k:v for d in all_results for k,v in d.items()}
    
    final_res = pd.DataFrame(all_results_dict).T
    final_res.columns = ["rank","std"]
    final_res = final_res.loc[~final_res.loc[:,"rank"].isna(),:]
    return final_res

def produce_nulldist(final_res : pd.DataFrame,test_df : pd.DataFrame,test_genes_all : np.ndarray,samplesize : int) -> pd.DataFrame:
    '''generate null distribution by switching labels randomly'''
    
    null_dist_df = pd.DataFrame(index=final_res.index)
    for counter in range(samplesize):
        new_idx = np.random.permutation(test_df.index)
        shuffled_df = get_shuffled_df(test_df,new_idx)
        tmp = run_dLFCZ_rank_generation(test_genes_all,shuffled_df)
        null_dist_df = null_dist_df.merge(tmp,left_index=True,
                                          right_index=True,
                                          how="left",suffixes=(None,"_"+str(counter)))
    
    mylen = null_dist_df.melt().dropna().shape[0]
    final_res["pvalue_leftside"] = calc_pvalue(final_res.loc[:,"rank"].values,
                              null_dist_df.melt().dropna().loc[:,"value"],
                              mylen,
                              "leftside")
    final_res["pvalue_rightside"] = calc_pvalue(final_res.loc[:,"rank"].values,
                                           null_dist_df.melt().dropna().loc[:,"value"],
                                           mylen,
                                           "rightside")
    
    #remove significnat gene-combinations from first run
    #repeate sampling to form null distribution
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5518132/
    mask = (final_res["pvalue_leftside"] <= 0.05) | (final_res["pvalue_rightside"] <= 0.05)
    mask = final_res.index[mask]
    mask2 = test_df.sortby1.isin(mask) | test_df.sortby2.isin(mask) #filter out gene pairs which were significant in first round
    
    null_dist_df = pd.DataFrame(index=final_res.index)
    for counter in range(samplesize):
        new_idx = np.random.permutation(test_df.loc[~mask2,:].index)
        shuffled_df = get_shuffled_df(test_df.loc[~mask2,:],new_idx)
        tmp = run_dLFCZ_rank_generation(test_genes_all,shuffled_df)
        null_dist_df = null_dist_df.merge(tmp,left_index=True,
                                          right_index=True,
                                          how="left",
                                          suffixes=(None,"_"+str(counter)))
    null_dist_df = null_dist_df.loc[:,[x for x in null_dist_df.columns if "rank" in x]]
    return null_dist_df

@click.command()
@click.option("--samplesize",type = click.INT, default = 100,help = "Sample size for permutation testing.")
@click.option("--workingdir",type = click.Path(exists = True),required = True,help = "Path to directory of screen.")
@click.option("--comparison",required = True, help = "Comparison of timepoints.")
@click.option("--gene_exclusion", default=None,help = "List of genes to exclude.")
@click.option("--cellline",required = True, help = "Name of cellline.")
def main(cellline : str,comparison : str,workingdir : str,gene_exclusion : list[str] = None,samplesize : int = 100) -> None:
    dicts_file = "all_" + cellline + "_" + comparison
    
    with open(workingdir + '/dLFC_analysis/all_{0}_{1}/expected_vs_observed_dict_{0}_{1}.pickle'.format(cellline,comparison), 'rb') as handle:
        expected_vs_observed = pickle.load(handle)
    

    df = pd.DataFrame.from_dict(expected_vs_observed,orient="index")
    anno = pd.read_csv(workingdir + "metainfo/guide_anno.csv",sep="\t",index_col=0)
    df["dLFCZ_sum"] = df["observed"] - df["exp_sum"]
    df["guide1"] = [x.split(";")[0] for x in df.index]
    df["guide2"] = [x.split(";")[1] for x in df.index]
    df = df.merge(anno,left_index=True,right_index=True,how="left")
    print("#Loaded data.")

    #remove all rows with at least one control
    mask = df.spot1.str.contains("control") | df.spot2.str.contains("control")
    #test_df has no control guides in it
    test_df = df.loc[~mask,:].copy()
    test_df["sortby1"] = test_df["spot1"] + ";" +test_df["spot2"]
    test_df["sortby2"] = test_df["spot2"] + ";" +test_df["spot1"]
    
    if gene_exclusion is not None:
        print("#Removing specified genes.")
        mask = test_df.spot1.isin(gene_exclusion) | test_df.spot2.isin(gene_exclusion)
        test_df = test_df.loc[~mask,:]
    
    #array of all unique genes
    test_genes_all = np.unique(test_df.loc[:,["spot1","spot2"]].values.flatten())
    
    print("#Starting first rank test.")
    final_res = run_dLFCZ_rank_generation(test_genes_all,test_df)
    
    print("#Started producing null dist.")
    null_dist_df = produce_nulldist(final_res,test_df,test_genes_all,samplesize=samplesize)
    
    print("#Savin null distribution to file")
    null_dist_df.to_csv(workingdir + '/dLFC_analysis/all_{0}_{1}/nulldist_{0}_{1}.tsv'.format(cellline,comparison),sep="\t")
    
    mylen = null_dist_df.melt().dropna().shape[0]
    final_res["pvalue_leftside"] = calc_pvalue(final_res.loc[:,"rank"].values,
                              null_dist_df.melt().dropna().loc[:,"value"],
                              mylen,
                              "leftside")
    final_res["pvalue_rightside"] = calc_pvalue(final_res.loc[:,"rank"].values,
                                           null_dist_df.melt().dropna().loc[:,"value"],
                                           mylen,
                                           "rightside")
    final_res["FDR_leftside"] = fdrcorrection(final_res["pvalue_leftside"].values)[1]
    final_res["FDR_rightside"] = fdrcorrection(final_res["pvalue_rightside"].values)[1]
    
    final_res.to_csv(workingdir + '/dLFC_analysis/all_{0}_{1}/ranktest_{0}_{1}.tsv'.format(cellline,comparison),sep="\t")

    
if __name__=="__main__":
    main()
