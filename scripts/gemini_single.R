library(gemini)
library(configr)

better_merge = function(x,y){
  m = merge(x,y, all=TRUE,by = 0)
  rownames(m) = m$Row.names
  m = m[,names(m)!="Row.names"]
  return(m)
}


#arguments from command line
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
	stop("Exactly 7 arguments must be specified.", call.=FALSE)
}else{
	guide_annotation_file = args[1] #path to guide annotation file
	sample_annotation_file = args[2] #path to sample annotation file
	counts_file = args[3] #path to file with counts
	exclusion_file = args[4] #path to file with guide combnations to exclude
	t = args[5] #number of cores to use
	current_wd = args[6] #current working dir , the directory where the merged count files is located
	config_file_path = args[7]
}

#guide_annotation_file = "/b06x-isi/b062/a-c/Braun_CRISPR_2D/metainfo/guide_anno.csv"
#sample_annotation_file = "/b06x-isi/b062/a-c/Braun_CRISPR_2D/metainfo/sample_replicate_anno.csv"
#counts_file = "/b06x-isi/b062/a-c/Braun_CRISPR_2D/GEMINI_input/X204SC20113232-Z01-F010/merged_counts_PATU_raw.csv"
#exclusion_file = "/b06x-isi/b062/a-c/Braun_CRISPR_2D/metainfo/combinations_affected_by_TTTTgate.tsv"
#t = 1
#current_wd =  "/b06x-isi/b062/a-c/Braun_CRISPR_2D/"
#config_file_path = "/b06x-isi/b062/a-c/Braun_CRISPR_2D/GEMINI_input/X204SC20113232-Z01-F010/config_PATU.toml"

  
  
config = read.config(config_file_path) #red config.toml file

#pass values from config file to variables
ETP_col_index = config$GEMINI$ETP_col_index
LTP_col_index = config$GEMINI$LTP_col_index
pc_genes = config$GEMINI$pc_genes
sample_anno_rows = config$GEMINI$sample_anno_rows
ETP_col_index = config$GEMINI$ETP_col_index
LTP_col_index = config$GEMINI$LTP_col_index
mean_x = config$GEMINI$mean_x
sd_x = config$GEMINI$sd_x
mean_xx = config$GEMINI$mean_xx
sd_xx = config$GEMINI$sd_xx
mean_y = config$GEMINI$mean_y
sd_y = config$GEMINI$sd_y
mean_s = config$GEMINI$mean_s
sd_s = config$GEMINI$sd_s

setwd(current_wd)
	

#load annotation tables
guide.annotation = read.table(guide_annotation_file,sep="\t",header=T)
sample.replicate.annotation = read.table(sample_annotation_file,sep="\t",header=T)
sample.replicate.annotation = sample.replicate.annotation[sample_anno_rows,]


#load counts from celllines
counts = read.table(counts_file,sep="\t",header=T,row.names = 1)
counts = counts[,c(ETP_col_index,LTP_col_index)]

######
#TTTTgate exclusion or general exclusion of gene combination
TTTTgate = read.table(exclusion_file,header=T)
counts_mask = rownames(counts) %in% TTTTgate$combined_guides
counts = counts[!counts_mask,]
anno_mask =guide.annotation$X %in% TTTTgate$combined_guides
guide.annotation = guide.annotation[!anno_mask,]


Input <- gemini_create_input(counts.matrix = counts,
                             sample.replicate.annotation = sample.replicate.annotation,
                             guide.annotation = guide.annotation,
                             ETP.column = c(1,2),
                             LTP.column = c(3,4),
                             gene.column.names = c("spot1", "spot2"),
                             sample.column.name = "samplename",
                             verbose = TRUE)

Input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                CONSTANT = 32)



nc_genes_spot1=unique(grep("control_guide",guide.annotation$spot1,value=TRUE))
nc_genes_spot2=unique(grep("control_guide",guide.annotation$spot2,value=TRUE))

nc_genes_for_initialization = c(sample(nc_genes_spot1,(length(nc_genes_spot1)-1)),sample(nc_genes_spot2,(length(nc_genes_spot2)-1)))

Model = gemini_initialize(Input = Input, 
                           nc_gene = nc_genes_for_initialization, 
                           pattern_join = ';',
                           pattern_split = ';', 
                           cores = t,
                           verbose = TRUE)


Model %<>% gemini_inference(cores = t,
                            mean_x = mean_x,
                            sd_x = sd_x,
                            mean_xx = mean_xx,
                            sd_xx = sd_xx,
                            mean_y = mean_y,
                            sd_y = sd_y,
                            mean_s = mean_s,
                            sd_s = sd_s,
                            force_results = FALSE,
                            verbose = TRUE,
                            save_iterations = FALSE)

nc_pairs = grep("control_guide", rownames(Model$s), value = TRUE)
Score = gemini_score(Model = Model,
                      pc_gene = pc_genes,
                      nc_pairs = nc_pairs)

df_list = list()
df_list[[1]] = as.data.frame(Score[["strong"]][order(Score[["strong"]],decreasing = T),])
colnames(df_list[[1]]) = c("score_strong")
df_list[[2]] = as.data.frame(Score[["pvalue_strong"]][order(Score[["pvalue_strong"]],decreasing = F),])
colnames(df_list[[2]]) = c("pvalue_strong")
df_list[[3]] = as.data.frame(Score[["fdr_strong"]][order(Score[["fdr_strong"]],decreasing = F),])
colnames(df_list[[3]]) = c("FDR_strong")

df_list[[4]] = as.data.frame(Score[["sensitive_lethality"]][order(Score[["sensitive_lethality"]],decreasing = T),])
colnames(df_list[[4]]) = c("score_sensitive_lethality")
df_list[[5]] = as.data.frame(Score[["pvalue_sensitive_lethality"]][order(Score[["pvalue_sensitive_lethality"]],decreasing = F),])
colnames(df_list[[5]]) = c("pvalue_sensitive_lethality")
df_list[[6]] = as.data.frame(Score[["fdr_sensitive_lethality"]][order(Score[["fdr_sensitive_lethality"]],decreasing = F),])
colnames(df_list[[6]]) = c("FDR_sensitive_lethality")

df_list[[7]] = as.data.frame(Score[["sensitive_recovery"]][order(Score[["sensitive_recovery"]],decreasing = T),])
colnames(df_list[[7]]) = c("score_sensitive_recovery")
df_list[[8]] = as.data.frame(Score[["pvalue_sensitive_recovery"]][order(Score[["pvalue_sensitive_recovery"]],decreasing = F),])
colnames(df_list[[8]]) = c("pvalue_sensitive_recovery")
df_list[[9]] = as.data.frame(Score[["fdr_sensitive_recovery"]][order(Score[["fdr_sensitive_recovery"]],decreasing = F),])
colnames(df_list[[9]]) = c("FDR_sensitive_recovery")

out_table = Reduce(function(x, y) better_merge(x, y), df_list)
write.table(out_table,file.path(dirname(counts_file),"GEMINI_results.tsv"),quote=FALSE,sep="\t")



