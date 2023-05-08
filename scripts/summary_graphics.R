args = commandArgs(trailingOnly=TRUE)
count_matrix = args[1] #path to file with count_Results
run = args[2]
type = args[3]

#call with summary_graphics.R /b06x-isi/.../X204SC20113232-Z01-F007/count_results_A1_1_4.csv X204SC20113232-Z01-F007 A1_1_4

library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)

gini=function(data){
  data_sorted=sort(data)
  N=length(data)
  i=seq(1:N)
  upper=2*sum(data_sorted*i)
  lower=as.numeric(N)*sum(data_sorted)
  gini=(upper/lower)-((N-1)/N)
  return(gini)
}

lorenz=function(data){
  N=length(data)
  x=seq(1:N)/N
  data_sorted=sort(data)
  data_normalized=data_sorted/sum(data_sorted)
  y=cumsum(data_normalized)
  return(data.frame(cbind(x,y)))
}

hist_by_count=function(data,col,xlab,title){
  hist=ggplot(data,aes(count))+
    geom_histogram(binwidth = 10,color=col,fill=col)+
    xlab(xlab)+
    ggtitle(title)+
    geom_vline(xintercept=mean(data$count),col="red")+
    annotate("label",x=mean(data$count),y=50,label="mean",col="red")
  return(hist)
}

box_by_count=function(data,col,xlab,title){
  box=ggplot(data,aes(count))+
    geom_boxplot(color=col)+
    xlab(xlab)+
    ggtitle(title)
  return(box)
}

lorenz_by_count=function(data,col,xlab,ylab,title){
  lorenz=ggplot(lorenz(data$count),aes(x=x,y=y))+
    geom_point(color=col)+
    geom_abline(intercept = 0, slope = 1)+
    annotate("label", x = 0.3, y = 0.7, label = paste("GINI : ",round(gini(data$count),4)))+
    ggtitle(title)+
    xlab(xlab)+
    ylab(ylab)+
    scale_y_continuous(limits = c(0,1), expand = c(0, 0))+
    scale_x_continuous(limits = c(0,1), expand = c(0, 0))
  return(lorenz)
}

graphics=list()


df1=read.csv(count_matrix,sep="\t",comment.char = "#")  #count results as input
df1["count"]=rowSums(df1[,c(7,8,9,10)])

histdf1=df1 %>% group_by(gene1,gene2) %>% summarise(count=mean(count))
graphics[[1]]=hist_by_count(histdf1,"navy","counts per gene",paste("Histogram unnormalized counts, mean grouped by genes,",type))

graphics[[2]]=box_by_count(histdf1,"navy","counts per gene",paste("Boxplot unnormalized counts, mean grouped by genes,",type))

graphics[[3]]=lorenz_by_count(histdf1,"navy","Cumulative fraction of genes","Cumulative fraction of reads",paste("Lorenz curve unnormalized counts, mean grouped by genes,",type))

graphics[[4]]=hist_by_count(df1,"limegreen","counts per sgRNA",paste("Histogram unnormalized counts,",type))

graphics[[5]]=box_by_count(df1,"limegreen","counts per sgRNA",paste("Boxplot unnormalized counts,",type))

graphics[[6]]=lorenz_by_count(df1,"limegreen","Cumulative fraction of sgRNAs","Cumulative fraction of reads",paste("Lorenz curve unnormalized counts,",type))



#########################################################################
#########################################################################
#########################################################################

target_file = paste("/b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/",run,"/",type,"/",type,"_L1_1_trimmed.fq.gz",collapse = "",sep="")
total_lines_fq = as.integer(system2("zcat",args = c(target_file," | wc -l"),stdout = TRUE))
processed = (total_lines_fq/4)
total_success = sum(df1[,11])
discarded_reads = processed - total_success



total_count=data.frame(type=c(total_success,discarded_reads,processed),row.names = c("total_success","discarded","processed"))
colnames(total_count)=type
total_count=total_count/total_count[3,type]
total_count_df=total_count[1:2,,drop=F]
total_count_df=melt(as.matrix(total_count_df))



split_success_df=data.frame(colSums(df1[,c(7,8,9,10)]))
colnames(split_success_df)=type
split_success_df=split_success_df/sum(df1[,11])
split_success_df=melt(as.matrix(split_success_df))


graphics[[16]]= ggplot(total_count_df, aes(x = "", y=value, fill = Var1)) + 
  geom_col(width = 1, position="fill") +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5))+ 
  labs(fill="type", 
       x=NULL, 
       y=NULL, 
       title="Pie Chart total counts")+
  facet_grid(~Var2)+
  geom_label(aes(label=round(value,2)),nudge_x = 1,show.legend = F)+
  coord_polar(theta = "y", start=0)

graphics[[17]] = ggplot(split_success_df, aes(x = "", y=value, fill = Var1)) + 
  geom_col(width = 1, position="fill") +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5),
        axis.text = element_blank())+
  labs(fill="type", 
       x=NULL, 
       y=NULL, 
       title="Types of success")+
  facet_grid(~Var2)+
  geom_label_repel(aes(label=round(value,3)),nudge_x = 1,show.legend = T)+
  coord_polar(theta = "y", start=0)

#########################################################################
#########################################################################
#########################################################################
output = paste("/b06x-isi/b062/a-c/Braun_CRISPR_2D/processed_sequencing_data/",run,"/",type,"/graphics/summary_graphics_",type,".pdf",collapse = "",sep="")
pdf(output)
list(graphics)
dev.off()

