#determinate direction
##get direction for time series data
##@ xingsheng yang

rm(list=ls())
##load data
##please run identify_molecular_transformation.R first
link=read.table(".../data/all_reaction_mic0.820none.txt",head=1,sep = "\t")
data1=read.table(".../data/data_st.txt",head=1,sep = "\t",row.names = 1)

source(".../data/get_dir_function.R")#load the function

##set parameters
method="mic"
n_expansion=1000000#The calculation of mic requires that the value is not too low, and that the multiplier does not affect the result.

##calculate
result=get_dir(link,data1,method,n_expansion)

#output result
write.table(result,".../data/all_direction_mic0.820none.txt",sep="\t",quote = F,row.names = F)
