#identify molecular transformation
##@ xingsheng yang
rm(list=ls())
##load the main code of moteam
source(".../function/moteam_function.R")

##load the data
## To ensure robust time-series data analysis, the average value of each molecule in all replicate bioreactors at each time point was used
## The original source of the data set can be available at https://github.com/yedeng-lab/Microbial_metabolism_during_anaerobic_digestion
data=read.table(".../data/data_st.txt",sep="\t",row.names=1,head=1)

##set up parameters
min_occur=20
p_adj_method="none"
select_num=5
threshold=0.8
method="mic"

##use motem to identify molecular transformation
cal_moteam=moteam(data,n.majority=select_num,threshold=threshold,min_occur=min_occur,p_adj_method=p_adj_method,method = method)

##output the result
setwd(".../data")
write.table(cal_moteam$MTE,paste0("all","_MTE_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,row.names = F) #1 transformation form
write.table(cal_moteam$network_matrix,paste0("all","_network_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,col.names =NA) #2 network matrix
write.table(cal_moteam$relation_table,paste0("all","_reaction_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,row.names = F) #3 pairwise molecules with transformation relationship