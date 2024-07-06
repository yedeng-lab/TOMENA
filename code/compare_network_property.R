#compare network property with random networks (n=1000)
##@ xingsheng yang
rm(list=ls())
##load data which was generated from iNAP pipeline (https://inap.denglab.org.cn/)
pro=read.table(".../data/Galaxy3-[Individual_nodes__centrality_result].tabular",sep="\t",head=1,row.names = 1)
random=read.table(".../data/random.txt",sep="\t",head=1,row.names = 1)# Reorganized from Galaxy8-[Random_network_properties_with_fast_greedy].tabular
random$pvalue=NA
##statistical test
for (ii in 1:nrow(random)) {
  random$pvalue[ii]=tsum.test(mean.x=random$Mean[ii], s.x=random$Sd[ii], n.x=1000, mu=random$Indexes[ii])$p.value
}
