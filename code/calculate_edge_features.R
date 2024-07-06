#calculate edge features
##@ xingsheng yang
rm(list=ls())

library(igraph)
##load data
##please run identify_molecular_transformation.R first
node=read.table(".../data/all_network_mic0.820none.txt")
##please run determine_direction.R first
edge_o=read.table(".../data/all_reaction_mic0.820none.txt",header = T)#original edge

##1 calculate edge betweenness
##Average betweenness centrality was defined as the number of shortest paths passing through a link

##transform data type
igraph=graph_from_adjacency_matrix(as.matrix(node),mode="undirected",diag =FALSE)
edge=data.frame(as_edgelist(igraph))
edge_be=data.frame(Small_molecule=edge[[1]],Large_molecule=edge[[2]],Betweenness_centrality=log10(edge.betweenness(igraph)))#get log fold betweenness
##merge betweenness and edge information
edge_m=merge(edge_o,edge_be)
##calculate average betweenness for each transformation
mte_sum=data.frame(aggregate(edge_m$Betweenness_centrality,by=list(Molecular_transformation=edge_m$Molecular_transformation),mean),aggregate(edge_m$Betweenness_centrality,by=list(Molecular_transformation=edge_m$Molecular_transformation),sd)[2])
names(mte_sum)[c(2,3)]=c("mean","sd")

##2 calculate edge efficiency contribution
##Efficiency contribution was calculated as the decrease in overall efficiency after removing all edges corresponding to a transformation.

##get transformation list
result=data.frame(transformation=unique(edge_o$Molecular_transformation),link_lack_effect=NA)

##a loop for calculate efficiency contributions
source(".../function/link_effect_function.R")#load function
for (ii in 1:nrow(result)) {
  trans=result$transformation[ii]
  edge2=edge_o[edge_o$Molecular_transformation==trans,]#match edges with transformation
  
  node2=node# make a copy of node
  rid=match(edge2$Small_molecule,rownames(node))
  cid=match(edge2$Large_molecule,colnames(node))
  node2[cbind(c(rid,cid),c(cid,rid))]<-0#delete matched edges in the notwork matrix
  graph_lack=graph_from_adjacency_matrix(as.matrix(node2),mode="undirected",diag =FALSE)
  result$link_lack_effect[ii]=link.effect(igraph,graph_lack)
}

