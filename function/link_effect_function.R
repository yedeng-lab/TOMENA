#function for calculating edge efficiency contribution
##Efficiency contribution was calculated as the decrease in overall efficiency after removing some edges
##@xingsheng yang
library(igraph)
link.effect <- function(graph_original,graph_lack){
  dd1 <- 1/shortest.paths(graph_original)
  diag(dd1) <- NA
  dd2 <- 1/shortest.paths(graph_lack)
  diag(dd2) <- NA
  efficiency1 <- mean(dd1, na.rm=T)
  efficiency2 <- mean(dd2, na.rm=T)
  link_efficiency_effect=efficiency1-efficiency2
  return(link_efficiency_effect)
}