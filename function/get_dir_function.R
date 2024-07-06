##function for getting direction for time series data
##@ xingsheng yang
library(minerva)

##link: link data from moteam
##data1: time-series data
##method: mic, pearson, spearman are available now
##n_expansion: Only for 'mic'. The calculation of mic requires that the value is not too low, and that the multiplier does not affect the result

get_dir=function(link,data1,method,n_expansion) {
  link2=link #a copy of link
  link2$small2large=0#original set
  link2$large2small=0
  link2$direction="uncertain"
  data2=data1*n_expansion#
  
  ##compare correlation values to determinate the directions
  ##small2large means smaller molecule changed first
  ##large2small means larger molecule changed first
  if (method=="pearson"|method=="spearman") {
    small_1=data1[link$Small_molecule,1:(ncol(data1)-1)]
    large0=data1[link$Large_molecule,2:ncol(data1)]
    small0=data1[link$Small_molecule,2:ncol(data1)]
    large_1=data1[link$Large_molecule,1:(ncol(data1)-1)]
    
    for (ii in 1:nrow(link2)) link2$small2large[ii]= cor(t(small_1[ii,]),t(large0[ii,]),method = method)
    
    for (ii in 1:nrow(link2)) {
      link2$large2small[ii]= cor(t(small0[ii,]),t(large_1[ii,]),method = method)
      if (abs(link2$large2small[ii])>abs(link2$small2large[ii]) ) link2$direction[ii]="decrease"
      else if (abs(link2$large2small[ii])<abs(link2$small2large[ii]) ) link2$direction[ii]="increase"
    }
  } else if (method=="mic") {
    small_1=data2[link$Small_molecule,1:(ncol(data1)-1)]
    large0=data2[link$Large_molecule,2:ncol(data1)]
    small0=data2[link$Small_molecule,2:ncol(data1)]
    large_1=data2[link$Large_molecule,1:(ncol(data1)-1)]
    for (ii in 1:nrow(link2)) link2$small2large[ii]=minerva::mine(as.numeric(small_1[ii,]),as.numeric(large0[ii,]),alpha=1,na.rm =F)$MIC
    for (ii in 1:nrow(link2)) {
      link2$large2small[ii]= minerva::mine(as.numeric(small0[ii,]),as.numeric(large_1[ii,]),alpha=1,na.rm =F)$MIC
      if (abs(link2$large2small[ii])>abs(link2$small2large[ii]) ) link2$direction[ii]="decrease"
      else if (abs(link2$large2small[ii])<abs(link2$small2large[ii]) ) link2$direction[ii]="increase"
    }
  }
  return (link2)
}
