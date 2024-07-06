#Molecular Transformation Enumerator based on Advantageous Mass Difference (MoTEAM)
##function for Molecular Transformation Enumerator
##"C","H","O","N","P","S"
##@xingsheng yang
library(vegan)
library(minerva)
library(reshape2)

##n.majority: Molecules with a frequency of occurrence greater than or equal to this value are retained
##threshold: Correlations greater than this value are considered to have a relationship
##p_adj_method: Methods for calibrating p-values
##min_occur: Molecular difference with frequency lower than this value do not enter subsequent consideration
##method: mic, pearson, spearman are available now
##n_expansion: Only for 'mic'. The calculation of mic requires that the value is not too low, and that the multiplier does not affect the result
##only_generate_cor: If TRUE, only generate correlation matrix

moteam=function (data,n.majority=5,threshold=0.6,p_adj_method="none",min_occur=10,method="mic",n_expansion=1000000,only_generate_cor=FALSE){
  output=list()
  #step 1 majority
  #rownames(data)=toupper(rownames(data))
  test=data
  test[data!=0]=1
  data.maj=data[rowSums(test)>=n.majority,]
  #data.maj=data[rowSums(vegan::decostand(data,"pa"))>=n.majority,]
  
  #step 2 get element composition
  ## modified from https://github.com/danczakre/ICRTutorial
  CHEMICAL_ELEMENTS = c("C","H","O","N","P","S")
  molecularFormula <- row.names(data.maj)
  numericalFormula <- array(0, dim=c(length(molecularFormula), length(CHEMICAL_ELEMENTS)))
  
  for (k in 1:length(molecularFormula)){
    formula <- molecularFormula[k]
    ge <- gregexpr("[A-Z]\\d*", formula, perl=TRUE)
    s_index <- ge[[1]]
    s_len <- attr(s_index, "match.length")
    for (i in 1:length(s_len)){
      token <- substr(formula, s_index[i], s_index[i] + s_len[i] - 1)
      element <- substr(token, 1, 1)
      if (grepl(element, paste(CHEMICAL_ELEMENTS,collapse = ""))) {
        idx = which(CHEMICAL_ELEMENTS %in% element)
        if (numericalFormula[k, idx] > 0) next 
        if (s_len[i] == 1) {
          numericalFormula[k, idx] = 1
        } else {
          numElement <- try(strtoi(substr(formula, s_index[i] + 1, s_index[i] + s_len[i] - 1)))
          if (class(numElement)=="integer"){
            numericalFormula[k, idx] = numElement
          } else {
            print(paste("[ERROR] an unknown chemical element found:", token, "in", formula))
          }
        }
      } else {
        print(paste("[ERROR] an unknown chemical element found:", element, "in", formula))
      }
    }
  }
  
  colnames(numericalFormula)=CHEMICAL_ELEMENTS
  rownames(numericalFormula)=molecularFormula
  
  #step 3 calculate molecular weight
  options(digits=10)
  element_mass=c(12.000000,1.007825,15.994915,14.003074,30.973762,31.972071)
  names(element_mass)=c("C","H","O","N","P","S")
  masscalFormula <- as.data.frame(t(apply(numericalFormula, 1, function(row) row * element_mass)))
  masscalFormula$molecular_weight=rowSums(masscalFormula)
  
  #step 4 order according to molecular weight
  orderdecrease=order(masscalFormula$molecular_weight,decreasing = F)
  data.maj_o=data.maj[orderdecrease,]
  numericalFormula_o=numericalFormula[orderdecrease,]
  
  #step 5 calculate correlation matrix
  if (method== "mic") {
    data.maj_o_e=data.maj_o*n_expansion #expansion for MIC
    mic=data.frame(matrix(NA,nrow=nrow(data.maj_o_e),ncol=nrow(data.maj_o_e)))
    for (i in 1:(nrow(data.maj_o_e)-1)) {
      for (j in (i+1):nrow(data.maj_o_e)) {mic[i,j]=minerva::mine(as.numeric(data.maj_o_e[i,]),as.numeric(data.maj_o_e[j,]),alpha=1,na.rm =F)$MIC}}
    rownames(mic)=colnames(mic)=rownames(data.maj_o_e)
    plot1=plot(density(mic[!is.na(mic)]),main="(A) MIC value distribution")
  } else if (method == "pearson") {
    data.maj_o_e=data.maj_o
    mic=cor(t(data.maj_o_e),method = method)
    mic[lower.tri(mic,diag = TRUE)]=NA
    mic=as.data.frame(mic)
    plot1=plot(density(abs(mic[!is.na(mic)])),main="(A) Pearson value distribution")
    } else if (method == "spearman") {
      data.maj_o_e=data.maj_o
      mic=cor(t(data.maj_o_e),method = method)
      mic[lower.tri(mic,diag = TRUE)]=NA
      mic=as.data.frame(mic)
      plot1=plot(density(abs(mic[!is.na(mic)])),main="(A) spearman value distribution")
  } else stop(warning(paste0 ("method ",method," unavailable"))) #only tree method are available
  
  out_matrix=mic #for output correlation matrix
  if (!only_generate_cor) {
    mic[upper.tri(mic,diag = FALSE)]=abs(mic[upper.tri(mic,diag = FALSE)])
    #step 6 get mass difference matrix
    dif_mass=data.frame(matrix(NA,nrow=nrow(data.maj_o_e),ncol=nrow(data.maj_o_e)))
    rownames(dif_mass)=colnames(dif_mass)=rownames(data.maj_o_e)
    for (i in 1:(nrow(data.maj_o_e)-1)) {
      for (j in (i+1):nrow(data.maj_o_e)) {
        dif_table=numericalFormula_o[j,]-numericalFormula_o[i,]
        dif_table=dif_table[!dif_table==0]
        dif_mass[i,j]=paste(paste0(names(dif_table),dif_table),collapse = "")
      }}
    
    #step 7 filter transformation with less frequency
    node_matrix=matrix(0,nrow=nrow(data.maj_o_e),ncol=nrow(data.maj_o_e))#save node message 0/1
    node_matrix[upper.tri(node_matrix)]=1
    
    extractfromtable=function (tableresult) {
      tableresult=as.data.frame(tableresult, row.names = NULL,
                                responseName = "Freq", stringsAsFactors = F,
                                sep = "", base = list(LETTERS))
      tableresult[2]=as.numeric(unlist(tableresult[2]))
      return(tableresult)
    }#a function resort data frame after table()
    
    trans_table1=extractfromtable(table(reshape2::melt(dif_mass,measure.vars=colnames(dif_mass),na.rm = T)[2]))
    plot2=plot(density(trans_table1[,2]),main="(B) Mass difference distribution")
    
    obtain_table1=trans_table1[trans_table1[2]>=min_occur,1]
    search_1=apply (dif_mass, c(1,2), function(element) mapply(function(x) x %in% obtain_table1, element))#contain nodes coordinate
    
    node_matrix[!search_1]=0#get matrix delete transformation with less frequency
    mic[!search_1]=NA
    dif_mass[!search_1]=NA
    
    #step 8 evaluate Advantageous Mass Difference based on hypergeometric distribution
    trans_table_o=extractfromtable(table(reshape2::melt(dif_mass,measure.vars=colnames(dif_mass),na.rm = T)[2]))#overall
    
    dif_mass_mic=dif_mass
    dif_mass_mic[!mic>=threshold]=NA
    trans_table_s=extractfromtable(table(reshape2::melt(dif_mass_mic,measure.vars=colnames(dif_mass_mic),na.rm = T)[2]))#sample with higher MIC
    
    knumber=sum(trans_table_s[2])#number of sample with higher MIC (k)
    onumber=sum(trans_table_o[2])#over number of sample with higher MIC (m+n)
    high_proportion=sum(trans_table_s[2])/sum(trans_table_o[2])#proportion of sample with higher MIC
    
    trans_table_s$pvalue=0
    for (ii in 1:nrow(trans_table_s)) {
      qnumber=trans_table_s[ii,2]#sample number for a transformation
      mnumber=trans_table_o[match(trans_table_s[ii,1],trans_table_o[,1]),2]#number for the transformation in overall
      nnumber=onumber-mnumber#number not the transformation in overall
      trans_table_s$pvalue[ii]=phyper(qnumber-1,mnumber,nnumber,knumber,lower.tail = F)
    }
    trans_table_s$adjust_pvalue=p.adjust(trans_table_s$pvalue,method = p_adj_method)#use selected method adjust p-value
    
    obtain_table2=trans_table_s[trans_table_s$adjust_pvalue<0.05,1]#adjusted p value < 0.05
    
    if (length(obtain_table2)<1) stop (warning("Not found Molecular Transformation"))
    
    search_2=apply (dif_mass_mic, c(1,2), function(element) mapply(function(x) x %in% obtain_table2, element))#contain nodes coordinate
    
    dif_mass_mic[!search_2]=NA#only Advantageous Mass Difference contain
    node_matrix[is.na(dif_mass_mic)]=0
    mic[is.na(dif_mass_mic)]=NA
    
    #step 9 output file generation
    ##Molecular Transformation Enumerator result
    MTE=extractfromtable(table(reshape2::melt(dif_mass_mic,measure.vars=colnames(dif_mass_mic),na.rm = T)[2]))
    colnames(MTE)=c("Molecular_transformation","Frequnecy")
    
    ##network matrix
    rownames(node_matrix)=colnames(node_matrix)=rownames(mic)#add name
    node_matrix=node_matrix+t(node_matrix)#Symmetric matrix
    network_matrix=node_matrix[rowSums(node_matrix)>0,rowSums(node_matrix)>0]#delete node without correlation
    
    ##relation table (edge information)
    relation_table=cbind(reshape2::melt(cbind(row.names(dif_mass_mic),dif_mass_mic),measure.vars=colnames(dif_mass_mic),na.rm = F),reshape2::melt(cbind(row.names(mic),mic),measure.vars=colnames(dif_mass_mic),na.rm = F)[3])#reshape and merge data
    colnames(relation_table)=c("Small_molecule","Large_molecule","Molecular_transformation","Correlation_index")
    relation_table=relation_table[!is.na(relation_table$Molecular_transformation),]
    
    ###output
    output$MTE=MTE
    output$network_matrix=network_matrix
    output$relation_table=relation_table
  }
  
  output$correlation_matrix=out_matrix
  return(output)
}




