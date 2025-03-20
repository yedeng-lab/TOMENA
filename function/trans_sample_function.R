##function for counting transformations in each sample
##@ xingsheng yang

##data: data with sample name
##reaction: "_reaction" file generated using "moteam_function"

trans_sample= function (data,reaction) {
  trans_form=unique(reaction$Molecular_transformation)
  nrow=length(trans_form)
  ncol=ncol(data)
  
  result_table=data.frame(matrix(NA, nrow = nrow, ncol = ncol),row.names =trans_form)
  colnames(result_table)=colnames(data)
  
  for (ii in 1:ncol) {
    molecule_list=row.names(data[(data[,ii]>0),])
    reaction_sample=reaction[(reaction$Small_molecule %in% molecule_list&reaction$Large_molecule %in% molecule_list),]
    result_pre=as.data.frame(table(reaction_sample$Molecular_transformation))
    result_table[,ii]=result_pre[match(row.names(result_table),result_pre$Var1),2]  
  }
  result_table[is.na(result_table)]=0
  return(result_table)
}