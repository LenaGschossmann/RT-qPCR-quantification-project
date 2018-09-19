
order_specific = function(df, colName, order_df){
  
  idxC = which(colnames(df)== colName)
  df_ordered = data.frame(matrix(ncol=ncol(df), nrow=nrow(df)))
  colnames(df_ordered) = colnames(df)
  idxR1 = 1
  for(idxO in order_df){
    idxR2 = idxR1 + nrow(df[df[,idxC] == idxO,]) - 1
    df_ordered[idxR1:idxR2,] = df[df[,idxC] == idxO,]
    idxR1=idxR1 + nrow(df[df[,idxC]==idxO,])
  }
  
  return(df_ordered)
}