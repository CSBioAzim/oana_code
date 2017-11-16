require(cba)

#split column by a delimiter
split_keep_one_part=function(v,delimiter,keep){
   desired_v=data.frame(do.call('rbind', strsplit(as.character(v),delimiter,fixed=TRUE)))
   return(desired_v[,keep])
}

keep_best=function(dataset,column,decreasing,use_abs_value,dup_col){
   data[,column]=as.numeric(as.character(data[,column]))
   if (use_abs_value==TRUE) {
      data[,column]=abs(as.numeric(as.character(data[,column])))
   }
   sorted=data[order(data[,column]),]
   if (decreasing==TRUE){
      sorted=data[order(data[,column],decreasing=TRUE),]
   }
   dup=which(duplicated(sorted[,dup_col]))
   if (length(dup)>0){
      sorted=sorted[-dup,]
   }
   return(sorted)
}

optimal_ordering=function(m,meth){
  if (dim(m)[1]<=2){
  return (m)
  }
  require(cba)
  d=dist(as.matrix(m),method=meth)
  hc=hclust(d)
  co=order.optimal(d, hc$merge)
  m.optimalRows=as.matrix(m)[co$order,]
  return(m.optimalRows)
}
					  