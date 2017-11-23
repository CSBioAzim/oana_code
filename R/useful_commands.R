require(cba)
require(ggplot2)

#=============
# nice plotting for ggplot
#=============
niceggplot=theme(panel.border = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

niceggplot_fullborders=theme(axis.text=element_text(size=20),
                             axis.title=element_text(size=20))

#=============

window_bed=function(bed1,bed2,w,out){
  write.table(bed1,file=paste(out,'1',sep=''),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
  write.table(bed2,file=paste(out,'2',sep=''),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
  cmd=paste("bedtools window -w ",w," -a ",paste(out,'1',sep='')," -b ",paste(out,'2',sep=''),sep='')
  df=data.frame(fread(cmd,sep='\t'))
  system(paste("rm ",out,'1',sep=''))
  system(paste("rm ",out,'2',sep=''))
  return(df)
}

closest_bed=function(bed1,bed2,w,out){
  write.table(bed1,file=paste(out,'1',sep=''),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(bed2,file=paste(out,'2',sep=''),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
      cmd=paste("bedtools window -w ",w," -a ",paste(out,'1',sep='')," -b ",paste(out,'2',sep=''),sep='')
        df=data.frame(fread(cmd,sep='\t'))
	  system(paste("rm ",out,'1',sep=''))
	    system(paste("rm ",out,'2',sep=''))
	      return(df)
	      }
	      

#split column by a delimiter
split_keep_one_part=function(v,delimiter,keep){
   desired_v=data.frame(do.call('rbind', strsplit(as.character(v),delimiter,fixed=TRUE)))
   return(desired_v[,keep])
}

keep_best=function(data,column,decreasing,use_abs_value,dup_col){
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
					  