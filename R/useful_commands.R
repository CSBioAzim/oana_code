
#split column by a delimiter
x=c('test,test2','test3,test4')
data.frame(do.call('rbind', strsplit(as.character(x),',',fixed=TRUE)))