import gzip

def file_to_dictionary_by_column(fname,col,delim='\t'):
    d={}
    f_opened=open_generic(fname,"r")
    indexed_cols,num_cols=index_columns(f_opened,delim)
    f_opened.close()

    f_opened=open_generic(fname,"r")
    for line in f_opened.readlines():
        if line[0]!="#":
            items=line.strip().split(delim)
            if len(items)!=num_cols:
                print "warning: encountered different number of columns. expected "+str(num_cols)+" but got "+str(len(items))
            item=items[col]
            if item not in d:
                d[item]=[]
            d[item].append(line.strip())
    f_opened.close()
    return d,num_cols,indexed_cols
        
def index_columns(opened_file,delim):
    column_dict={}
    
    first_line=opened_file.readlines()[0]
    num_cols=len(first_line.strip().split(delim))
    column_list=first_line.strip().split(delim)
    for i in range(num_cols):
        if first_line[0]!="#":
            column_dict[i]="col"+str(i)
        else:
            column_dict[i]=column_list[i]
    return column_dict,len(column_list)
    
def open_generic(fname,type_of_open):
    if fname.endswith(".gz"):
        return gzip.open(fname,type_of_open)
    else:
        return open(fname,type_of_open)
