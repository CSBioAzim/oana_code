
import argparse
import sys
import os
import copy
import re
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+'/utils/')
import oana_utils 

def main():
    parser = argparse.ArgumentParser(description='Join files by a column. Similar to unix join, but with some extra functionality.')
    parser.add_argument('--filelist',help='Comma-delimited list of files')
    parser.add_argument('--join_columns',help='Comma-delimited list of columns by which all files should be joined. 1-based.')
    parser.add_argument('--empty',help='What to write for entries that are missing from a given file')
    parser.add_argument('--inner_or_outer',default='inner',help='"inner" or "outer" join. By default, this returns the inner join')
    parser.add_argument('--delim',default='\t')
    parser.add_argument('--out')
    args = parser.parse_args()

    dictionaries={}
    dictionary_of_num_cols={}
    dictionary_of_cols={}
    filelist=args.filelist.split(',')
    columns=args.join_columns.split(',')
    entries=set()
    for f_idx in range(len(filelist)):
        cur_file=filelist[f_idx]
        cur_col=int(columns[f_idx])-1
        d,num_cols,indexed_cols=oana_utils.file_to_dictionary_by_column(cur_file,cur_col,args.delim)
        dictionaries[cur_file]=d
        dictionary_of_num_cols[cur_file]=num_cols
        dictionary_of_cols[cur_file]=indexed_cols
        cur_entries=set(dictionaries[cur_file].keys())
        if args.inner_or_outer=='inner':
            if f_idx==0:
                entries=copy.deepcopy(cur_entries)
            else:
                entries=copy.deepcopy(entries.intersection(cur_entries))
        if args.inner_or_outer=='outer':
            entries=copy.deepcopy(entries.union(cur_entries))

    out=oana_utils.open_generic(args.out,"w")

    #write column names
    collist=[]
    for f_idx in range(len(filelist)):
        cur_file=filelist[f_idx]
        for c_idx in range(dictionary_of_num_cols[cur_file]):
            collist.append(dictionary_of_cols[cur_file][c_idx])
    out.write("#"+args.delim.join(collist)+'\n')
            
    #write rows
    for entry in entries:
        list_of_lines_for_this_entry=[""]
        for f_idx in range(len(filelist)):
            cur_file=filelist[f_idx]
            new_list_of_lines_for_this_entry=[]
            if entry in dictionaries[cur_file]:
                times=len(dictionaries[cur_file][entry])
                for t in range(times):
                    for available_line in list_of_lines_for_this_entry:
                        new_list_of_lines_for_this_entry.append(args.delim.join([available_line,dictionaries[cur_file][entry][t]]))
            else:
                for available_line in list_of_lines_for_this_entry:
                    new_list_of_lines_for_this_entry.append(args.delim.join([available_line,args.delim.join([args.empty] * dictionary_of_num_cols[cur_file])]))
            list_of_lines_for_this_entry=copy.deepcopy(new_list_of_lines_for_this_entry)

        for entryline in list_of_lines_for_this_entry:
            out.write(re.sub("^"+args.delim,'',entryline,flags=re.M)+'\n')
    out.close()
    
main()


