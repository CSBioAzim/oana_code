
import numpy as np
import pandas as pd
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')
sc.logging.print_versions()

def analyze_feature_vs_feature(adata_here,f1='louvain',f2='batch',corrected_pval=0.05):
    sc.pl.umap(adata_here, color=[f1,f2])
    import scipy

    f1s=list(set(adata_here.obs[f1]))
    f2s=list(set(adata_here.obs[f2]))

    oddsratios=np.zeros((len(f1s),len(f2s)))
    pvals=np.zeros((len(f1s),len(f2s)))
    proportions=np.zeros((len(f1s),len(f2s)))

    for f1_idx in range(len(f1s)):
        f1_here=f1s[f1_idx]
        cells_in_f1=list(adata_here.obs_names[adata_here.obs[f1]==f1_here])
        for f2_idx in range(len(f2s)):
            f2_here=f2s[f2_idx]
            cells_in_f2=list(adata_here.obs_names[adata_here.obs[f2]==f2_here])

            total=list(adata_here.obs_names)
            overlap=list(set(cells_in_f1).intersection(set(cells_in_f2)))
            contingency_table=np.array([[len(overlap),len(cells_in_f1)-len(overlap)],[len(cells_in_f2)-len(overlap),0]])
            contingency_table[1,1]=len(total)-contingency_table[0,0]-contingency_table[1,0]-contingency_table[0,1]
            oddsratio, pvalue = scipy.stats.fisher_exact(contingency_table)
            if pvalue<=corrected_pval/(1.0*(len(f1s)*len(f2s))):
                ps=0.0000000000001
                oddsratios[f1_idx,f2_idx]=np.log2(oddsratio+ps)
                pvals[f1_idx,f2_idx]=-1.0*np.sign(oddsratio)*np.log10(pvalue+ps)
            proportion_cells_in_f1_from_f2=1.0*len(overlap)/len(cells_in_f1)
            proportions[f1_idx,f2_idx]=proportion_cells_in_f1_from_f2

    oddsratios_df=pd.DataFrame(oddsratios)
    oddsratios_df.index=f1s
    oddsratios_df.columns=f2s

    pvals_df=pd.DataFrame(pvals)
    pvals_df.index=f1s
    pvals_df.columns=f2s

    proportions_df=pd.DataFrame(proportions)
    proportions_df.index=f1s
    proportions_df.columns=f2s

    import seaborn as sns
    x=4
    sns.clustermap(proportions_df.T,vmin=0, vmax=1,cmap='bwr')
    sns.clustermap(oddsratios_df.T,vmin=-x, vmax=x,cmap='bwr')


def annotate_cell_cycle_scores_human(adata_results_file,cell_cycle_file='/ahg/regevdata/users/oursu/code/general_data/cellcycle/regev_lab_cell_cycle_genes.txt'):
    cell_cycle_genes = [x.strip() for x in open(cell_cycle_file)]

    s_genes  =  cell_cycle_genes [:43]
    g2m_genes = cell_cycle_genes[43:]
    
    adata_cellcycle = sc.read(adata_results_file+'.basic.h5ad')
    adata_cellcycle
    
    sc.pp.log1p(adata_cellcycle)
    sc.pp.scale(adata_cellcycle)
    
    sc.tl.score_genes_cell_cycle(adata_cellcycle, s_genes=s_genes, g2m_genes=g2m_genes)
    
    adata_annotated=sc.read(adata_results_file+'.basic.h5ad')
    #now, assign the cell cycle scores from the adata_cellcycle to adata_annotated
    s_scores=adata_cellcycle.obs['S_score']
    g2m_scores=adata_cellcycle.obs['G2M_score']
    adata_annotated.obs['S_score_added']=s_scores.loc[adata_annotated.obs_names]
    adata_annotated.obs['G2M_score_added']=g2m_scores.loc[adata_annotated.obs_names]
    adata_annotated

    adata_annotated.write(adata_results_file+'.basic.cc.h5ad')

def compute_TPT(gbcs_dataset):
    
    '''
    input: pandas data frame with the columns "cbc", "umi", "gbc", "r2" where every row is a read 
    output: pandas data frame with the columns "gbc", "cbc", "umi", "cbc-umi-r2-count", "cbc-umi-count", "TPT"
    NOTE: for the input, multiple reads corresponding to the same cbc-umi combination should be listed as separate lines!
    '''

    import copy
    import re
    import time

    print("======== annotating cbc-umi pairs, and cbc-umi-r2")
    gbcs_dataset['cbcumi']=[x+'-'+y for x,y in zip(gbcs_dataset['cbc'],gbcs_dataset['umi'])]
    cbcumir2=list([x+'+'+y for x,y in zip(gbcs_dataset['cbcumi'],gbcs_dataset['r2'])])
    gbcs_dataset['cbcumir2']=list([x+'_gbc_'+y for x,y in zip(cbcumir2,gbcs_dataset['gbc'])])

    print("======== counting the numbers of reads supporting cbc-umi-r2 and for denominator cbc-umi")
    cbcumi_group=gbcs_dataset.groupby('cbcumi').size()
    cbcumi_r2_group=gbcs_dataset.groupby('cbcumir2').size()
    cbcumi_from_grouped_reads=[x.split('+')[0] for x in cbcumi_r2_group.index]

    print("======== computing TPT")
    #divide every cbc-umi-r2 value by the cbc-umi values 
    combo_counts=pd.DataFrame({'cbc-umi-r2':cbcumi_r2_group,
                              'cbc-umi-r2-name':cbcumi_r2_group.index,
                              'cbc-umi':cbcumi_from_grouped_reads})
    combo_counts['cbc-umi-total']=copy.deepcopy(list(cbcumi_group.loc[combo_counts['cbc-umi']]))
    combo_counts['TPT']=1.0*combo_counts['cbc-umi-r2']/combo_counts['cbc-umi-total']
    
    combo_counts['gbc']=list([x.split('_gbc_')[1] for x in list(combo_counts.loc[:,'cbc-umi-r2-name'])])
    combo_counts['umi']=list([x.split('-')[1] for x in list(combo_counts['cbc-umi'])])
    combo_counts['cbc']=list([x.split('-')[0] for x in list(combo_counts['cbc-umi'])])
    
    print("======== compiling the final result")
    to_return=pd.DataFrame({'gbc':combo_counts['gbc'],
                           'cbc':combo_counts['cbc'],
                           'umi':combo_counts['umi'],
                           'cbc-umi-r2-count':combo_counts['cbc-umi-r2'],
                           'cbc-umi-count':combo_counts['cbc-umi-total'],
                           'TPT':combo_counts['TPT']})
    to_return=to_return.reset_index(drop=True)
    return(to_return)


