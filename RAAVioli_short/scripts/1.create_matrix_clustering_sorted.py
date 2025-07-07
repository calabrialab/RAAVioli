import pandas as pd
import numpy as np

from scipy.cluster import *
import argparse


def pairDistance(s,l):
    dist = abs(s[0]-l[0])
    for i in range(1,len(s)):
        dist = max(dist,abs(s[i]-l[i]))
    return dist


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file')
parser.add_argument('-b', '--basename', help='Output basename')
parser.add_argument('-o', '--output', help='output path')

args = parser.parse_args()
file = args.input
path = args.output
outputbn = args.basename


df = pd.read_csv(file, sep="\t")
df['shs_end'] = df.apply(lambda x: x['pos_end_r2'] if x['strand_r2']=="-" else x['pos_r2'],axis=1)
merge_id_col = "AddedField1"
cluster_columns = ['target_chr','n_aav_aln', 'target_strand', merge_id_col]
df_comb = df[cluster_columns].drop_duplicates()
df_comb = df_comb.sort_values(cluster_columns).reset_index(drop=True)
df_group = pd.DataFrame(columns=['name','group'])
list_ids_df = df[merge_id_col].unique()
dict_ids = dict(zip(list_ids_df,range(len(list_ids_df))))

for index, row in df_comb.iterrows():
    targ_chr = row['target_chr']
    n_aav_aln = row['n_aav_aln']
    target_strand = row['target_strand']
    id_col = row[merge_id_col]
    dfx = df[(df['target_chr']==targ_chr) & (df['n_aav_aln']==n_aav_aln) & (df[merge_id_col]==id_col)
             & (df['target_strand']==target_strand)].copy()
   # dfx = dfx.drop('group',axis=1)
    dfx['cluster_list'] = dfx[['integration_locus','junction_locus', 'gap']].values.tolist()
    '''
    Here we create a matrix of unique cluster_list values. 
    Since we have a lot of PCR duplicates, the unique here speed up the cluster.
    We then assign back the cluster to the original matrix
    '''
    dfx_unique = dfx.drop_duplicates(['integration_locus','junction_locus', 'gap']).reset_index(drop=True)
    ll = np.array(dfx_unique['cluster_list'].to_list())
    #ll = np.array(dfx['cluster_list'].to_list())
    if ll.shape[0]>1:
        res = hierarchy.fclusterdata(ll, criterion='distance', metric=pairDistance, method='complete', t=20)
    else:
        res=[1]
    dfx_unique['group'] = res
    dfx = dfx.merge(dfx_unique[['integration_locus','junction_locus', 'gap','group']],on=['integration_locus','junction_locus', 'gap'])
    dfx['group'] = dfx['group'].apply(lambda x: f'chr:{targ_chr}_strand:{target_strand}_naav:{n_aav_aln}_idcol:{dict_ids[id_col]}_gid:{x}')
    df_group = pd.concat([df_group,dfx[['name','group']]])
    #print(targ_chr,n_aav_aln,dfx.shape[0],max(res))

df_shape = df.shape
df = df.merge(df_group[['name', 'group']], on='name', how='left')
df_group[['name','group']].to_csv(path+outputbn+".cluster.groups.namereads.tsv",sep="\t",index=False)
gr_cols = ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","CompleteAmplificationID",
 "AddedField1", "aav_last_strand","n_aav_aln"]
dfgroupis = df.groupby(gr_cols+['group'])['shs_end'].nunique().reset_index()
#dfshscid = df.groupby(['group','CompleteAmplificationID'])['shs_end'].nunique().reset_index() #.agg({'seq_count':'size','shs_count':pd.Series.nunique})
df_shs = dfgroupis.groupby(['group'])['shs_end'].sum().reset_index()
df_shs.columns = ['group','shs_count']
df_seq_count = df.groupby(['group'])['shs_end'].size().reset_index()
df_seq_count.columns = ['group','seq_count']

df_maxg = df.groupby(['target_chr','target_strand','integration_locus','gap','junction_locus','n_aav_aln','group',merge_id_col]).size().reset_index(name='ssize').sort_values('ssize',ascending=False).drop_duplicates('group')
df_est = df_shs.copy()
df_est = df_est.merge(df_seq_count)
#df_maxg = df.sort_values(['total_matches','gap'], ascending=False).drop_duplicates('group').sort_index()
df_maxg['IS_genomicID'] = df_maxg.apply(lambda x: f'{x["target_chr"]}_{x["integration_locus"]}_{x["target_strand"]}', axis=1)
df_est = df_est.merge(df_maxg[['group','IS_genomicID','gap','junction_locus','n_aav_aln',merge_id_col]],on="group",how="left")
newdf_shs = df_est.pivot(index=['IS_genomicID','gap','junction_locus','n_aav_aln'], columns=merge_id_col, values='shs_count').reset_index()
df = newdf_shs.copy()
orig_columns =list(df.columns[:4])
pcc=df.columns[4:]
df.to_csv(path+f"ShsCount_{outputbn}.CLUSTER20.tsv",sep="\t",index=False)
df[pcc]=df[pcc].replace(1,np.nan)
df['sum']=df[pcc].sum(axis=1)
df = df[df['sum']!=0].copy()
df = df.drop('sum',axis=1)
cols = df[pcc].columns[df[pcc].sum()>0]
orig_columns.extend(cols)
df_final = df[orig_columns]

df_finalnc = ['IS_genomicID','gap','junction','aav_alignments_start_end']
df_final.columns = df_finalnc + df_final.columns[4:].to_list()

"""df_final = df_final.sort_values(['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction'])"""

df_final.to_csv(path+f"ShsCount_{outputbn}.CLUSTER20.cleaned1.tsv",sep="\t",index=False)


newdf_shs = df_est.pivot(index=['IS_genomicID','gap','junction_locus','n_aav_aln'], columns=merge_id_col, values='seq_count').reset_index()
df = newdf_shs.copy()
orig_columns =list(df.columns[:4])
pcc=df.columns[4:]
df.to_csv(path+f"SeqCount_{outputbn}.CLUSTER20.tsv",sep="\t",index=False)
#df[pcc]=df[pcc].replace(1,np.nan)
df['sum']=df[pcc].sum(axis=1)
df = df[df['sum']!=0].copy()
df = df.drop('sum',axis=1)
cols = df[pcc].columns[df[pcc].sum()>0]
orig_columns.extend(cols)
df_final = df[orig_columns]

df_finalnc = ['IS_genomicID','gap','junction','aav_alignments_start_end']
df_final.columns = df_finalnc + df_final.columns[4:].to_list()

"""df_final = df_final.sort_values(['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction'])"""

df_final.to_csv(path+f"SeqCount_{outputbn}.CLUSTER20.cleaned1.tsv",sep="\t",index=False)

