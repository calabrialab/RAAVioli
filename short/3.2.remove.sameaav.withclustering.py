import numpy as np
import pandas as pd
import scipy
from scipy.cluster import *
def pairDistance(s,l):
    dist = abs(s[0]-l[0])
    for i in range(1,len(s)):
        dist = max(dist,abs(s[i]-l[i]))
    return dist

path = "FINAL/post_mod/"
file = "PCR_COUNT.EXTENDED.SPARK.removed_targetgenes.tsv"
dfg = pd.read_csv(path+file,sep="\t")
dfg['list_alignment'] = dfg['aav_alignments_start_end_y'].apply(lambda x: [int(j) for xs in [
    l.split(",")[1:3] if l.split(",")[3] == "+" else l.split(",")[1:3][::-1] for l in x.split(";")] for j in xs])
dfg['list_alignment'] = dfg['list_alignment'].apply(
    lambda x: [x[i] for j in [list(range(0, len(x), 2)), list(range(1, len(x), 2))] for i in j])
dfg['group'] = np.NAN
sample_col = "AddedField1"
for af in dfg[sample_col].unique():
    dfg_af = dfg[dfg[sample_col] == af]

    for i in dfg_af['n_aav_aln'].unique():
        dfg_i = dfg_af[dfg_af['n_aav_aln'] == i].copy()
        k_af = f'Naav:{i};sample:{af}'
        print("K: ", k_af)
        if dfg_i.shape[0]>1:
            ll = np.array(dfg_i['list_alignment'].to_list())
            res = hierarchy.fclusterdata(ll, criterion='distance', metric=pairDistance, method='complete', t=3)
            dfg_i['group'] = res
            for index, row in dfg_i.iterrows():
                dfg.at[index,'final_group'] =f'{k_af};Group:{row["group"]}'
        else:
            for index, row in dfg_i.iterrows():
                dfg.at[index,'final_group'] =f'{k_af};Group:1'

dfg.to_csv(path+"Extend.reargroupcluster3.tsv",sep="\t",index=False)
#df1 = dfg[dfg.n_aav_aln==1].copy().reset_index(drop=True)
#dfrear=dfg[dfg.n_aav_aln>1].sort_values(['ssc','read_length'], ascending=False).drop_duplicates(['final_group']).reset_index(drop=True)
#df_ssc = pd.concat([df1, dfrear])
#df_ssc.to_csv(path+"EXTENDED.SPARK.cleaned_rearrangements.cluster3.highest_ssc.tsv",sep="\t", index=False)

df1 = dfg[dfg.n_aav_aln==1].copy().reset_index(drop=True)
dfrear=dfg[(dfg.n_aav_aln>1)].sort_values(['ssc','read_length'], ascending=False).drop_duplicates(['final_group']).reset_index(drop=True)
dfrear_tokeep =dfg[(dfg.n_aav_aln>1) & (dfg['PCR_count']>1)].copy().reset_index(drop=True)
dfrear_tokeep = dfrear_tokeep[~dfrear_tokeep['new_k'].isin(dfrear['new_k'])].copy().reset_index(drop=True)
df_ssc = pd.concat([df1, dfrear_tokeep, dfrear])
df_ssc.to_csv(path+"EXTENDED.SPARK.cleaned_rearrangements.cluster3.highest_ssc.removed_targetgenes.tsv",sep="\t", index=False)
idx_cols = list(df_ssc.columns[:8])
#idx_cols.extend(["aav_alignments_start_end_y"])
newdf = df_ssc.pivot(index=idx_cols,columns='ID', values='ssc').reset_index()
newdf.to_csv(path+"SPARK.cleaned_rearrangements.cluster3.highest_ssc.removed_targetgenes.tsv",sep="\t", index=False)
