import pandas as pd
import numpy as np
path_asso = ""
path = "FINAL/"
file = "results.SPARK.withrawread.tsv"
asso_file = "AF.allpools.SPARK.tsv"
asso = pd.read_csv(path_asso+asso_file,sep="\t")
asso['BinaryRN'] = asso['ReplicateNumber'].apply(lambda x: 1 if x==1 else 10 if x==2 else 100)

"""def getRN_id(x):
    res = 0
    for i in x['BinaryRN']:
        res = res | i
    return res

asso.groupby(['AddedField1']).apply(lambda x: getRN_id(x))"""

df = pd.read_csv(path+file,sep="\t")
inpf = df['input_file'].unique()
dict_inpf={}
dict_af1={}
dict_pool={}
dict_expid={}
dict_brn={}
for c in inpf:
    if c.startswith("/opt"):
        assof = c.split("/")[-2]
        asso_col = assof

        col = c.split("/")[-1]
        col = ".".join(col.split(".")[0:2])
        #record = asso[(asso['TagID'] == col)6 & (asso.CompleteAmplificationID.str.contains(asso_col))].iloc[0]
        record = asso[(asso['TagID'] == col) & (asso['concatenatePoolIDSeqRun'] == assof) ].iloc[0]
        colf = record['CompleteAmplificationID']
        colf_af1 = record['AddedField1']
        print(c, colf, colf_af1)
        dict_inpf[c] = colf
        dict_af1[c] = colf_af1
        dict_expid[c] = record['ExperimentID'].split(".")[0]
        dict_brn[c] = record['BinaryRN']

df['CompleteAmplificationID'] = df['input_file'].map(dict_inpf)
df['AddedField1'] = df['input_file'].map(dict_af1)
df['integration_locus'] = df.apply(lambda x: x['integration_locus'] if x['gap'] > 0 else (
    x['integration_locus'] - x['gap'] if x['target_strand'] == '+' else
    x['integration_locus'] + x['gap'])
                                   , axis=1)
df['Expo'] = df['input_file'].map(dict_expid)
df['System'] = df['AddedField1'].apply(lambda x: x.split("-")[-1])
df['BinaryRN'] = df['input_file'].map(dict_brn)
#df=df.replace({"input_file":dict_inpf})
df['target_start_r2']=df['pos_r2']
df['target_end_r2']=df['target_start_r2']+df['target_matches_r2']
df['delta']=df.apply(lambda x: x['target_end_r2']-x['target_end'] if x['target_strand']=='+' else x['target_start']-x['target_start_r2'],axis=1)
df['f_value']=df.apply(lambda x: x['target_matches']+x['delta'] if x['delta']>0 else x['target_matches'],axis=1)

df['partial_id'] = df.input_file.apply(lambda x: x.split(".",1)[1])
df = df[df.f_value>=30]
df = df[(df.gap<=120) & (df.gap>=-120)].copy()
#df2 = df[df.target_matches<=40]
df['old_junction_locus'] = df['junction_locus']

dict_nucleo = {"A":"T","T":"A","C":"G","G":"C"}
# df2 = df[df.target_matches<=40]
# Ferrari ITR AAV 1	120 2764	2894
itr5_start = 674
itr5_end = 729
itr3_start = 5297
itr3_end = 5352
"""
for i in range(729-674+1):
    if g[674+i]!=dict_nucleo[g[5352-i]]:
        print(674+i)
    else:
        print(g[674+i],":",g[5352-i])
"""
#676->5350
df['junction_locus'] = df.apply(
    lambda x: itr3_end - (x["junction_locus"] - itr5_start)  if (itr5_start <= x["junction_locus"] <= itr5_end) and (("Rev-expo" not in x["Expo"])  and ("Rexp" not in x["Expo"])) else x["junction_locus"], axis=1 )

itr5_start = 730
itr5_end = 792
itr3_start = 5212
itr3_end = 5274


#735->5269

df['junction_locus'] = df.apply(
    lambda x: itr3_end - (x["junction_locus"] - itr5_start)  if (itr5_start <= x["junction_locus"] <= itr5_end) and (("Rev-expo" not in x["Expo"])  and ("Rexp" not in x["Expo"])) else x["junction_locus"], axis=1 )


df['old_aav_last_strand'] = df['aav_last_strand']
inverse_strand = {'+': '-', '-': '+'}
df['aav_last_strand'] = df.apply(
    lambda x: inverse_strand[x['aav_last_strand']] if x['junction_locus'] != x['old_junction_locus'] else x[
        'aav_last_strand'], axis=1)


df = df[df.aav_matches>=30].copy()
df = df[df.start_aav==0].copy()

df.to_csv(path+"results.SPARK.withrawread.modif.tsv",sep="\t", index=False)

