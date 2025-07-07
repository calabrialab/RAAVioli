import sys

import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input tsv file with all IS')
parser.add_argument('-a', '--assofile', help='assofile')
parser.add_argument('-o', '--output', help='Output File')
parser.add_argument('-d', '--dfitr', help='tsv with itr starts and ends',default=None)
parser.add_argument('-n', '--poolname', help='Pool Name',default=None)


args = parser.parse_args()
file = args.input
asso_file = args.assofile
outputbn = args.output
itr_file = args.dfitr
pool_name = args.poolname


asso = pd.read_csv(asso_file,sep="\t")



df = pd.read_csv(file,sep="\t")
df['TagID'] = df['input_file'].apply(lambda x: ".".join(x.split("/")[-1].split(".")[:2]))
df['concatenatePoolIDSeqRun'] = pool_name

df_shape_bm = df.shape
df = df.merge(asso[['TagID','CompleteAmplificationID','concatenatePoolIDSeqRun','AddedField1','ExperimentID','ReplicateNumber']],
         on=['TagID','concatenatePoolIDSeqRun'])

df_shape_new = df.shape
if df_shape_new[0]!=df_shape_bm[0] or (df_shape_new[1]-4)!=df_shape_bm[1]:
    print("SOMETHING WRONG IN THE MERGE, NEW SHAPE DOES NOT MATCH, check poolname is equal to concatenateseqidrun")
    sys.exit(-1)

df['integration_locus'] = df.apply(lambda x: x['integration_locus'] if x['gap'] > 0 else (
    x['integration_locus'] - x['gap'] if x['target_strand'] == '+' else
    x['integration_locus'] + x['gap']), axis=1)



df['System'] = df['AddedField1'].apply(lambda x: x.split("-")[-1])

df['target_start_r2']=df['pos_r2']
df['target_end_r2']=df['target_start_r2']+df['target_matches_r2']
df['delta']=df.apply(lambda x: x['target_end_r2']-x['target_end'] if x['target_strand']=='+' else x['target_start']-x['target_start_r2'],axis=1)
df['f_value']=df.apply(lambda x: x['target_matches']+x['delta'] if x['delta']>0 else x['target_matches'],axis=1)

df['partial_id'] = df.input_file.apply(lambda x: x.split(".",1)[1])
df = df[df.f_value>=30]
df = df[(df.gap<=120) & (df.gap>=-120)].copy()
df['old_junction_locus'] = df['junction_locus']

dict_nucleo = {"A":"T","T":"A","C":"G","G":"C"}

if itr_file:
    itr_df = pd.read_csv(itr_file,sep="\t")
    for index_itr, row_itr in itr_df.iterrows():
        itr5_start = row_itr['itr5_start']
        itr5_end = row_itr['itr5_end']
        itr3_start = row_itr['itr3_start']
        itr3_end = row_itr['itr3_end']
        df['junction_locus']=df['junction_locus'].map(lambda x: itr3_end-(x-itr5_start) if itr5_start<=x<=itr5_end else x)




df['old_aav_last_strand'] = df['aav_last_strand']
inverse_strand = {'+': '-', '-': '+'}
df['aav_last_strand'] = df.apply(
    lambda x: inverse_strand[x['aav_last_strand']] if x['junction_locus'] != x['old_junction_locus'] else x[
        'aav_last_strand'], axis=1)


df = df[df.aav_matches>=30].copy()
df = df[df.start_aav==0].copy()

df.to_csv(outputbn,sep="\t", index=False)


