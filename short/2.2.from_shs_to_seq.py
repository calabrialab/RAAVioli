import pandas as pd
path = "FINAL/"
pool_id = "SPARK"
index_columns = 8
end_name = ".minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated."
file_shs = "ShsCount_"+pool_id+end_name+"removed1.tsv"
file_seq = "SeqCount_"+pool_id+end_name+"tsv"
df_seq_count = pd.read_csv(path+file_seq,sep="\t")
df_shs_count = pd.read_csv(path+file_shs,sep="\t")
key_cols = df_shs_count.columns[:index_columns]
pcr_cols = df_shs_count.columns[index_columns:]
df_seq_count['k']=df_seq_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_shs_count['k']=df_shs_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_seq_count.index=df_seq_count['k']
df_shs_count.index=df_shs_count['k']
for index, row in df_shs_count.iterrows():
    for c in pcr_cols:
        if not pd.isna(row[c]):
            df_shs_count.loc[[index],[c]]=df_seq_count.loc[[index],[c]]

df_shs_count = df_shs_count.drop('k',axis=1)
df_shs_count.to_csv(path+"SeqCount_from_"+file_shs, sep="\t",index=False)