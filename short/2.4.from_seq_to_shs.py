import pandas as pd
path = "FINAL/"
pool_id = "SPARK"
index_columns = 8
end_name = ".minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1."
file_shs = "ShsCount_"+pool_id+end_name+"tsv"
file_seq = "SeqCount_from_ShsCount_"+pool_id+end_name+"nocoll.tsv"
df_seq_count = pd.read_csv(path+file_seq,sep="\t")

df_shs_count = pd.read_csv(path+file_shs,sep="\t")
key_cols = df_seq_count.columns[:index_columns]
pcr_cols = df_seq_count.columns[index_columns:]
df_seq_count['k']=df_seq_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_shs_count['k']=df_shs_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_seq_count.index=df_seq_count['k']
df_shs_count.index=df_shs_count['k']
for index, row in df_seq_count.iterrows():
    for c in pcr_cols:
        if not pd.isna(row[c]):
            df_seq_count.loc[[index],[c]]=df_shs_count.loc[[index],[c]]

df_seq_count = df_seq_count.drop('k',axis=1)
df_seq_count = df_seq_count.sort_values(['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction'])
df_seq_count = df_seq_count[abs(df_seq_count.gap)<=30].copy()
df_seq_count.to_csv(path+"Shs_from_"+file_seq, sep="\t",index=False)
df = df_seq_count.copy()
infcols = list(df.columns[:8])
pccols = list(df.columns[8:])
pccols.sort()
infcols.extend(pccols)
df2 =  df[infcols]

df2.to_csv(path+"sorteccols.Shs_from_"+file_seq, sep="\t",index=False)

idx_cols = list(df2.columns[:8])
melted_is  =df2.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
dfvc=pd.DataFrame(melted_is['integration_locus'].value_counts())
dfvc.columns=['count_is']
dfvc['integration_locus'] = dfvc.index
mergedf=melted_is.merge(dfvc, on='integration_locus')
mergedf.to_csv(path+"meltedf.withcounts.tsv",sep="\t",index=False)

"""
import pandas as pd
path = "fixedjun/"
pool_id = "SPARK"
index_columns = 8
end_name = ".minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1."
file_shs = "ShsCount_"+pool_id+end_name+"tsv"
file_seq = "SeqCount_from_ShsCount_"+pool_id+end_name+"nocoll.tsv"
df_seq_count = pd.read_csv(path+file_seq,sep="\t")
df_seq_count =df_seq_count.astype({"chr":str})

df_shs_count = pd.read_csv(path+file_shs,sep="\t")
df_shs_count =df_shs_count.astype({"chr":str})

key_cols = df_seq_count.columns[:index_columns]
pcr_cols = df_seq_count.columns[index_columns:]
df_seq_count['k']=df_seq_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_shs_count['k']=df_shs_count.apply(lambda x: "_".join([str(x[c]) for c in key_cols]),axis=1)
df_seq_count.index=df_seq_count['k']
df_shs_count.index=df_shs_count['k']
for index, row in df_seq_count.iterrows():
    for c in pcr_cols:
        if not pd.isna(row[c]):
            df_seq_count.loc[[index],[c]]=df_shs_count.loc[[index],[c]]

df_seq_count = df_seq_count.drop('k',axis=1)
df_seq_count = df_seq_count.sort_values(['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction'])
ord_cols = ['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction','aav_alignments_start_end']
ord_cols.extend(sorted(pcr_cols))
df_seq_count = df_seq_count[abs(df_seq_count.gap)<=30].copy()
df_seq_count = df_seq_count[ord_cols]
df_seq_count.to_csv(path+"Shs_from_"+file_seq, sep="\t",index=False)





df = df_seq_count.copy()
infcols = list(df.columns[:8])
pccols = list(df.columns[8:])
pccols.sort()
infcols.extend(pccols)
df2 =  df[infcols]

df2.to_csv(path+"sorteccols.Shs_from_"+file_seq, sep="\t",index=False)

idx_cols = list(df2.columns[:8])
melted_is  =df2.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
dfvc=pd.DataFrame(melted_is['integration_locus'].value_counts())
dfvc.columns=['count_is']
dfvc['integration_locus'] = dfvc.index
mergedf=melted_is.merge(dfvc, on='integration_locus')
mergedf.to_csv(path+"meltedf.withcounts.tsv",sep="\t",index=False)"""