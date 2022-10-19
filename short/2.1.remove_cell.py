import pandas as pd
import numpy as np
path = "FINAL/"
file="ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.tsv"
df = pd.read_csv(path+file,sep="\t")
orig_columns =list(df.columns[:8])
pcc=df.columns[8:]
df[pcc]=df[pcc].replace(1,np.nan)
df['sum']=df[pcc].sum(axis=1)
df = df[df['sum']!=0].copy()
df = df.drop('sum',axis=1)
cols = df[pcc].columns[df[pcc].sum()>0]
orig_columns.extend(cols)
df_final = df[orig_columns]
file = file.replace(".tsv","")
df_final = df_final.sort_values(['chr','integration_locus','strand','GeneName',
                                         'GeneStrand','gap','junction'])
df_final.to_csv(path+file+".removed1.tsv",sep="\t",index=False)
