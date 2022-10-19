import pandas as pd

path_is = "FINAL/post_mod/"
file_is = "SPARK.removed_targetgenes.tsv"
path_all = "FINAL/"

file_all = "results.SPARK.withrawread.modif.checked.tsv"

df_is = pd.read_csv (path_is+file_is, sep="\t")
df_all = pd.read_csv(path_all+file_all, sep="\t")
idx_cols = list(df_is.columns[:8])
melted_is  =df_is.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
cols_k = [x for x in melted_is.columns[:-1] if "Gene" not in x]
#melted_is['chr'] = melted_is['chr'].apply(lambda x: "chr"+str(x))
df_extend = pd.DataFrame(columns=list(df_all.columns).extend(['new_k']))
for index,row in melted_is.iterrows():
    df_sub = df_all[(df_all['AddedField1']==row['ID']) & (df_all['target_chr']==row['chr']) & (df_all['n_aav_aln']==row['aav_alignments_start_end']) &
                    (abs(row['integration_locus']-df_all['integration_locus'])<=10) &
                    (abs(row['junction']-df_all['junction_locus'])<=10) & (abs(row['gap']-df_all['gap'])<=10)].copy()
    new_k = "_".join([str(row[i]) for i in cols_k])
    df_sub['new_k'] = new_k
    df_extend = pd.concat([df_extend,df_sub])
df_extend.to_csv(path_is+"EXTENDED.WITHERRMARG."+file_is,sep="\t",index=False)
df_extend.merge(df_extend.groupby("new_k")['CompleteAmplificationID'].nunique(), on ="new_k")

file_extended = "EXTENDED.SPARK.removed_targetgenes.tsv"
df_precise_extend = pd.read_csv(path_is+file_extended,sep="\t")
df_precise_extend['new_k'] = df_precise_extend.apply(lambda row: "_".join([str(row[i]) for i in cols_k]),axis=1)
df_merged = df_precise_extend.merge(df_extend.groupby("new_k")['CompleteAmplificationID'].nunique(), on ="new_k")
df_merged.columns = [x if x!="CompleteAmplificationID_y" else "PCR_count" for x in df_merged.columns]
df_merged.to_csv(path_is+"PCR_COUNT."+file_extended,sep="\t", index=False)
