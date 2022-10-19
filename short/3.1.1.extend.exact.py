import pandas as pd


path_all = "FINAL/"

file_all = "results.SPARK.withrawread.modif.checked.tsv"

path_is = "FINAL/post_mod/"
file_is = "MODIFIED.Shs_from_SeqCount_from_ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1.nocoll.tsv"

df_is = pd.read_csv(path_is+file_is,sep="\t")
df_all = pd.read_csv(path_all+file_all, sep="\t")
idx_cols = list(df_is.columns[:8])
melted_is  =df_is.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
cols_k = [x for x in melted_is.columns[:-1] if "Gene" not in x]
melted_is['chr'] = melted_is['chr'].apply(lambda x: "chr"+str(x))
melted_is['k']=melted_is.apply(lambda x: "_".join(str(x[i]) for i in cols_k), axis=1)
cols_all_k = ['target_chr','integration_locus','target_strand','gap',
                  'junction_locus','n_aav_aln','AddedField1']
df_all['k'] = df_all.apply(lambda x: "_".join(str(x[i]) for i in cols_all_k), axis=1)
melted_is2 = melted_is.merge(df_all,on=['k'],suffixes=('', '_y'))
final_df = melted_is2.sort_values('total_matches', ascending=False).drop_duplicates('k').sort_index()
final_df = final_df.sort_values(['chr','integration_locus','gap','junction'])
final_df.to_csv(path_is+"EXTENDED."+file_is,sep="\t", index=False)
