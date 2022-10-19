import pandas as pd
import numpy as np
import pandas as pd
import scipy

def getOverlap_bed(annot,row_is,start_col,end_col,k=5):
    nearest_feat = "None"
    nearest_dist_start = 10000000000
    for index, row in annot.iterrows():
        if row_is[start_col]>=row['start']-k and row_is[end_col]<=row['end']+k:
            return row_is[end_col]-row_is[start_col], "Full", row['name']
        if row_is[start_col]<=row['start']-k and ( row['start']<=row_is[end_col]<=row['end']):
            return row['start']-row_is[start_col], "Partial", row['name']
        if (row['start']<=row_is[start_col]<=row['end']) and row_is[end_col]>=row['end']+k:
            return row_is[end_col]-row['end'], "Partial", row['name']
        if abs(row_is[start_col]-row['start'])<nearest_dist_start:
            nearest_dist_start = abs(row_is[start_col]-row['start'])
            nearest_feat = row['name']

    return nearest_dist_start, "None", nearest_feat

path = "FINAL/post_mod/"
file = "EXTENDED.MODIFIED.Shs_from_SeqCount_from_ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1.nocoll.tsv"
df = pd.read_csv(path+file,sep="\t")
cols_k = ['chr', 'integration_locus', 'strand', 'gap', 'junction','aav_alignments_start_end', 'ID']
df['new_k'] = df.apply(lambda row: "_".join([str(row[i]) for i in cols_k]),axis=1)


df_gaa = df[df.GeneName=="GAA"].copy()
bed_gaa = pd.read_csv(path+"tb_GAA_table_browser_Exons.bed",sep="\t")
df_gaa[['r1_bases_overlap','r1_overlap_type','r1_feature_overlap']]=df_gaa.apply(lambda x: getOverlap_bed(bed_gaa,x,"target_start","target_end"),axis=1,result_type="expand")
df_gaa[['r2_bases_overlap','r2_overlap_type','r2_feature_overlap']]=df_gaa.apply(lambda x: getOverlap_bed(bed_gaa,x,"target_start_r2","target_end_r2"),axis=1,result_type="expand")
df_gaa_tonotexclude = df_gaa[~((df_gaa['r1_overlap_type']=="Full") & (df_gaa['r2_overlap_type']=="Full"))].copy()
df_gaa_toexclude = df_gaa[((df_gaa['r1_overlap_type']=="Full") & (df_gaa['r2_overlap_type']=="Full"))].copy()
#df_gaa_tonotexclude.to_csv(path+"gaa_is.tsv",sep="\t",index_label="original_index")


df_serpina = df[df.GeneName=="SERPINA1"].copy()
bed_serpina = pd.read_csv(path+"tb_promoter_serpina_table_browser.bed",sep="\t")
df_serpina[['r1_bases_overlap','r1_overlap_type','r1_feature_overlap']]=df_serpina.apply(lambda x: getOverlap_bed(bed_serpina,x,"target_start","target_end"),axis=1,result_type="expand")
df_serpina[['r2_bases_overlap','r2_overlap_type','r2_feature_overlap']]=df_serpina.apply(lambda x: getOverlap_bed(bed_serpina,x,"target_start_r2","target_end_r2"),axis=1,result_type="expand")
df_serpina_tonotexclude = df_serpina[~((df_serpina['r1_overlap_type']=="Full") & (df_serpina['r2_overlap_type']=="Full"))].copy()
df_serpina_toexclude = df_serpina[((df_serpina['r1_overlap_type']=="Full") & (df_serpina['r2_overlap_type']=="Full"))].copy()

#df_serpina_tonotexclude.to_csv(path+"SERPINA_is.tsv",sep="\t",index_label="original_index")

df_apoc = df[df.GeneName=="APOC1"].copy()
bed_apoc = pd.read_csv(path+"tb_enhancer_apoc_table_browser.bed",sep="\t")
df_apoc[['r1_bases_overlap','r1_overlap_type','r1_feature_overlap']]=df_apoc.apply(lambda x: getOverlap_bed(bed_apoc,x,"target_start","target_end"),axis=1,result_type="expand")
df_apoc[['r2_bases_overlap','r2_overlap_type','r2_feature_overlap']]=df_apoc.apply(lambda x: getOverlap_bed(bed_apoc,x,"target_start_r2","target_end_r2"),axis=1,result_type="expand")
df_apoc_tonotexclude = df_apoc[~((df_apoc['r1_overlap_type']=="Full") & (df_apoc['r2_overlap_type']=="Full"))].copy()
df_apoc_toexclude = df_apoc[((df_apoc['r1_overlap_type']=="Full") & (df_apoc['r2_overlap_type']=="Full"))].copy()
df_apoc_tonotexclude['target_read'] = df_apoc_tonotexclude.apply(lambda x: x['raw_read'][x['end_aav_position']:],axis=1)
with open(path+"apoc_seqs.fa","w") as f:
    for i,r in df_apoc_tonotexclude.iterrows():
        f.write(f">{r['name']}\n{r['target_read']}\n")

#df_apoc_tonotexclude.to_csv(path+"APOC_is.tsv",sep="\t",index_label="original_index")
df_hbb = df[df.GeneName=="HBB"].copy()
bed_hbb = pd.read_csv(path+"tb_Bglobin_intron_table_browser.bed",sep="\t")
df_hbb[['r1_bases_overlap','r1_overlap_type','r1_feature_overlap']]=df_hbb.apply(lambda x: getOverlap_bed(bed_hbb,x,"target_start","target_end"),axis=1,result_type="expand")
df_hbb[['r2_bases_overlap','r2_overlap_type','r2_feature_overlap']]=df_hbb.apply(lambda x: getOverlap_bed(bed_hbb,x,"target_start_r2","target_end_r2"),axis=1,result_type="expand")
df_hbb_tonotexclude = df_hbb[~((df_hbb['r1_overlap_type']=="Full") & (df_hbb['r2_overlap_type']=="Full"))].copy()
df_hbb_toexclude = df_hbb[((df_hbb['r1_overlap_type']=="Full") & (df_hbb['r2_overlap_type']=="Full"))].copy()

df_to_exclude = pd.concat([df_hbb_toexclude, df_apoc_toexclude, df_serpina_toexclude, df_gaa_toexclude])
df_to_exclude.to_csv(path+"EXCLUDEDE.target_genes.tsv",sep="\t", index=False)
df_ssc = df[~(df['new_k'].isin(df_hbb_toexclude['new_k'])) & ~(df['new_k'].isin(df_serpina_toexclude['new_k'])) &
   ~(df['new_k'].isin(df_apoc_toexclude['new_k']))  & ~(df['new_k'].isin(df_gaa_toexclude['new_k']))].copy()
df_ssc.to_csv(path+"EXTENDED.SPARK.removed_targetgenes.tsv",sep="\t", index=False)



idx_cols = list(df_ssc.columns[:8])
#idx_cols.extend(["aav_alignments_start_end_y"])
newdf = df_ssc.pivot(index=idx_cols,columns='ID', values='ssc').reset_index()
newdf.to_csv(path+"SPARK.removed_targetgenes.tsv",sep="\t", index=False)