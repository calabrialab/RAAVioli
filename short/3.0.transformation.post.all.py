# se punto IS<=3 abs se stesso topo e stesso sistema, ignora GAP e Junction e portale al 5' (SISTEMA P)
import pandas as pd
import os
path = "FINAL/"
file = "Shs_from_SeqCount_from_ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1.nocoll.tsv"
df = pd.read_csv(path+file, sep="\t")
#df = df.drop(columns=["RANDY-LIVER-1000-I"]).copy()
df['count'] = df[df.columns[8:]].count(axis=1)
df = df[df['count']!=0].copy()
df = df.drop(columns=["count"]).copy()
idx_cols = list(df.columns[:8])
melted_is = df.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
df_sorted = melted_is.sort_values(['chr','integration_locus','gap','junction','aav_alignments_start_end','ID']).reset_index(drop=True)
df_sorted['SubjectID'] =df_sorted['ID'].apply(lambda x: x.split("-")[0])
df_sorted['System'] =df_sorted['ID'].apply(lambda x: x.split("-")[-1])
"""
condivisione giunzione IS in topi differenti	teniamo solo quella con ssc più alto
stesso topo 2 ITR	siccome sistema è il P e oligo expo è sulla porzione al 5' se non ci sono riarrangiamenti va mappata sulla 5' ITR sommando gli ssc (somma anche se ci sono I riarrangiamenti)
stesso topo 2 ITR	siccome sistema è ha oligo expo è sulla porzione al 5' se non ci sono riarrangiamenti va mappata sulla 5' ITR sommando gli ssc (somma anche se ci sono I riarrangiamenti)

stesso topo diversi riarrangimenti	le IS si sommano e si tengono le features di quella con il riarrangiamento
stesso topo diversi sistemi	NON DOBBIAMO FAR NIENTE PERCHè è UNA SOLA IS
"""

itr3_end = 5352
itr3_start = 5212

for i in range(df_sorted.shape[0]-1):
    if df_sorted.loc[i]['GeneName'] not in ['APOC1','GAA','SERPINA1','HBB']:
        junc_diff = abs(df_sorted.loc[i]['junction'] - df_sorted.loc[i + 1]['junction'])
        gap_diff = abs(df_sorted.loc[i]['gap'] - df_sorted.loc[i + 1]['gap'])
        is_diff = abs(df_sorted.loc[i]['integration_locus']-df_sorted.loc[i+1]['integration_locus'])
        if not (junc_diff==0 and gap_diff==0 and is_diff==0):
            if is_diff<=3 and df_sorted.loc[i]['chr']==df_sorted.loc[i+1]['chr']:
                if df_sorted.loc[i]['ID'] == df_sorted.loc[i+1]['ID'] and df_sorted.loc[i]['System'] in (['P','F','G','H'])\
                        and df_sorted.loc[i + 1]['aav_alignments_start_end']==df_sorted.loc[i]['aav_alignments_start_end']:
                    if itr3_start<=df_sorted.loc[i+1]['junction']<=itr3_end:
                        df_sorted.at[i+1,'ssc'] += df_sorted.loc[i]['ssc']
                        df_sorted.at[i,'ssc'] = 0
                    elif itr3_start<=df_sorted.loc[i]['junction']<=itr3_end:
                        df_sorted.at[i,'ssc'] += df_sorted.loc[i+1]['ssc']
                        df_sorted.at[i+1, 'ssc'] = 0
                elif df_sorted.loc[i]['ID'] == df_sorted.loc[i+1]['ID']:
                    print("rear ", i)
                    if df_sorted.loc[i + 1]['aav_alignments_start_end']>df_sorted.loc[i]['aav_alignments_start_end']:
                        df_sorted.at[i + 1, 'ssc'] += df_sorted.loc[i]['ssc']
                        df_sorted.at[i, 'ssc'] = 0
                    elif df_sorted.loc[i]['aav_alignments_start_end']>df_sorted.loc[i+1]['aav_alignments_start_end']:
                        df_sorted.at[i, 'ssc'] += df_sorted.loc[i+1]['ssc']
                        df_sorted.at[i + 1, 'ssc'] = 0
                    else:
                        print("wffID",i)
                else:
                    junc_diff = abs(df_sorted.loc[i]['junction'] - df_sorted.loc[i + 1]['junction'])
                    gap_diff = abs(df_sorted.loc[i]['gap'] - df_sorted.loc[i + 1]['gap'])
                    if junc_diff<=3 and gap_diff<=3:
                        if df_sorted.loc[i]['ssc'] > df_sorted.loc[i+1]['ssc']:
                            df_sorted.at[i+1, 'ssc'] = 0
                        else:
                            df_sorted.at[i, 'ssc'] = 0
                    else:
                        print("Not Same ID: ",i, junc_diff, gap_diff)

    #max_diff = max([df_sorted.loc[i][col] - df_sorted.loc[i+1][col] for col in ['gap','junction'] ])

df_sorted_new = df_sorted[df_sorted['ssc']!=0].copy()
df_sorted_new.drop(['SubjectID','System'],axis=1, inplace=True)
idx_cols = df_sorted_new.columns[:8]
newdf = df_sorted_new.pivot(index=idx_cols,columns='ID', values='ssc').reset_index()

path = "FINAL/post_mod2/"
os.makedirs(path, exist_ok=True)

newdf.to_csv(path+"MODIFIED."+file,sep="\t", index=False)
df2 = newdf
idx_cols = list(df2.columns[:8])
melted_is  =df2.melt(id_vars=idx_cols,
            var_name="ID",
            value_name="ssc").dropna().reset_index(drop=True)
dfvc=pd.DataFrame(melted_is['integration_locus'].value_counts())
dfvc.columns=['count_is']
dfvc['integration_locus'] = dfvc.index
mergedf=melted_is.merge(dfvc, on='integration_locus')
mergedf.sort_values(idx_cols,inplace=True)
mergedf.to_csv(path+"meltedf.MODIFIED.withcounts.tsv",sep="\t",index=False)

