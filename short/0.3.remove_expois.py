import pandas as pd
import os

path = "FINAL/"
df = pd.read_csv(path + "results.SPARK.withrawread.modif.tsv", sep="\t")

path_flex = path+"flexOut/"
list_out = [x for x in os.listdir(path_flex) if "flexbarOut" in x]
dict_reads = {}
for file in list_out:
    reads_trimmed = open(path_flex + file, "r").read().splitlines()
    name_r = reads_trimmed[0]
    list_seqs = []
    for i in range(1, len(reads_trimmed)):
        r = reads_trimmed[i]
        if r.startswith(">"):
            dict_reads[name_r] = "".join(list_seqs)
            list_seqs = []
            name_r = r
        else:
            list_seqs.append(r)
    dict_reads[name_r] = "".join(list_seqs)

df_reads = pd.DataFrame.from_dict(dict_reads, orient="index", columns=["trimmed_read"])
df_reads["name"] = df_reads.index
df_reads['name'] = df_reads.name.apply(lambda x: x.replace(">", ""))
df_reads['trimmed_length'] = df_reads.trimmed_read.apply(lambda x: len(x))
df_reads.reset_index(inplace=True, drop=True)

mergedf = df.merge(df_reads, on="name", how="left")
mergedf['diff_length'] = mergedf['read_length'] - mergedf['trimmed_length']
mergedf['end_aav_position'] = mergedf['alignment_on_query'].apply(lambda x: int(x.split(";")[-2].split(",")[-2]))
mergedf['start_targe_position'] = mergedf['alignment_on_query'].apply(lambda x: int(x.split(";")[-1].split(",")[-3]))
mergedf['expo_start'] = mergedf['diff_length'] - 20
mergedf['distance_expo_aavend'] = mergedf['expo_start'] - mergedf['end_aav_position']
mergedf['dist_expoaav_gap'] = mergedf['distance_expo_aavend'] - mergedf['gap']
mergedf['dist_expoend_targ_start'] = mergedf['diff_length'] - mergedf['start_targe_position']
mergedf['last_aav_length'] = mergedf['aav_alignments_start_end'].apply(
    lambda x: int(x.split(";")[-1].split(",")[-2]) - int(x.split(";")[-1].split(",")[-3]))

mergedf_towrite = mergedf[~((mergedf['distance_expo_aavend'] > -20) | (
            (-30 <= mergedf['distance_expo_aavend']) & (mergedf['distance_expo_aavend'] <= -20) & (
                mergedf['last_aav_length'] <= 30)))].copy()
mergedf_towrite.to_csv(path + "results.SPARK.withrawread.modif.checked.tsv", sep="\t", index=False)
