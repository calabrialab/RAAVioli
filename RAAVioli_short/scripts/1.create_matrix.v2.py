import pandas as pd
import numpy as np

path = "hLSP1/"
df = pd.read_csv(path+"results.hLSP1.withrawread.modif.checked.tsv",sep="\t")

dfg = pd.DataFrame({'count': df.groupby(
    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","input_file","AddedField1", "aav_last_strand", "aav_alignments_start_end","n_aav_aln"]).size()}).reset_index()

dfg = pd.DataFrame({'count': dfg.groupby(
    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","AddedField1", "aav_last_strand", "aav_alignments_start_end","n_aav_aln"])['count'].sum()}).reset_index()
#dfg = pd.DataFrame({'count': dfg.groupby(
#    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","input_file"]).size()}).reset_index()
#dfg['junction_locus']=dfg['junction_locus'].map(lambda x: itr3_end-(x-itr5_start) if itr5_start<=x<=itr5_end else x)

results = {}
"""
key=chr_integration_strand_gap
results={chr1_3671051_+_-6:{junction:{input_file:count}}
"""
for index, row in dfg.iterrows():
    k=f"{row.target_chr}_{row.integration_locus}_{row.target_strand}_{row.gap}_{row.aav_alignments_start_end}_{row.n_aav_aln}_{row.aav_last_strand}"
    if k not in results.keys():
        results[k]={}
        results[k][row.junction_locus]={}
        results[k][row.junction_locus][row.AddedField1]=row['count']
    else:
        junction_present = False
        for other_junction in results[k]:
            if abs(other_junction-row.junction_locus)==0:
                junction_present = True
                if row.AddedField1 in results[k][other_junction].keys():
                    results[k][other_junction][row.AddedField1]+= row['count']
                else:
                    results[k][other_junction][row.AddedField1] = row['count']
                break
        if not junction_present:
            results[k][row.junction_locus]={}
            results[k][row.junction_locus][row.AddedField1] = row['count']



final_results = {}
for k in results.keys():
    for junction in results[k].keys():
        new_k = f'{k}_{junction}'
        final_results[new_k] = results[k][junction]
new_df = pd.DataFrame.from_dict(final_results)
r = new_df.T

r['k'] = r.index
r['IS_genomicID'] = r.apply(lambda x: "_".join(x['k'].split('_')[0:3]), axis=1)
r['gap'] = r.apply(lambda x: x['k'].split('_')[3], axis=1)
r['aav_alignments_start_end'] = r.apply(lambda x: x['k'].split('_')[4], axis=1)
r['n_aav_aln'] = r.apply(lambda x: x['k'].split('_')[5], axis=1)
r['aav_last_strand'] = r.apply(lambda x: x['k'].split('_')[6], axis=1)
r['junction'] = r.apply(lambda x: x['k'].split('_')[7], axis=1)

cols = r.columns[:-7]
old_df_col = ['IS_genomicID']
key_cols = ['IS_genomicID', 'gap', 'junction', 'aav_alignments_start_end', 'n_aav_aln', 'aav_last_strand']

key_cols.extend(cols)
old_df_col.extend(cols)
r_orig = r[old_df_col]
r_new = r[key_cols]
r_new.to_csv(path+"SeqCount_hLSP1.minfvalue30.gapabs120.aggregAF1.tsv",sep="\t",index=False)
df_seq_count = r_new.copy()




dfg = pd.DataFrame({'count': df.groupby(
    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","input_file","AddedField1", "aav_last_strand", "aav_alignments_start_end","n_aav_aln", "total_matches"
                    ]).size()}).reset_index()
dfg = pd.DataFrame({'count': dfg.groupby(
    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","input_file","AddedField1", "aav_last_strand", "aav_alignments_start_end","n_aav_aln"]).size()}).reset_index()

dfg = pd.DataFrame({'count': dfg.groupby(
    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","AddedField1", "aav_last_strand", "aav_alignments_start_end","n_aav_aln"])['count'].sum()}).reset_index()
#dfg = pd.DataFrame({'count': dfg.groupby(
#    ["target_chr", "integration_locus","target_strand", "junction_locus", "gap","input_file"]).size()}).reset_index()
#dfg['junction_locus']=dfg['junction_locus'].map(lambda x: itr3_end-(x-itr5_start) if itr5_start<=x<=itr5_end else x)

results = {}
"""
key=chr_integration_strand_gap
results={chr1_3671051_+_-6:{junction:{input_file:count}}
"""
for index, row in dfg.iterrows():
    k=f"{row.target_chr}_{row.integration_locus}_{row.target_strand}_{row.gap}_{row.aav_alignments_start_end}_{row.n_aav_aln}_{row.aav_last_strand}"
    if k not in results.keys():
        results[k]={}
        results[k][row.junction_locus]={}
        results[k][row.junction_locus][row.AddedField1]=row['count']
    else:
        junction_present = False
        for other_junction in results[k]:
            if abs(other_junction-row.junction_locus)==0:
                junction_present = True
                if row.AddedField1 in results[k][other_junction].keys():
                    results[k][other_junction][row.AddedField1]+= row['count']
                else:
                    results[k][other_junction][row.AddedField1] = row['count']
                break
        if not junction_present:
            results[k][row.junction_locus]={}
            results[k][row.junction_locus][row.AddedField1] = row['count']



final_results = {}
for k in results.keys():
    for junction in results[k].keys():
        new_k = f'{k}_{junction}'
        final_results[new_k] = results[k][junction]
new_df = pd.DataFrame.from_dict(final_results)
r = new_df.T

r['k'] = r.index
r['IS_genomicID'] = r.apply(lambda x: "_".join(x['k'].split('_')[0:3]), axis=1)
r['gap'] = r.apply(lambda x: x['k'].split('_')[3], axis=1)
r['aav_alignments_start_end'] = r.apply(lambda x: x['k'].split('_')[4], axis=1)
r['n_aav_aln'] = r.apply(lambda x: x['k'].split('_')[5], axis=1)
r['aav_last_strand'] = r.apply(lambda x: x['k'].split('_')[6], axis=1)
r['junction'] = r.apply(lambda x: x['k'].split('_')[7], axis=1)

cols = r.columns[:-7]
old_df_col = ['IS_genomicID']
key_cols = ['IS_genomicID', 'gap', 'junction', 'aav_alignments_start_end', 'n_aav_aln', 'aav_last_strand']

key_cols.extend(cols)
old_df_col.extend(cols)
r_orig = r[old_df_col]
r_new = r[key_cols]
r_new.to_csv(path+"ShsCount_hLSP1.minfvalue30.gapabs120.aggregAF1.tsv",sep="\t",index=False)

# redoing for both seq_count and shs{
r_new = r_new[~r_new['IS_genomicID'].str.contains("chrUn")].copy()
r_new = r_new[~r_new['IS_genomicID'].str.contains("_random")].copy()
r_new[['chr', 'integration', 'strand']] = r_new.IS_genomicID.str.split('_', expand=True)
r_new = r_new.astype({'junction': 'int', 'gap': 'int', 'integration': 'int', 'n_aav_aln': 'string'})
r_new['old_gap']=r_new['gap']
r_new['gap']=r_new.apply(lambda x: x['gap'] if x['old_gap']>0 else 0, axis = 1)
r_new['total'] = r_new.apply(lambda x: x['integration'] + x['junction'] if x['gap'] < 0 else
                                        (x['integration'] + x['junction'] + x['gap'] if x['aav_last_strand'] == '+' else
                                        (x['integration'] + x['junction'] - x['gap'] if x['strand'] == '+' else
                                         x['integration'] - x['junction'] + x['gap'])), axis=1)
r_new['gap']=r_new['old_gap']
r_new = r_new.drop('old_gap',axis=1)
results = {}
results_seq_count = {}
results_index = {}
max_diff = 8
max_group_diff = 20
xyz = list(set((i,j) for i,j in zip(r_new["chr"], r_new["n_aav_aln"])))
xyz = sorted(xyz, key = lambda kxyz: (kxyz[0], kxyz[1]))
for chr,aav_s_e in xyz:
    print(chr,aav_s_e)
    r_new_chr=r_new[(r_new.chr==chr) & (r_new.n_aav_aln==aav_s_e)]
    r_new_chr = r_new_chr.sort_values(by=['total'])
    total = None

    group = []
    group_index = []
    for index,row in r_new_chr.iterrows():
        if total==None:
            total = row['total']
            start_total = total
            group.append(row)
            group_index.append(index)
        else:
            if abs(total-row['total'])<=max_diff and abs(start_total-row['total'])<=max_group_diff:
                group.append(row)
                group_index.append(index)
                total = row['total']
            else:
                group_is =[]
                group_is_seq_count =[]
                df_group_is = pd.DataFrame(group)
                df_group_is = df_group_is.sort_values(by=['integration'])
                start_IS = df_group_is.iloc[0]['integration']
                for index_is,row_is in df_group_is.iterrows():
                    if abs(row_is['integration']-start_IS)<=20:
                        group_is.append(row_is)
                        group_is_seq_count.append(index_is)
                    else:
                        df_group = pd.DataFrame(group_is)
                        df_group_sorted = df_group.sort_values(by=cols.to_list(), ascending=False, na_position='last')
                        df_group_sorted = df_group_sorted.dropna(axis=1, how="all")
                        df_group_seq_count = df_seq_count.loc[group_is_seq_count]
                        for tag in cols:
                            gs_cols = df_group_sorted.columns
                            if tag in gs_cols:
                                df_group_sorted = df_group_sorted.sort_values(by=[tag], ascending=False, na_position='last')
                                key_of_gp = df_group_sorted.index[0]
                                if df_group_sorted.shape[0]>1 and not np.isnan(df_group_sorted.iloc[1][tag]) and df_group_sorted.iloc[0][tag]==df_group_sorted.iloc[1][tag]:
                                    df_group_seq_count_sorted = df_group_seq_count.sort_values(by=[tag],
                                                                                               ascending=False,
                                                                                               na_position='last')

                                    key_of_gp = df_group_seq_count_sorted.index[0]

                                df_group_sorted_tag = df_group_sorted[df_group_sorted[tag].notnull()]
                                tag_sum = df_group_sorted[tag].sum()
                                tag_sum_seq_count = df_group_seq_count[tag].sum()
                                int_site = df_group_sorted.loc[key_of_gp]['IS_genomicID']
                                gap = df_group_sorted.loc[key_of_gp]['gap']
                                junction = df_group_sorted.loc[key_of_gp]['junction']
                                aav_alignments_start_end = df_group_sorted.loc[key_of_gp]['n_aav_aln']
                                k = f'{int_site}_{gap}_{junction}_{aav_alignments_start_end}'
                                if k not in results.keys():
                                    results[k] = {}
                                    results_seq_count[k]={}

                                results[k][tag] = tag_sum
                                results_seq_count[k][tag]=tag_sum_seq_count
                                k_with_tag = f'{k}_{tag}'
                                results_index[k_with_tag]=list(df_group_sorted_tag.index)
                        group_is = []
                        group_is_seq_count = []
                        start_IS = row_is['integration']
                        group_is.append(row_is)
                        group_is_seq_count.append(index_is)
                df_group = pd.DataFrame(group_is)
                df_group_sorted = df_group.sort_values(by=cols.to_list(), ascending=False, na_position='last')
                df_group_sorted = df_group_sorted.dropna(axis=1, how="all")
                df_group_seq_count = df_seq_count.loc[group_is_seq_count]
                for tag in cols:
                    gs_cols = df_group_sorted.columns
                    if tag in gs_cols:
                        df_group_sorted = df_group_sorted.sort_values(by=[tag], ascending=False, na_position='last')
                        key_of_gp = df_group_sorted.index[0]
                        if df_group_sorted.shape[0]>1 and not np.isnan(df_group_sorted.iloc[1][tag]) and df_group_sorted.iloc[0][tag] == \
                                df_group_sorted.iloc[1][tag]:
                            df_group_seq_count_sorted = df_group_seq_count.sort_values(by=[tag],
                                                                                       ascending=False,
                                                                                       na_position='last')

                            key_of_gp = df_group_seq_count_sorted.index[0]
                        df_group_sorted_tag = df_group_sorted[df_group_sorted[tag].notnull()]
                        tag_sum = df_group_sorted[tag].sum()
                        tag_sum_seq_count = df_group_seq_count[tag].sum()
                        int_site = df_group_sorted.loc[key_of_gp]['IS_genomicID']
                        gap = df_group_sorted.loc[key_of_gp]['gap']
                        junction = df_group_sorted.loc[key_of_gp]['junction']
                        aav_alignments_start_end = df_group_sorted.loc[key_of_gp]['n_aav_aln']
                        k = f'{int_site}_{gap}_{junction}_{aav_alignments_start_end}'
                        if k not in results.keys():
                            results[k] = {}
                            results_seq_count[k] = {}

                        results[k][tag] = tag_sum
                        results_seq_count[k][tag] = tag_sum_seq_count
                        k_with_tag = f'{k}_{tag}'
                        results_index[k_with_tag] = list(df_group_sorted_tag.index)

                group=[]
                group_index=[]
                group.append(row)
                group_index.append(index)
                total = row['total']
                start_total = total
    # taking last value of each chr
    group_is = []
    group_is_seq_count = []
    df_group_is = pd.DataFrame(group)
    df_group_is = df_group_is.sort_values(by=['integration'])
    start_IS = df_group_is.iloc[0]['integration']
    for index_is, row_is in df_group_is.iterrows():
        if abs(row_is['integration'] - start_IS) <= 20:
            group_is.append(row_is)
            group_is_seq_count.append(index_is)
        else:
            df_group = pd.DataFrame(group_is)
            df_group_sorted = df_group.sort_values(by=cols.to_list(), ascending=False, na_position='last')
            df_group_sorted = df_group_sorted.dropna(axis=1, how="all")
            df_group_seq_count = df_seq_count.loc[group_is_seq_count]
            for tag in cols:
                gs_cols = df_group_sorted.columns
                if tag in gs_cols:
                    df_group_sorted = df_group_sorted.sort_values(by=[tag], ascending=False, na_position='last')
                    key_of_gp = df_group_sorted.index[0]
                    if df_group_sorted.shape[0]>1 and not np.isnan(df_group_sorted.iloc[1][tag]) and df_group_sorted.iloc[0][tag] == \
                            df_group_sorted.iloc[1][
                                tag]:
                        df_group_seq_count_sorted = df_group_seq_count.sort_values(by=[tag],
                                                                                   ascending=False,
                                                                                   na_position='last')

                        key_of_gp = df_group_seq_count_sorted.index[0]

                    df_group_sorted_tag = df_group_sorted[df_group_sorted[tag].notnull()]
                    tag_sum = df_group_sorted[tag].sum()
                    tag_sum_seq_count = df_group_seq_count[tag].sum()
                    int_site = df_group_sorted.loc[key_of_gp]['IS_genomicID']
                    gap = df_group_sorted.loc[key_of_gp]['gap']
                    junction = df_group_sorted.loc[key_of_gp]['junction']
                    aav_alignments_start_end = df_group_sorted.loc[key_of_gp]['n_aav_aln']
                    k = f'{int_site}_{gap}_{junction}_{aav_alignments_start_end}'
                    if k not in results.keys():
                        results[k] = {}
                        results_seq_count[k] = {}

                    results[k][tag] = tag_sum
                    results_seq_count[k][tag] = tag_sum_seq_count
                    k_with_tag = f'{k}_{tag}'
                    results_index[k_with_tag] = list(df_group_sorted_tag.index)
            group_is = []
            group_is_seq_count = []
            start_IS = row_is['integration']
            group_is.append(row_is)
            group_is_seq_count.append(index_is)
    df_group = pd.DataFrame(group_is)
    df_group_sorted = df_group.sort_values(by=cols.to_list(), ascending=False, na_position='last')
    df_group_sorted = df_group_sorted.dropna(axis=1, how="all")
    df_group_seq_count = df_seq_count.loc[group_is_seq_count]
    for tag in cols:
        gs_cols = df_group_sorted.columns
        if tag in gs_cols:
            df_group_sorted = df_group_sorted.sort_values(by=[tag], ascending=False, na_position='last')
            key_of_gp = df_group_sorted.index[0]
            if df_group_sorted.shape[0]>1 and not np.isnan(df_group_sorted.iloc[1][tag]) and df_group_sorted.iloc[0][tag] == df_group_sorted.iloc[1][
                tag]:
                df_group_seq_count_sorted = df_group_seq_count.sort_values(by=[tag],
                                                                           ascending=False,
                                                                           na_position='last')

                key_of_gp = df_group_seq_count_sorted.index[0]
            df_group_sorted_tag = df_group_sorted[df_group_sorted[tag].notnull()]
            tag_sum = df_group_sorted[tag].sum()
            tag_sum_seq_count = df_group_seq_count[tag].sum()
            int_site = df_group_sorted.loc[key_of_gp]['IS_genomicID']
            gap = df_group_sorted.loc[key_of_gp]['gap']
            junction = df_group_sorted.loc[key_of_gp]['junction']
            aav_alignments_start_end = df_group_sorted.loc[key_of_gp]['n_aav_aln']
            k = f'{int_site}_{gap}_{junction}_{aav_alignments_start_end}'
            if k not in results.keys():
                results[k] = {}
                results_seq_count[k] = {}

            results[k][tag] = tag_sum
            results_seq_count[k][tag] = tag_sum_seq_count
            k_with_tag = f'{k}_{tag}'
            results_index[k_with_tag] = list(df_group_sorted_tag.index)

collapsed_df = pd.DataFrame.from_dict(results)
r = collapsed_df.T

r['k']=r.index
r['integration_chr']=r.apply(lambda x: x['k'].split('_')[0], axis=1)
r['integration_locus']=r.apply(lambda x: x['k'].split('_')[1], axis=1)
r['integration_strand']=r.apply(lambda x: x['k'].split('_')[2], axis=1)
r['gap']=r.apply(lambda x: x['k'].split('_')[3], axis=1)
r['junction']=r.apply(lambda x: x['k'].split('_')[4], axis=1)
r['aav_alignments_start_end']=r.apply(lambda x: x['k'].split('_')[5], axis=1)

r=r.astype({"integration_locus": int, "gap": int})



r['IS_genomicID']=r.apply(lambda x: "_".join([x['integration_chr'],str(x['integration_locus']),x['integration_strand']]), axis=1)

cols = r.columns[:-8]
old_df_col = ['IS_genomicID']
key_cols = ['IS_genomicID','gap','junction','aav_alignments_start_end']
key_cols.extend(cols)
old_df_col.extend(cols)
cdf = r[key_cols].copy()

cdf.to_csv(path+"ShsCount_hLSP1.minfvalue30.gapabs120.aggregAF1.window8_20.tsv",sep="\t",index=False)
cols_sum = cdf.columns[4:]
cdf['sum'] = cdf[cols_sum].sum(axis=1)
cdf_cleaned = cdf[cdf['sum']>1].copy()
cdf_cleaned = cdf_cleaned.drop('sum',axis=1)
cdf_cleaned.to_csv(path+"ShsCount_hLSP1.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.tsv",sep="\t",index=False)



collapsed_df = pd.DataFrame.from_dict(results_seq_count)
r = collapsed_df.T

r['k']=r.index
r['integration_chr']=r.apply(lambda x: x['k'].split('_')[0], axis=1)
r['integration_locus']=r.apply(lambda x: x['k'].split('_')[1], axis=1)
r['integration_strand']=r.apply(lambda x: x['k'].split('_')[2], axis=1)
r['gap']=r.apply(lambda x: x['k'].split('_')[3], axis=1)
r['junction']=r.apply(lambda x: x['k'].split('_')[4], axis=1)
r['aav_alignments_start_end']=r.apply(lambda x: x['k'].split('_')[5], axis=1)

r=r.astype({"integration_locus": int, "gap": int})



r['IS_genomicID']=r.apply(lambda x: "_".join([x['integration_chr'],str(x['integration_locus']),x['integration_strand']]), axis=1)

cols = r.columns[:-8]
old_df_col = ['IS_genomicID']
key_cols = ['IS_genomicID','gap','junction','aav_alignments_start_end']
key_cols.extend(cols)
old_df_col.extend(cols)
cdf = r[key_cols].copy()

cdf.to_csv(path+"SeqCount_hLSP1.minfvalue30.gapabs120.aggregAF1.window8_20.tsv",sep="\t",index=False)
cols_sum = cdf.columns[4:]
cdf['sum'] = cdf[cols_sum].sum(axis=1)
cdf_cleaned = cdf[cdf['sum']>1].copy()
cdf_cleaned = cdf_cleaned.drop('sum',axis=1)
cdf_cleaned.to_csv(path+"SeqCount_hLSP1.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.tsv",sep="\t",index=False)