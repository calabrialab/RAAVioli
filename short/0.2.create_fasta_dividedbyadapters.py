import pandas as pd
path = "FINAL/"
file = "results.SPARK.withrawread.modif.tsv"
df = pd.read_csv(path+file, sep="\t")

df['pcr_system'] = df.AddedField1.apply(lambda x: x.split("-")[-1])
df_systems = pd.read_csv("system_sequences.tsv",sep="\t")
df=df.merge(df_systems, left_on="pcr_system", right_on="System")
df['Expo'] = df['Expo_y']
for expo in df.Expo.unique():
    dfe = df[df['Expo']==expo].copy()
    with open(f'check/{expo}_reads.fa','w+') as f:
        for index, row in dfe.iterrows():
            f.write(f'>{row["name"]}\n')
            f.write(row['raw_read']+"\n")

def reverseComplement(seq):
    dict_nucleo = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join([dict_nucleo[i] for i in seq[::-1]])

for index, row in df_systems.iterrows():
    with open(f'{path}{row["Expo"]}.fa','w+') as f:
        f.write(">"+row['Expo']+"\n")
        f.write(row['Seq_expo']+"\n")
        f.write(">" + row['Expo'] + "_revcomp\n")
        f.write(reverseComplement(row['Seq_expo']) + "\n")
