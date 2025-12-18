import argparse
import pandas as pd
from os import listdir
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input directory with results')
parser.add_argument('-o', '--output', help='output path+bn', default="")
parser.add_argument('-a', '--assofile', help='asso file', default="")
args = parser.parse_args()
input_path = args.input
outputbn = args.output
assofile = args.assofile

assodf = pd.read_csv(assofile,sep="\t")
list_files = [x for x in listdir(input_path) if x.endswith("results.R1.tsv")]
assodf_tagid = assodf['TagID'].to_list()
res_df = []
for infile in list_files:
    tag = Path(infile).name.replace(".results.R1.tsv", "")
    if tag in assodf_tagid:
        data = pd.read_csv(f'{input_path}/{infile}',sep="\t")
        res_df.append(data)


appended_data = pd.concat(res_df)

appended_data.to_csv(outputbn,sep="\t",index=False)



