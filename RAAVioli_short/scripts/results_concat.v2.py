import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input directory with results')
parser.add_argument('-o', '--output', help='output path+bn', default="")
parser.add_argument('-a', '--assofile', help='asso file', default="")
args = parser.parse_args()

input_path = Path(args.input)
outputbn = Path(args.output)
assofile = Path(args.assofile)

# Load association file
assodf = pd.read_csv(assofile, sep="\t")
tag_order = assodf['TagID'].tolist()

# Find input files
list_files = sorted(input_path.glob("*results.R1.tsv"))

res_df = []
for tag in tag_order:
    # Expected filename: "<tag>results.R1.tsv"
    expected_filename = f"{tag}results.R1.tsv"
    full_path = input_path / expected_filename

    if full_path.exists():
        df = pd.read_csv(
            full_path,
            sep="\t",
            )
        res_df.append(df)

if res_df:
    final_df = pd.concat(res_df, ignore_index=True)
else:
    final_df = pd.DataFrame()


final_df.to_csv(outputbn, sep="\t", index=False)
