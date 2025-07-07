#!/usr/bin/env python3

import os
import pandas as pd
from datetime import datetime
import sys

# ---- Usage: python append_stats.py <output_tsv_path> ----

# Expected environment variables
fields = [
    "RUN_ID", "RUN_NAME", "DISEASE", "PATIENT", "POOL", "TAG",
    "LTR_ID", "LC_ID", "RAW_READS", "QUALITY_PASSED", "PHIX_MAPPING",
    "PLASMID_MAPPED_BYPOOL", "RAW_NO_PLASMID", "BARCODE_MUX",
    "LTR_IDENTIFIED", "TRIMMING_FINAL_LTRLC", "LV_MAPPED",
    "BWA_MAPPED_OVERALL", "ISS_MAPPED_OVERALL", "ISS_MAPPED_PP"
]

# Output path passed from bash
out_path = sys.argv[1]

# Create row dict with timestamp
now = datetime.now()
row = {
    "DATE": now.strftime("%Y-%m-%d"),
    "TIME": now.strftime("%H:%M")
}

# Add each field from env
for field in fields:
    row[field] = os.environ.get(field, "NA")

# Load or create TSV
if os.path.exists(out_path):
    df = pd.read_csv(out_path, sep="\t")
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
else:
    df = pd.DataFrame([row])

# Save back to TSV
df.to_csv(out_path, sep="\t", index=False)

