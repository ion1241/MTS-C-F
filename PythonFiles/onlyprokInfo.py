#!/usr/bin/env python3
# This script extracts  prokaryotic ORFs from a All ORFs file and saves them to a ProkORF_Info.txt file.

import pandas as pd

# Read the file into a pandas DataFrame
df = pd.read_csv("ORFInfo_Clean.txt", sep='\t')

# Filter the DataFrame
filtered_df = df[(df['Tax'].str.contains('k_B|k_A', na=False)) & (df['KEGG'].notna()) & (df['KEGG'] != '')]

# Rewrite the file with the filtered DataFrame
filtered_df.to_csv("ProkORF_Info.txt", sep='\t', index=False)
