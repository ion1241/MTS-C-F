#!/usr/bin/env python3
# This script extracts ORF IDs from a filtered file and saves them to a Prok_IDs.txt file.
#This is used to create a list of prokaryotic ORFs for further analysis.

import pandas as pd

# Read the filtered file into a pandas DataFrame
filtered_df = pd.read_csv("ProkORF_Info.txt", sep='\t')

# Extract the ORF ID column
orf_ids = filtered_df['ORF ID']

# Save the ORF IDs to a new .txt file
orf_ids.to_csv("Prok_IDs.txt", index=False, header=False)
