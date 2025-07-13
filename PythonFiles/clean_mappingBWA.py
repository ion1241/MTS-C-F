# This script is used to clean the output files from the mapped tables.
# the output files are in the format of "filename.txt" and the script will remove the columns that are not needed, modify them apropiately and
# and rename the second column with the start of the filename. The cleaned files will be saved as "filename_clean.txt".

import pandas as pd
import sys

def clean_and_save_file(filename):
    # Read the file into a pandas DataFrame
    df = pd.read_csv(filename, delim_whitespace=True, header=0)
   
    # Drop the unusfull columns
    dfmod=df.drop(df.columns[[1,3,4,5,8,9,10]], axis=1)

    #Sum columns 2 and 3 into a new column
    dfmod['sum'] = dfmod['Minus_reads'] + dfmod['Plus_reads']

    #remove rows with sum=0 from dfmod
    dfmod = dfmod[dfmod['sum'] != 0]

    #divide sum by Length to get coverage in Dfmod into a new column and remove the rest
    dfmod['coverage'] = dfmod['sum'] / dfmod['Length']

    dfmod2=dfmod.drop(dfmod.columns[[1,2,3,4]], axis=1)
    
    # Rename the second column with the start of the filename
    new_filename = filename.split('_')[0] + "_clean.txt"

    dfmod2.rename({'coverage': new_filename.split('_')[0]}, axis=1, inplace=True, errors='raise')
    
    # Save the cleaned DataFrame to a new file
    dfmod2.to_csv(new_filename, sep='\t', index=False)

if __name__ == "__main__":
    # Check if filename is provided as an argument
    if len(sys.argv) < 1:
        print("Usage: python script.py <filename1> <filename2> ...")
        sys.exit(1)
    
    # Retrieve the filename from the command line argument
    for filename in sys.argv[1:]:
        clean_and_save_file(filename)
    
        # Call the function to clean and save the file
        clean_and_save_file(filename)

