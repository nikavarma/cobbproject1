import pandas as pd
from Bio import SeqIO
import re
import statistics
import matplotlib.pyplot as plt





def dataFiltering(df:pd.DataFrame, missense_count: int)->pd.DataFrame:
    """My instinct is to have a function that generates a dataframe for me to apply PCA to, this dataframe's
    only adjustable parameter is the missense_count. Does changing the missense count change anything about the data?
    if it does, or doesn't, it shouldn't matter this function is still useful"""

    missense_counts = df[df["record_type"] == "missense"].groupby('Position').size()

    print(missense_counts.iloc[0])

    pass




reference_sequence = ""
for record in SeqIO.parse("pten-reference-protein.fasta", "fasta"):   #gives the referencesequence. 
    reference_sequence += str(record.seq)



abudanceDF = pd.read_csv("pten-variant-abundance.csv")

#step 1: we need to add a column that parses the hgvs_pro column and creates a position column out of it

abudanceDF["Position"] = abudanceDF["hgvs_pro"].str.extract(r"(\d+)")

#step 2: we want to keep only positions that have a missense count >= 5

dataFiltering(abudanceDF,5)


