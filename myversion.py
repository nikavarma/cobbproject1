import pandas as pd
from Bio import SeqIO
import re
import statistics
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



#TODO make all three into one function later tonight
def dataFiltering(df:pd.DataFrame, missense_count: int)->pd.DataFrame:
    """My instinct is to have a function that generates a dataframe for me to apply PCA to, this dataframe's
    only adjustable parameter is the missense_count. Does changing the missense count change anything about the data?
    if it does, or doesn't, it shouldn't matter this function is still useful"""

    missense_counts = df[df["record_type"] == "missense"].groupby('Position').size()

    valid_positions = missense_counts[missense_counts >= missense_count].index
    
    filtered_df = df[df["Position"].isin(valid_positions)].copy()
    

    return filtered_df

def dataFiltering2(df: pd.DataFrame, SynCount: int, NonCount: int) -> pd.DataFrame:

    non_counts = df[df["record_type"] == "nonsense"].groupby("Position").size()

    valid_positions = non_counts[non_counts >= NonCount].index
    
    filtered_df = df[df["Position"].isin(valid_positions)].copy()
    
    syn_counts = filtered_df[filtered_df["record_type"] == "synonymous"].groupby("Position").size()
    
    valid_positions = syn_counts[syn_counts >= SynCount].index

    filtered_df2 = filtered_df[filtered_df["Position"].isin(valid_positions)].copy()
    
    return filtered_df2


def getRandomMissenseNumber(df: pd.DataFrame)->int:
    """Iterate through a column and find out which position has the minimum number of missense mutations
    Input: DF should have a position and record_typecolumn, and it should be sorted by position(ascending=true)"""
    lowest = 0
    last_position = 0
    count = 0
    for row in df.itertuples():
        
        
        if(last_position == 0):
            last_position = row.Position
        
        if((row.Position == last_position) & (row.record_type == "missense")):
            count += 1
        
        if(row.Position != last_position):
            
            if(lowest == 0):
                lowest = count
            
            if(lowest > count):
                lowest = count
                        
            count = 1 #reset count
            last_position = row.Position
            
    print(f"This is the last position for debug purpose {last_position}")
    return lowest
        


def generateRandomMissenseDF(df: pd.DataFrame)->int:
    lowest_missense_count = getRandomMissenseNumber(df)
    unique_positions = list(set(df["Position"].to_list())).sort()
    newDF = pd.DataFrame()
    for pos in unique_positions:
        tempDF = df[(df["record_type"] == "missense") & (df["Position"] == pos)] 
        randomsliceDF = tempDF.sample(n=lowest_missense_count,replace=True)
        synrow = df[(df["record_type"] == "synonymous") & (df["Position"] == pos)]
        nonrow = df[(df["record_type"] == "nonsense") & (df["Position"] == pos)]
        
        
        randomsliceDF = pd.concat([randomsliceDF,synrow], ignore_index=True)
        randomsliceDF = pd.concat([randomsliceDF,nonrow], ignore_index=True)
        
        newDF = pd.concat([newDF,randomsliceDF], ignore_index=True)
        
    
    return newDF

reference_sequence = ""
for record in SeqIO.parse("pten-reference-protein.fasta", "fasta"):   #gives the referencesequence. 
    reference_sequence += str(record.seq)



abudanceDF = pd.read_csv("pten-variant-abundance.csv")

#step 1: we need to add a column that parses the hgvs_pro column and creates a position column out of it

abudanceDF["Position"] = abudanceDF["hgvs_pro"].str.extract(r"(\d+)").astype("Int64")

#step 2: we want to keep only positions that have a missense count >= 5

#the number at the end represents how many missense we are counting 

filteredDF5 = dataFiltering(abudanceDF,5)
print(filteredDF5)


#step 3: lets isolate the phosphatase domain


phosphataseDF5 = filteredDF5[(13 < filteredDF5["Position"]) & (filteredDF5["Position"] < 186)]






phosphataseDF5 = phosphataseDF5.sort_values(
    by=["Position", "score"]
)



phosphataseDF5.to_csv("phosphatase5.csv")

print(f'This is the min missense score: {getRandomMissenseNumber(phosphataseDF5)}')



forPCADF = phosphataseDF5.loc[:,["score", "Position"]]


excludeSynandMissDF = dataFiltering2(phosphataseDF5,1,1)



usedpositions_list = set(excludeSynandMissDF["Position"].to_list())
print(f'These are the used positions: {usedpositions_list}\n Number of unique positions: {len(usedpositions_list)}')

excludeSynandMissDF.to_csv("phosphatase5v2.csv")



randomlySelectedDF = generateRandomMissenseDF(excludeSynandMissDF)


    
randomlySelectedDF.sort_values(
    by=["score","Position"],ascending=True)

randomlySelectedDF.to_csv("Randomlyselected.csv")



"""Then, randomly select a number of missense mutations to be kept from each position.
  This number of randomly selected mutations is equal to the number of missense mutations that the position with the lowest number of missense mutations has.
"""




# print(forPCADF)


# scaler = StandardScaler()
# scaled_data = scaler.fit_transform(forPCADF)

# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(scaled_data)

# pcaDF = pd.DataFrame(
#     data=principal_components, columns= ["PC1","PC2"]
# )





# phosphataseDF5.plot.scatter(
#     x="Position",
#     y="score",
#     title="PhosphataseScorePlot"
# )


pcaDF.plot.scatter(
    x="PC1",
    y="PC2",
    title="PCA for Phosphatase Region of PTEN"
    
    
)

# plt.show()