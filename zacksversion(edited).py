import pandas as pd
from Bio import SeqIO
import re
import statistics
import matplotlib.pyplot as plt


"""This code sucks bro, Its basically a summary of me learning how to abuse pandas"""


# aminoacid_dict = { #this was for old parser i THink
#     "Ala": "A",
#     "Arg": "R",
#     "Asn": "N",
#     "Asp": "D",
#     "Cys": "C",
#     "Gln": "Q",
#     "Glu": "E",
#     "Gly": "G",
#     "His": "H",
#     "Ile": "I",
#     "Leu": "L",
#     "Lys": "K",
#     "Met": "M",
#     "Phe": "F",
#     "Pro": "P",
#     "Ser": "S",
#     "Thr": "T",
#     "Trp": "W",
#     "Tyr": "Y",
#     "Val": "V"
# }




def assumptionChecks(df: pd.DataFrame, override: bool=False) -> None:
    '''Ensures assumptions are met before allowing computation
    
    Parameters:
    df: a dataframe
    
    Returns:
    None
    '''

    if override:
        return

    assert sum(df.columns.str.contains("_count")) > 0, "File must contain counts for a greater than zero number of bins"
    assert "replicate" in df.columns, "File must contain a column designating replicate numbers for the variants"
    assert "variant_type" in df.columns, "File must contain a column designating the type of the variant for each replicate"
    assert df.groupby('replicate')['variant_type'].apply(lambda x: x.eq('wild_type').any()).all(), "File must contain a wildtype variant for each replicate"
    return

def getFrequencies(df:pd.DataFrame, override: bool=False) -> pd.DataFrame:
    '''Computes the frequencies of the variant in each bin as defined by the following formula:
    Ci / ∑(Cn)

    These variables are defined as follows:
    - Ci = the count of the variant in that bin
    - Cn = the count of the nth variant in that bin

    Parameters:
    df: a dataframe containing counts of cells in different fluorescent bins

    Returns:
    Dataframe: an updated dataframe containing the frequencies of the variant in each bin
    '''

    assumptionChecks(df, override)
    
    for i in range(sum(df.columns.str.contains("_count"))):
        df[f"bin{i+1}_freq"] = df[f"bin{i+1}_count"] / df.groupby("replicate")[f"bin{i+1}_count"].transform("sum")
    return df

def getWeightedScores(freqDF: pd.DataFrame, weights:list[int], override: bool=False) -> pd.DataFrame:
    '''Computes the weighted scores as defined by the following formula:
    ∑(Wn*Fn)/∑(Fn)

    These variables are defined as follows:
    - Fn = The frequency of the variant in the nth bin
    - Wn = The weight of the variant in the nth bin

    Parameters:
    freqDF: a dataframe already with the frequency calculations present
    weights: a list of the weights corresponding to each bin in order

    Returns:
    DataFrame: an updated dataframe containing the weighted scores
    '''

    assumptionChecks(freqDF, override)
    assert sum(freqDF.columns.str.contains("_freq")) == sum(freqDF.columns.str.contains("_count")), "File must have the same number of frequency columns as its count columns"

    freqDF["weighted_score"] = (freqDF.loc[:, freqDF.columns.str.contains("_freq")].mul(weights, axis="columns").sum(axis=1) / 
                                freqDF.loc[:, freqDF.columns.str.contains("_freq")].sum(axis=1))

    return freqDF

def getNormalizedScores(weightedDF: pd.DataFrame, override: bool=False) -> pd.DataFrame:
    '''Computes the normalized weighted scores as defined by the following formula:
    (Wv - Mns) / (Wwt - Mns)
    
    These variables are defined as follows:
    - Wv  = the weighted score of the variant
    - Mns = the median of the weighted scores of the nonsense variants
    - Wwt = the weighted score of the wildtype

    Parameters:
    weightedDF: a dataframe already with weighted scores present

    Returns:
    DataFrame: an updated dataframe containing the normalized scores
    '''

    assumptionChecks(weightedDF, override)
    assert "weighted_score" in weightedDF.columns, "File must contain a column with weighted scores for its variants"

    weightedDF["normalized_score"] = ((weightedDF["weighted_score"] - 
                                       weightedDF.groupby("replicate")["weighted_score"].transform(lambda x: weightedDF.loc[x.index, "weighted_score"]
                                                                                                   [weightedDF.loc[x.index, "variant_type"]=='nonsense'].median())) / 
                                      (weightedDF.groupby("replicate")["weighted_score"].transform(lambda x: weightedDF.loc[x.index, "weighted_score"]
                                                                                                   [weightedDF.loc[x.index, "variant_type"]=='wild_type'].values[0]) - 
                                       weightedDF.groupby("replicate")["weighted_score"].transform(lambda x: weightedDF.loc[x.index, "weighted_score"]
                                                                                                   [weightedDF.loc[x.index, "variant_type"]=='nonsense'].median())))
    return weightedDF


def generate_coverage(ref_seq,df: pd.DataFrame)->pd.DataFrame:
    """Return a df? Yes we can"""
    #NOTE: the abundanceDF is globally updated to have a new row I believe...

    
    df["reference_position"] = df["hgvs_pro"].str.extract(r'(\d+)') #only captures position

    missense_list = [0]*403
    print(df["reference_position"].iloc[3789])

    two_count = 0
    #iter_count = 2

    for position,mut_type in zip(df["reference_position"], df["record_type"]):
        
        if(mut_type == "missense"):
            missense_list[int(position)-1] += 1
            

    

        
        # if(position == "2"):
        #     two_count += 1
        #     print(f'At index({iter_count} there is a 2)')

        # iter_count += 1

    missenseDF = pd.DataFrame(missense_list, columns= ["missense_count_at_index"])
    
    #print(f'This is two_count: {two_count}')

    return missenseDF #has two columns: one with a position index
    

#idk what it should output yet
def get_all_statistics(df_abundance:pd.DataFrame, list_positions:list)->dict:
    """The input dataframe should be the scorestats dataframe. Specifically the one which has missensecount >= 5. 
    So we only want rows that meet that condition. Lets try to update df_statistics"""

    
    statsdict = dict(zip((list_positions),([[]*len(list_positions)])))
    #i apologize in advance, there is probably a better way of doing this. I am just ass at oneliners :/ 
    #if a row has a missense value and is present in the list
    

    for position in list_positions:

        tempDF = df_abundance[
        (df_abundance["record_type"] == "missense") &
        (df_abundance["reference_position"] == str(position))
        ]
        #print(tempDF)
        statsdict[int(position)] = tempDF["score"].to_list()
        
    

    
        
    

        


    




    return statsdict



def getVariantCountPerRegion(region_dict:dict)->list:

    
    list_to_return = []
    string_to_print=""
    for key in region_dict:
        item_count = 0
        for item in region_dict[key]:
            list_to_return.append(item)
            item_count += 1
        string_to_print = string_to_print + str(item_count) + "  "
    
    #print(string_to_print) #for debug
    return list_to_return


#---------------------------------------------------------------------------------------------------------Main---------------------------------------------------------------------------------------------------#



df2 = pd.read_csv("toy_vampseq_bin_counts.csv")

freqs = getFrequencies(df2, override=True)
wts = getWeightedScores(freqs, [0.25, 0.5, 0.75, 1], override=True)
normscores = getNormalizedScores(wts, override=True)
normscores.to_csv("normalized_scores.csv", index=False)


#Just getting referenceseq from the fasta file. 
reference_sequence = ""
for record in SeqIO.parse("pten-reference-protein.fasta", "fasta"):
    reference_sequence += str(record.seq)

abundanceDF = pd.read_csv("pten-variant-abundance.csv")

labelDF = abundanceDF["hgvs_pro"]

print(ord("B"))
print(ord("9")) 


"""This is was my parser lmao"""
# itercount = 1
# for label in labelDF:
#     if label[0] == "_":
#         print(f"WT label at index {itercount} of the csv")
#         itercount += 1
#         continue

#     #rint(itercount)
#     aminoacidWildtype = label[2:5]
#     substring = label[5:-1]

#     count = 0
#     for chr in substring:
#         if(ord(chr) > 57):
#             break
#         count += 1
#     position = (label[5: 5+count])
#     aminoacidMutant = label[5+count:]

#     #print(f"{aminoacidWildtype}     {position}     {aminoacidMutant}")
    
#     if(aminoacid_dict[aminoacidWildtype] != reference_sequence[int(position)-1]):
#         print(f"{aminoacidWildtype}     {position}     {aminoacidMutant}")
#         print("Doesn't match somehow?")


#     itercount += 1

print("goodbyeworld")


position_list = [x+1 for x in range(len(reference_sequence))] 


df3 = generate_coverage(reference_sequence,abundanceDF)
df3["position"] = position_list




scorestatisticsDF = df3[df3["missense_count_at_index"] >= 5]

positions_tocompute_list = scorestatisticsDF["position"].to_list()


#print(scorestatisticsDF)
dict_to_plot = get_all_statistics(abundanceDF,positions_tocompute_list)

mean_list = [statistics.mean(dict_to_plot[x]) for x in positions_tocompute_list]
#print(mean_list)


#separating it into regions

unannotated_region = {} #1-20
phosphatase_domain = {}#21-184
c2_domain = {}#185-350
cterminal_tail = {}#351-403

for key in dict_to_plot:
    if (0 < key <= 14):
        unannotated_region[key] = dict_to_plot[key]
    if (14 < key <= 184):
        phosphatase_domain[key] = dict_to_plot[key]
    if(184 < key <= 350):
        c2_domain[key] = dict_to_plot[key]
    if(350 < key <= 403):
        cterminal_tail[key] = dict_to_plot[key]



#1. summarize the number of included variants
#my educated guess is that this is purely the number of variants? We can get it by lening all the lists?

print(unannotated_region.keys())
print(len(getVariantCountPerRegion(unannotated_region)))
print("\n")


print("SUMMARY OF INCLUDED VARIANTS")
print(f'number of included variants in the unannotated region = {len(getVariantCountPerRegion(unannotated_region))}')
print(f'number of included variants in the phosphatase_domain = {len(getVariantCountPerRegion(phosphatase_domain))}')
print(f'number of included variants in the c2_domain = {len(getVariantCountPerRegion(c2_domain))}')
print(f'number of included variants in the cterminal_tail = {len(getVariantCountPerRegion(cterminal_tail))}')
print("\n")
#2. total amount of sites per region
#which positions are in the 
print("TOTAL AMOUNT OF SITES PER REGION")
print(f'number of  sites in unannotated region = {len(unannotated_region.keys())}')
print(f'number of  sites in phosphatase_domain = {len(phosphatase_domain.keys())}')
print(f'number of  sites in c2domain = {len(c2_domain.keys())}')
print(f'number of  sites in cterminal_tail = {len(cterminal_tail.keys())}')
print("\n")





#median of all variants per region:
#maybe I convert it back to dataframes?

position_list2 = []
score_list = []
region_list = []
for key in dict_to_plot:
    position_list2 = position_list2 + [key] * len(dict_to_plot[key])
    score_list = score_list + dict_to_plot[key]
    if key in unannotated_region:
        region_list += ["unannotatedRegion"] * len(dict_to_plot[key])
    if key in phosphatase_domain:
        region_list += ["phosphataseRegion"] * len(dict_to_plot[key])
    if key in c2_domain:
        region_list += ["c2Domain"] * len(dict_to_plot[key])
    if key in cterminal_tail:
        region_list += ["cterminalTail"] * len(dict_to_plot[key])



print(len(position_list2))
print(len(score_list))


annotatedDF = pd.DataFrame()
annotatedDF["site_position"] = position_list2
annotatedDF["score"] = score_list
annotatedDF["Region"] = region_list

#print(annotatedDF)

#ok now I want a df only from the phosphatase region

phosphataseDF = annotatedDF.loc[annotatedDF["Region"] == "phosphataseRegion"]

print(phosphataseDF)

phosphataseDF.to_csv("phosphatase_data.csv")


'''This is the scatterplot figure stuff'''
# fig = plt.figure(figsize=(10, 5))
# fig.canvas.manager.set_window_title("Variant Mean Scatter Plot")
# plt.scatter(positions_tocompute_list, mean_list, alpha=0.7)

# plt.xlabel("Position")
# plt.ylabel("Mean score")
# plt.title("Mean Score by Position")

# plt.ylim(0, 1.4)
# plt.grid(alpha=0.2)

# plt.show()










    
    




