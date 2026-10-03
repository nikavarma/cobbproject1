import pandas as pd

"""
This file  "should" contain all the functions and code used for part1 of the project 
that will be reviewed in class on 9/22. It computes the weighted scores for each variant
and computes normalized scores for said variants. As for the critical thinking, um there is none
in this file. Enjoy


"""


def getFrequencyList(variantIndex, dataframe): 
    """This function should return a list of four frequencies, FOR A SPECIFIC VARIANT. A frequency of a variant in a specific bin, 
    is the count of that variant in the specific  bin divided by the total count of all variants in that bin. For example, the WT variant has 
    a bin1_count of 9 and we would divide that by the total count of all variants in bin1 which is (9 + 12 + 221 + 217...+0+4 = ~1200 )."""
    

    
    #i think this is the best way of getting the freq
    binList = ["bin1_count", "bin2_count", "bin3_count", "bin4_count"]
    freqList = []
    for somebin in binList:
        specificBin = dataframe[somebin]
        variantCount = specificBin[variantIndex]
        count = variantCount #we are going to divide by the count so we want to include the variant
        
        for i in range(len(dataframe)): #len of dataframe is always 12
            
            if i == variantIndex:
                continue
            else:
                count += specificBin[i]
        #print(count)
        freq = variantCount / count
        freqList.append(freq)
    
    
    return freqList



def calculateVariantWeightedScore(freqList):
    """This function implements the weighted score formula from the slides. W = f1*0.25 + f2*0.50 + f3*0.75 + f4*1.0 / (f1+f2+f3+f4)"""

    #list must be of length 4
    f1,f2,f3,f4 = freqList #idk if i can do this
    w1,w2,w3,w4 = [0.25,0.50,0.75,1.0]
    denom = f1+f2+f3+f4
    weightedScore = (w1*f1) + (w2*f2) + (w3*f3) + (w4*f4)
    weightedScore = weightedScore / denom
    return weightedScore


#idk if this actually normalizes, this is my interpretation of the formula on the last slide
def normalize(weightedScore, medianNonsense, wildType):
    """This function implements the normalization formula from the slides. N = (W - Wmedian) / (Wwt - Wmedian)"""

    normalizedScore = weightedScore - medianNonsense
    normalizedScore = normalizedScore / (wildType - medianNonsense)
    return normalizedScore
    

#ok i think if i just switch to interpreter it should work?



#-------------------------------------main------------------------------------------------------#


df = pd.read_csv("toy_vampseq_bin_counts.csv")
print("Hello World")





