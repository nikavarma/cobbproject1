import pandas as pd
import kagglehub
import os # Imported to combine folder paths smoothly
import numpy as np

# Download latest version
path = kagglehub.dataset_download("shuofxz/titanic-machine-learning-from-disaster")
print("Path to dataset files:", path)

# 1. Combine the download folder path with the filename 'train.csv'
train_file_path = os.path.join(path, "train.csv")

# 2. Load the training dataset into pandas
titanic = pd.read_csv(train_file_path)

# 3. Look at the first 5 rows to confirm the 'Survived' column is there
print(titanic.head())

#print(titanic.info)

ages = titanic["Age"]

# .shape is a pandas attribute. Do not use parenthesis for attributes

olderthan35 = titanic[titanic["Age"] < 35]
#.shape is similar to a len operation for a specific df

print(titanic)

# i might be able to find which columns have these NaN values

#1. How many rows and columns are there, and which columns have missing values? How many are missing in each?

print(titanic.shape) #gives rows and columns
print(titanic.isna().sum()) #literally gives u the run down its insane ansh



#2.1 overall survival rate by gender

#woman
survived_women = titanic.loc[(titanic["Sex"] == "female") & (titanic["Survived"] > 0) ]
total_women = titanic.loc[(titanic["Sex"] == "female")]
print(f"This is the survival rate of women on the titanic: {len(survived_women) / len(total_women)}")

print("\n")
#man

survived_men = titanic.loc[(titanic["Sex"] == "male") & (titanic["Survived"] > 0) ]
total_men = titanic.loc[(titanic["Sex"] == "male")]
print(f"This is the survival rate of men on the titanic: {len(survived_men) / len(total_men)}")
print("\n")

#2.2 survival rate by class and gender

#women 1st class survival rate

survivedfirstclasswomenSR = len(survived_women.loc[survived_women["Pclass"] == 1]) / len(total_women)
survivedsecondclasswomenSR = len(survived_women.loc[survived_women["Pclass"] == 2]) /len(total_women)
survivedthirdclasswomenSR = len(survived_women.loc[survived_women["Pclass"] == 3])/ len(total_women)

print(f'These are the survival for first, second, and third class women respectively: {survivedfirstclasswomenSR:.2f}, {survivedsecondclasswomenSR:.2f}, {survivedthirdclasswomenSR:.2f}')

survivedfirstclassmenSR = len(survived_men.loc[survived_men["Pclass"] == 1]) / len(total_men)
survivedsecondclassmenSR = len(survived_men.loc[survived_men["Pclass"] == 2]) /len(total_men)
survivedthirdclassmenSR = len(survived_men.loc[survived_men["Pclass"] == 3])/ len(total_men)


print(f'These are the survival for first, second, and third class men respectively: {survivedfirstclassmenSR:.2f}, {survivedsecondclassmenSR:.2f}, {survivedthirdclassmenSR:.2f}')
print("\n")

#2.3 overall survival rate if u are a person I guess:

overallsurvivalrate = len(titanic.loc[titanic["Survived"] >  0]) / len(titanic)
print(f'This is the chances of a person surviving the titanic sinking: {overallsurvivalrate}')



#3. Fill the missing Age values with the median age. Then, how many passengers over 60 survived?



#median age

medianAge = titanic["Age"].aggregate(["median"])
medianAge = medianAge.iloc[0] #median age is an it now

titanic["Age"] = titanic["Age"].fillna(medianAge)


#4. Create an AgeGroup column (Child, Teen, Adult, Senior) and find the survival rate for each group.

conditions = [
    titanic["Age"] < 12, #u are child
    titanic["Age"] < 21, #u are teen
    titanic["Age"] < 65,   #u are adult
    titanic["Age"] >= 65    #u are senior

]

choices = ["Child", "Teen", "Adult", "Senior"]


titanic["AgeGroup"] = np.select(conditions, choices)

titanic.to_csv("titanic_withagegroup.csv")


