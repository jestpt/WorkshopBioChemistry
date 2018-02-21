import pandas as pd 

# So students don't have to deal with different data types

dataset = pd.read_csv('kag_risk_factors_cervical_cancer.csv')
for column in dataset:
	if dataset[column].dtypes == object:
		dataset[column] = dataset[column].loc[dataset[column] != '?']

dataset = dataset.astype(float)
dataset.to_csv('cervical_cancer.csv')