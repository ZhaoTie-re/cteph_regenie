# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse

parser = argparse.ArgumentParser(description='Prepare total covariate file (add age, sex and PCs)')
parser.add_argument('--covariate_file', type=str, help='Path to the covariate file')
parser.add_argument('--age_sex_file', type=str, help='Path to the age_sex file')

args = parser.parse_args()

# %%
covariate_path = args.covariate_file
age_sex_path = args.age_sex_file

# %%
import pandas as pd
import numpy as np

covariate = pd.read_csv(covariate_path, sep="\t")
age_sex = pd.read_csv(age_sex_path, sep=" ")

# %%
age_sex.rename(columns={'ID': 'IID'}, inplace=True)

# %%
merged_df = pd.merge(covariate, age_sex, on='IID', how='left')

# %%
cols = merged_df.columns.tolist()
pc1_index = cols.index('PC1_AVG')
cols.remove('SEX_M1_F2')
cols.remove('AGE')
new_cols = cols[:pc1_index] + ['SEX_M1_F2', 'AGE'] + cols[pc1_index:pc1_index+1] + cols[pc1_index+1:]
merged_df = merged_df[new_cols]

merged_df['AGE'] = (merged_df['AGE'] - merged_df['AGE'].mean()) / merged_df['AGE'].std()


# %%
merged_df.to_csv("covariate.txt", sep=" ", index=False)


