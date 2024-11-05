# usr/bin/env python3
# -- coding: utf-8 --
# %%
import argparse

parser = argparse.ArgumentParser(description='Phenotype file preparation')
parser.add_argument('--covariate_file', type=str, help='covariate file')

args = parser.parse_args()

# %%
covariate_file = args.covariate_file

import pandas as pd
import numpy as np

covariate_df = pd.read_csv(covariate_file, sep='\t')
pheno_df = covariate_df[['#FID', 'IID', 'PHENO1']]
pheno_df.columns = pheno_df.columns.str.replace('#', '')

pheno_df.to_csv('pheno_cteph.txt', sep=' ', index=False)


