# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse
import os

parser = argparse.ArgumentParser(description='Plot F coefficient distribution') 
parser.add_argument('--chr', type=str, help='chromosome number')
parser.add_argument('--het_path', type=str, help='het file path')
args = parser.parse_args()

# %%
chr = args.chr
het_path = args.het_path

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

het_pd = pd.read_csv(het_path, sep="\s+")

# %%
mean_F = het_pd['F'].mean()
std_F = het_pd['F'].std()

plt.style.use('default')
plt.figure(figsize=(8, 5))
plt.hist(het_pd['F'], bins=30, edgecolor='k', alpha=0.7)

plt.axvline(mean_F, color='red', linestyle='dashed', linewidth=1)
plt.axvline(mean_F + 3*std_F, color='green', linestyle='dashed', linewidth=1)
plt.axvline(mean_F - 3*std_F, color='green', linestyle='dashed', linewidth=1)
plt.axvline(0.1, color='orange', linestyle='dashed', linewidth=1)
plt.axvline(-0.1, color='orange', linestyle='dashed', linewidth=1)

plt.xlabel('F')
plt.ylabel('Count') 
plt.title('Distribution of F with 3 SD in Chr' + chr)
# plt.show()
plt.savefig(chr + ".het.pdf") 

# %%
outliers = het_pd[(het_pd['F'] > 0.1) | (het_pd['F'] < -0.1)]
outliers_to_remove = outliers[['FID', 'IID']]
outliers_to_remove.to_csv(f'{chr}.high_het.sample', sep='\t', index=False, header=False)


