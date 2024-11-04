# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse
import os

parser = argparse.ArgumentParser(description='MAF Distribution')
parser.add_argument('--frq_path', type=str, help='frq file path')
args = parser.parse_args()

# %%
maf_path = args.frq_path

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

maf = pd.read_csv(maf_path, sep='\s+')

# %%

threshold = 0.01
less_equal_threshold = maf['MAF'].dropna() <= threshold
greater_threshold = maf['MAF'].dropna() > threshold

count_less_equal = less_equal_threshold.sum()
count_greater = greater_threshold.sum()

plt.style.use('default')
plt.figure(figsize=(8, 5))
plt.bar(['<= 0.01', '> 0.01'], [count_less_equal, count_greater], color=['blue', 'orange'], edgecolor='k', alpha=0.7)
plt.xlabel('MAF')
plt.ylabel('Count') 
plt.title('MAF Distribution')
# plt.grid(True)
# plt.show()
plt.savefig("maf_distribution.pdf")

