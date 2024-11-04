# usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os

parser = argparse.ArgumentParser(description='HWE P-value Distribution')
parser.add_argument('--hwe_path', type=str, help='hwe file path')
args = parser.parse_args()

# %%
hwe_path = args.hwe_path

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

hwe = pd.read_csv(hwe_path, sep="\s+")

# %%
threshold = 1e-6
less_equal_threshold = hwe['P'].dropna() <= threshold
greater_threshold = hwe['P'].dropna() > threshold

count_less_equal = less_equal_threshold.sum()
count_greater = greater_threshold.sum()

plt.style.use('default')
plt.figure(figsize=(8, 5))
plt.bar(['<= 1e-6', '> 1e-6'], [count_less_equal, count_greater], color=['blue', 'orange'], edgecolor='k', alpha=0.7)
plt.xlabel('HWE P-value')   
plt.ylabel('Count') 
plt.title('HWE P-value distribution')
# plt.grid(True)
# plt.show()
plt.savefig("hwe_pvalue_distribution.pdf")


