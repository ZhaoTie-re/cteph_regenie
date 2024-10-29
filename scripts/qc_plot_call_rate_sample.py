# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse
import os

parser = argparse.ArgumentParser(description='Plot sample missing rate distribution')
parser.add_argument('--chr', type=str, help='chromosome number')
parser.add_argument('--sample_miss_path', type=str, help='sample missing rate file path')
args = parser.parse_args()

# %%
chr = args.chr
sample_miss_path = args.sample_miss_path

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sample_miss = pd.read_csv(sample_miss_path, sep="\s+")

# %%
plt.style.use('default')
plt.figure(figsize=(8, 5))
plt.hist(sample_miss['F_MISS'].dropna(), bins=30, edgecolor='k', alpha=0.7)
plt.title('Sample missing rate distribution of Chr' + chr)
plt.xlabel('Sample missing rate')
plt.ylabel('Count') 
plt.grid(True)
# plt.show()
plt.savefig("sample_missing_rate_distribution_chr" + chr + ".pdf") 


