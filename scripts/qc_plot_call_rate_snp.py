# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse
import os

parser = argparse.ArgumentParser(description='Plot snp missing rate distribution')
parser.add_argument('--chr', type=str, help='chromosome number')
parser.add_argument('--snp_miss_path', type=str, help='snp missing rate file path')
args = parser.parse_args()

# %%
chr = args.chr
snp_miss_path = args.snp_miss_path

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

snp_miss = pd.read_csv(snp_miss_path, sep="\s+")

# %%
plt.style.use('default')
plt.figure(figsize=(8, 5))
plt.hist(snp_miss['F_MISS'].dropna(), bins=30, edgecolor='k', alpha=0.7)
plt.title('SNP missing rate distribution of Chr' + chr)
plt.xlabel('SNP missing rate')
plt.ylabel('Count') 
plt.grid(True)
# plt.show()
plt.savefig("snp_missing_rate_distribution_chr" + chr + ".pdf") 