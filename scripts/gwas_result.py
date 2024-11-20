# usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description='GWAS result analysis')
parser.add_argument('--model_name', type=str, help='Model name')
parser.add_argument('--plink_result', type=str, help='plink result file')
args = parser.parse_args()

# %%
import gwaslab as gl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
model_name = args.model_name
plink_result_path = args.plink_result

# %%
sumstats = gl.Sumstats(plink_result_path, fmt="plink2")

# %% [markdown]
# # Get the significant SNP

# %%
sig_snp = sumstats.get_lead(sig_level=5e-8, anno=True, build="38")
sig_snp.to_csv(f'{model_name}_sig_snp.csv', index=False)

# %% [markdown]
# # Manhattan and QQ plot

# %%
plt.style.use('default')
sumstats.plot_mqq(sig_level=5e-8,anno="GENENAME", build="38",mode="qqm",save=f"{model_name}_mqq.pdf",save_args={"dpi":300})

