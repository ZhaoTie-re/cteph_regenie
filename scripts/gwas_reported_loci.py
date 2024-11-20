# usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

parser = argparse.ArgumentParser(description='GWAS reprted loci analysis')
parser.add_argument('--model_name', type=str, help='Model name')
parser.add_argument('--plink_result', type=str, help='plink result file')
parser.add_argument('--report_loci', type=str, help='Reported loci file')
args = parser.parse_args()

# %%
import gwaslab as gl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
plink_file_path = args.plink_result
report_loci_path = args.report_loci
model_name = args.model_name

# %%
report_loci = pd.read_excel(report_loci_path)
report_snp = report_loci['ID'].values

sumstats = gl.Sumstats(plink_file_path,fmt="plink2")
# %% [markdown]
# # Reported SNPs present in the results

# %%
data = sumstats.data
filtered_data = data[data['SNPID'].isin(report_snp)]
filtered_data = filtered_data.merge(report_loci[['ID', 'GENE']], left_on='SNPID', right_on='ID', how='left')
filtered_data.drop(columns=['ID'], inplace=True)
cols = list(filtered_data.columns)
cols.insert(1, cols.pop(cols.index('GENE')))
filtered_data = filtered_data[cols]
filtered_data.to_csv(f'{model_name}_reported_snp.csv', index=False)

# %%
for index, row in filtered_data.iterrows():

    snpid = row['SNPID']
    chrom = row['CHR']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    
    region_ref = [snpid]
    region = (chrom, pos - 100000, pos + 100000)
    
    sumstats.plot_mqq(
        mode="r",
        anno=True,
        region_ref=region_ref,
        region=region,
        region_grid=True,
        build="38",
        vcf_path=gl.get_path("1kg_eas_hg38"),
        save=f"{model_name}_region_{snpid}.pdf",
        save_args={"dpi":300}
    )


