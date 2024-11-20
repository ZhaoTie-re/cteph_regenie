# usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

parser = argparse.ArgumentParser(description='GWAS reprted gene analysis')
parser.add_argument('--model_name', type=str, help='Model name')
parser.add_argument('--plink_result', type=str, help='plink result file')
parser.add_argument('--report_gene', type=str, help='Reported gene file')
args = parser.parse_args()

# %%
reported_gene_path = args.report_gene
plink_file_path = args.plink_result
model_name = args.model_name

# %%
import gwaslab as gl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

reported_gene = pd.read_excel(reported_gene_path)
sumstats = gl.Sumstats(plink_file_path,fmt="plink2")

# %%
reported_gene['CHROM'] = reported_gene['CHROM'].apply(lambda x: int(x[3:]))

# %%
for index, row in reported_gene.iterrows():

    gene = row['GENE']
    chrom = row['CHROM']
    start = row['START']
    end = row['END']
    range = end - start
    
    region = (chrom, start-100000, end+100000)
    
    plt.style.use('default')
    sumstats.plot_mqq(
        mode="r",
        anno=True,
        region=region,
        region_grid=True,
        build="38",
        vcf_path=gl.get_path("1kg_eas_hg38"),
        save=f"{model_name}_region_{gene}.pdf",
        save_args={"dpi":300}
    )


