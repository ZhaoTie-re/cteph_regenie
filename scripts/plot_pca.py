# usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import argparse

parser = argparse.ArgumentParser(description='PCA plot')
parser.add_argument('--pca_path', type=str, help='Path to the PCA projection file')
parser.add_argument('--eigenval_path', type=str, help='Path to the eigenval file')
args = parser.parse_args()

# %%
pca_path = args.pca_path
eigenval_path = args.eigenval_path

# %%
import pandas as pd

pca = pd.read_csv(pca_path, sep="\t", header=0)

# %%
with open(eigenval_path, 'r') as file:
    eigenvals = [float(line.strip()) for line in file]

total_variance = sum(eigenvals)

explained_variance_ratio = [val / total_variance for val in eigenvals]
first3_explained_variance = sum(explained_variance_ratio[:3])


# %%
import plotly.express as px
import plotly.io as pio

pca['group'] = pca['#FID'].apply(lambda x: 'NAGAHAMA' if x.startswith('NAG') else ('CTEPH' if x.startswith('PHOM') else 'Other'))

fig = px.scatter_3d(pca, x='PC1_AVG', y='PC2_AVG', z='PC3_AVG', color='group',
                    labels={'PC1_AVG': 'PC1', 'PC2_AVG': 'PC2', 'PC3_AVG': 'PC3'},
                    title=f'PCA 3D Scatter Plot (Total Explained Variance: {first3_explained_variance:.2%})',
                    hover_name='#FID',
                    hover_data={'PC1_AVG': False, 'PC2_AVG': False, 'PC3_AVG': False, 'group': False},
                    color_discrete_map={'NAGAHAMA': 'blue', 'CTEPH': 'red', 'Other': 'gray'},
                    opacity=0.7)

fig.update_traces(marker=dict(size=4))

fig.update_layout(scene=dict(
                    xaxis_title=f'PC1 ({explained_variance_ratio[0]:.2%})',
                    yaxis_title=f'PC2 ({explained_variance_ratio[1]:.2%})',
                    zaxis_title=f'PC3 ({explained_variance_ratio[2]:.2%})'),
                  margin=dict(l=10, r=10, b=10, t=40),
                  legend=dict(title='', itemsizing='constant', x=1, y=0.5, xanchor='left', yanchor='middle'),
                  width=800,
                  height=600)

# fig.show()
pio.write_html(fig, file='pca_3d.html')


