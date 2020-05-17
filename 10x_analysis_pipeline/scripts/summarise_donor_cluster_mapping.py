from __future__ import print_function
import pandas as pd
import numpy as np
import os
import seaborn as sns
from matplotlib.pyplot import *
switch_backend('agg')
import glob
import argparse

sns.set_style("ticks")

parser = argparse.ArgumentParser()
parser.add_argument('exprs_file',action='store',type=str)
parser.add_argument('metadata_file',action='store',type=str)
parser.add_argument('outdir',action='store',type=str)

args = parser.parse_args()

exprs_file = args.exprs_file
outdir = args.outdir
metadata_file = args.metadata_file

figure_dir = outdir




pDF = pd.read_csv(metadata_file, sep='\t', index_col=0)
pDF['cluster_id'] = pDF['res.0.4']




donors = sorted(pDF['donor_id'].drop_duplicates().tolist())
clusters = sorted(pDF['cluster_id'].drop_duplicates().tolist())

df = pDF.groupby('donor_id').apply(lambda x: x['cluster_id'].value_counts())

donor_cluster_df = pd.DataFrame(index=donors, columns=clusters)

for donor in donors:
    cluster_counts_ds = df.loc[donor]
    clusters_with_cells = cluster_counts_ds.index
    donor_cluster_df.loc[donor, clusters_with_cells] = cluster_counts_ds

donor_cluster_df = donor_cluster_df.loc[sorted(donors), sorted(clusters)]
donor_cluster_df = donor_cluster_df.fillna(0)

print(donor_cluster_df.head(10))



# plot cell counts

figure_title = 'summary_of_cell_counts_across_cell_lines_and_clusters.png'
figure_title = os.path.join(figure_dir,figure_title)

FS = 16

fig = figure(figsize=(5,(len(donors) + len(clusters))*0.4))


positions = list(range(len(donors))) + list(range(len(donors) + 1, len(donors) + len(clusters) + 1))
positions = [len(positions)-x for x in positions]
data = donor_cluster_df.sum(axis=1).tolist() + donor_cluster_df.sum(axis=0).tolist()

barh(positions, data)
plot([0, max(data)], [len(clusters), len(clusters)], 'k--')

yticks(positions, donors + ['cluster {}'.format(x) for x in clusters])
xlabel('N cells', fontsize=FS)

fig.savefig(figure_title, bbox_inches='tight', dpi=100)





figure_title = 'heatmaps_of_cell_counts_across_cell_lines_and_clusters.png'
figure_title = os.path.join(figure_dir,figure_title)

fig = figure(figsize=(len(clusters)*1.4,len(donors)*2.1))

FS = 16

subplot(3,1,1)
sns.heatmap(donor_cluster_df, cbar_kws={'label':'N cells'})

subplot(3,1,2)
sns.heatmap(donor_cluster_df.divide(donor_cluster_df.sum(axis=1), axis=0), 
            vmin=0,vmax=1.0,  cbar_kws={'label':'Proportion of donor cells\n(rows sum to 1)'})
ylabel('Cell line ID', fontsize=FS)

subplot(3,1,3)
sns.heatmap(donor_cluster_df/donor_cluster_df.sum(axis=0), 
            vmin=0,vmax=1.0,  cbar_kws={'label':'Proportion of cluster cells\n(columns sum to 1)'})

xlabel('Cluster ID', fontsize=FS)
tight_layout()

fig.savefig(figure_title, bbox_inches='tight', dpi=100)
