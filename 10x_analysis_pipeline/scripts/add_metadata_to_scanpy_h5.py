import argparse
import os
import scanpy.api as sc
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--h5_input',
                    help="Input h5 file", type=str)
parser.add_argument('--h5_output',
                    help="Output results file", type=str)
parser.add_argument('--cluster_h5_input',
                    help="Input cluster h5 file", type=str)
parser.add_argument('--cluster_tsv_input',
                    help="Input cluster tsv file", type=str)
args = parser.parse_args()

data_file = args.h5_input
dataset_outfile = args.h5_output
clustered_data_file = args.cluster_h5_input
cluster_file = args.cluster_tsv_input


sample_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/metadata/sample_tracking.tsv'

if os.path.exists(clustered_data_file):
    cluster_df = pd.read_csv(cluster_file, sep='\t', index_col=0)
else:
    cluster_df = None

sample_df = pd.read_csv(sample_file, sep='\t', index_col=0)
sample_df = sample_df.set_index('sanger_sample_id')


# load cluster information, umap
clustered_data = sc.read(clustered_data_file, backed='r')

umap_mat = clustered_data.obsm['X_umap']
umap_cell_index = clustered_data.obs.index
umap_df = pd.DataFrame(data=umap_mat, index=umap_cell_index)

# free up memory
clustered_data = None



# load sc data

data = sc.read(data_file)

# augment obs with sample and cluster information
if cluster_df is not None:
    data.obs['cluster_id'] = cluster_df.loc[data.obs.index, 'louvain']

for col in ['time_point', 'pool_id', 'treatment']:
    data.obs[col] = data.obs['sample_id'].apply(lambda x: sample_df.loc[x,col])

# add umap to the file
data.obsm['X_umap'] = umap_df.loc[data.obs.index,:].values

# minimal normalisation
sc.pp.normalize_per_cell(data, key_n_counts='n_counts_all')

# write to file
data.write(dataset_outfile)
