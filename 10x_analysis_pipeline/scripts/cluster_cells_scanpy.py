import argparse
import pandas as pd
import scanpy.api as sc


parser = argparse.ArgumentParser()
parser.add_argument('--h5_input',
                    help="Input h5 file", type=str)
parser.add_argument('--pca_input',
                    help="Input dimensionality reduction file", type=str)
parser.add_argument('--h5_output',
                    help="Output results file", type=str)
parser.add_argument('--cluster_output',
                    help="Output cluster tsv file", type=str)
parser.add_argument('--resolution',
                    help="Louvain clustering resolution - higher values create smaller clusters.",
                    default=0.15, type=float)
args = parser.parse_args()

in_file = args.h5_input
pca_file = args.pca_input
dataset_outfile = args.h5_output
cluster_outfile = args.cluster_output
resolution = args.resolution

pca_df = pd.read_csv(pca_file, sep='\t', index_col=0)

adata = sc.read(in_file)

# replace PCA with the specified PCA
adata.obsm['X_pca'] = pca_df.values

print('clustering...')

sc.pp.neighbors(adata, n_neighbors=10)

sc.tl.louvain(adata, flavor='vtraag',resolution=resolution)

print('{} clusters identified'.format(adata.obs['louvain'].drop_duplicates().shape[0]))

print('computing UMAP...')
sc.tl.umap(adata, random_state=2)

adata.write(dataset_outfile)

cell_clustering_df = adata.obs[['louvain']]
cell_clustering_df.to_csv(cluster_outfile, sep='\t')
