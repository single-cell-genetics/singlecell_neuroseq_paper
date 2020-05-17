import os
import pandas as pd
import scanpy.api as sc
import argparse
from functools import reduce

raw_directory = '/hps/nobackup2/stegle/users/dseaton/hipsci/singlecell_neuroseq/data'
directory = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data'
sample_tracking_file = os.path.join(directory, 'metadata/sample_tracking.tsv')

parser = argparse.ArgumentParser()
parser.add_argument('--input_file',
                    help="Input file", type=str)
parser.add_argument('--dataset_outfile',
                    help="Output scanpy h5 results file", type=str)
parser.add_argument('--pca_outfile',
                    help="Output PCA file", type=str)
args = parser.parse_args()


input_file = args.input_file
dataset_outfile = args.dataset_outfile
pca_outfile = args.pca_outfile


filter_genes = True
n_top_genes = 3000
n_pcs = 50

adata = sc.read(input_file)

print('{} observations, {} genes'.format(adata.n_obs, adata.n_vars))

print('Modified recipe of Zheng 2017')


sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')


if filter_genes:
    print('Filter genes by counts')
    sc.pp.filter_genes(adata, min_cells=0.005*len(adata.obs.index))  # only consider genes expressed in more than 0.5% of cells
    print('{} observations, {} genes'.format(adata.n_obs, adata.n_vars))
    
    print('Filtering for variable genes...')
    filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
    adata = adata[:, filter_result.gene_subset]     # subset the genes
    
    print('{} observations, {} genes'.format(adata.n_obs, adata.n_vars))
    sc.pp.normalize_per_cell(adata)          # renormalize after filtering


sc.pp.scale(adata)                       # scale to unit variance and shift to zero mean


print('PCA')
sc.tl.pca(adata, n_comps=n_pcs)

pca_df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names, columns = ['PC{}'.format(x) for x in range(1,n_pcs+1)])

pca_df.to_csv(pca_outfile, sep='\t')

adata.write(dataset_outfile)
