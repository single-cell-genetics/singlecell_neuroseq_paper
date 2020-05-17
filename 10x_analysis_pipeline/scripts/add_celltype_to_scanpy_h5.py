import argparse
import os
import scanpy.api as sc
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--h5_input',
                    help="Input h5 file", type=str)
parser.add_argument('--h5_output',
                    help="Output results file", type=str)
parser.add_argument('--celltype_mapping_file',
                    help="Input mapping of clusters to celltypes tsv file", type=str)
args = parser.parse_args()

data_file = args.h5_input
dataset_outfile = args.h5_output
celltype_mapping_file = args.celltype_mapping_file

# load celltype mapping
mapping_df = pd.read_csv(celltype_mapping_file, sep='\t')
mapping_df = mapping_df.set_index('cluster_id')

# load sc data

data = sc.read(data_file)

data.obs['celltype'] = data.obs['cluster_id'].apply(lambda x: mapping_df.loc[x,'celltype'])

# write to file
data.write(dataset_outfile)
