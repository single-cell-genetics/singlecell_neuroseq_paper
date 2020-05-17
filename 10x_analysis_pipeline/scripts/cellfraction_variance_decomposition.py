from __future__ import print_function
import limix
import numpy as np
import pandas as pd
import re
import limix_tools
import os.path
import random
import argparse


input_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point.celltype_fractions_pivoted.tsv'.format(dataset='pool1_13_noddd_D11')
output_file = input_file.replace('.tsv','.var_decomp.tsv')
#metadata_file = ''
donor_filter = True
effects = None
variance_stabilizing_transform = 'none'

# parser = argparse.ArgumentParser()
# parser.add_argument('--infile')
# parser.add_argument('--outfile')
# parser.add_argument('--metadatafile')
# parser.add_argument('--effects',nargs='+',default=None)
# parser.add_argument('--transform')
# parser.add_argument('--no_donor_filter',action='store_true')
# args = parser.parse_args()

# input_file = args.infile
# output_file = args.outfile
# metadata_file = args.metadatafile
# effects = args.effects
# variance_stabilizing_transform = args.transform
# donor_filter = not args.no_donor_filter


data_df = pd.read_csv(input_file,sep='\t')

data_df['index'] = data_df['donor_id'] + ':' + data_df['pool_id']
data_df = data_df.set_index('index')
metadata_df = data_df[['donor_id','pool_id']]
data_df = data_df.drop(['donor_id','pool_id','time_point'], axis=1)
data_df = data_df.transpose()

#metadata_df = pd.read_csv(metadata_file,sep='\t',index_col=0)

if donor_filter:
    #select only lines from donors with at least 2 lines
    donor_list = metadata_df['donor_id'].unique()
    donor_counts = metadata_df['donor_id'].value_counts()
    donors_w_multiple_lines = [x for x in donor_list if donor_counts[x]>1]
    metadata_df = metadata_df[metadata_df['donor_id'].isin(donors_w_multiple_lines)]

#limit to effects specified by the user
if effects is not None:
    metadata_df = metadata_df.loc[:, effects]

samples = list(set(data_df.columns)&set(metadata_df.index))
print('Number of samples: {n_samples}'.format(n_samples=len(samples)))

data_df = data_df.loc[:, samples]
#data_df = data_df.iloc[:nProteins, :]

if variance_stabilizing_transform=='log2':
    transform_fcn = np.log2
elif variance_stabilizing_transform=='log2_offset':
    transform_fcn = lambda x: np.log2(x+1)
elif variance_stabilizing_transform=='none':
    transform_fcn = lambda x: x
else:
    print('Specified transform function should be log2, log2_offset, or none')
    raise(ValueError)

var_df = limix_tools.run_variance_analysis(data_df,metadata_df,transform_fcn=transform_fcn)

var_df.to_csv(output_file,sep='\t')
