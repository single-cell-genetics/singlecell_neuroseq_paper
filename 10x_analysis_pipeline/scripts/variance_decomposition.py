from __future__ import print_function
import limix
import numpy as np
import pandas as pd
import re
import limix_tools
import os.path
import random
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--outfile')
parser.add_argument('--metadatafile')
parser.add_argument('--effects',nargs='+',default=None)
parser.add_argument('--transform')
parser.add_argument('--ncells',type=int)
parser.add_argument('--seed',type=float)
args = parser.parse_args()

input_file = args.infile
output_file = args.outfile
metadata_file = args.metadatafile
effects = args.effects
variance_stabilizing_transform = args.transform
ncells = args.ncells
seed = args.seed

random.seed(seed)

print('Input: {}'.format(input_file))


df = pd.read_csv(input_file,sep='\t',index_col=0)

metadata_df = pd.read_csv(metadata_file,sep='\t',index_col=0)

#limit to effects specified by the user
if effects is not None:
    metadata_df = metadata_df.loc[:, effects]


samples = list(set(df.columns)&set(metadata_df.index))
samples = random.sample(samples, ncells)
df = df.loc[:, samples]

n_genes = df.shape[0]
n_cells = df.shape[1]

print('Number of samples: {n_samples}'.format(n_samples=len(samples)))
print('Number of features: {}'.format(df.shape[0]))

if variance_stabilizing_transform=='log2':
    transform_fcn = np.log2
elif variance_stabilizing_transform=='log2_offset':
    transform_fcn = lambda x: np.log2(x+1)
elif variance_stabilizing_transform=='none':
    transform_fcn = lambda x: x
else:
    print('Specified transform function should be log2, log2_offset, or none')
    raise(ValueError)

var_df = limix_tools.run_variance_analysis(df,metadata_df,transform_fcn=transform_fcn)

var_df.to_csv(output_file,sep='\t')
