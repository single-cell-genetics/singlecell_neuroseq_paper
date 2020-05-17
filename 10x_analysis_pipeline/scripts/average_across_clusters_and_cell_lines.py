import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--exprs_file')
parser.add_argument('--cell_metadata_file')
parser.add_argument('--outfile_prefix')

args = parser.parse_args()
exprs_file = args.exprs_file
metadata_file = args.cell_metadata_file
outfile_prefix = args.outfile_prefix

eDF = pd.read_csv(exprs_file, sep='\t', index_col=0)
pDF = pd.read_csv(metadata_file, sep='\t', index_col=0)

pDF = pDF.rename(columns={'res.0.4':'cluster_id', 'donor.id':'donor_id'})

assert(all([x in pDF.columns for x in ['cluster_id','donor_id']]))

df = pDF[['cluster_id','donor_id']].join(eDF.transpose(), how='inner')



cmean_df = df.drop(['donor_id'],axis=1).groupby(['cluster_id']).mean().transpose()
dmean_df = df.drop(['cluster_id'],axis=1).groupby(['donor_id']).mean().transpose()
bothmean_df = df.groupby(['cluster_id','donor_id']).mean().transpose()

cmean_df.to_csv(outfile_prefix+'.cellcluster_means.tsv', sep='\t')
dmean_df.to_csv(outfile_prefix+'.cellline_means.tsv', sep='\t')
bothmean_df.to_csv(outfile_prefix+'.cellline_cellcluster_means.tsv', sep='\t')
