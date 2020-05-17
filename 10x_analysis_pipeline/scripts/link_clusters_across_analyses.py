import pandas as pd
import os
import scanpy as sc

def reindex_obs_df(df):
    df['index'] = df.index
    df['cell_barcode'] = df['index'].apply(lambda x: x.split('-')[0])
    df['new_index'] = df['cell_barcode'] + '-' + df['sample_id'].astype('str')
    df = df.set_index('new_index')
    return df

ref_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/pool1_13_noddd_subsample/pool1_13_noddd_subsample.scanpy.w_metadata.scanpy.obs_df.tsv'
file_template = '../../data/data_processed/{dataset}/{dataset}.scanpy.w_metadata.scanpy.obs_df.tsv'
ref_h5_file = ref_file.replace('.obs_df.tsv','.h5')
out_file = ref_file.replace('.obs_df.tsv', '.mapped_to_new_clusters.obs_df.tsv')
out_h5_file = ref_h5_file.replace('.h5','.mapped_to_new_clusters.h5')

ref_df = pd.read_csv(ref_file, sep='\t', index_col=0)
ref_df = reindex_obs_df(ref_df)


datasets = ['pool1_13_noddd_'+x for x in ['D11','D30','D52']][:]


df_list = []
for dataset in datasets:
    df = pd.read_csv(file_template.format(dataset=dataset), sep='\t', index_col=0)
    df = reindex_obs_df(df)
    df['new_cluster_id'] = dataset + '_' + df['cluster_id'].astype('str')
    df_list.append(df)
#    df = df[['new_cluster_id']]
combined_df = pd.concat(df_list, axis=0)

obs_df = ref_df.join(combined_df[['new_cluster_id']],how='inner')

obs_df = obs_df.set_index('index')

obs_df.to_csv(out_file, sep='\t')


adata = sc.read(ref_h5_file)

assert(all(adata.obs.index==obs_df.index))
adata.obs = obs_df

adata.write(out_h5_file)
