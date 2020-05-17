import pandas as pd
import os
import glob

pool_df = pd.read_csv('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/metadata/pool_tracking.tsv', sep='\t')
pool_df['donor_id'] = pool_df['cell_line_id'].apply(lambda x: x.split('_')[0])

pool_id = 'pool7'
#sample_id = 'cellranger211_count_26746_5245STDY7619068_hg19-1_2_0'
sample_id = 'cellranger211_count_26746_5245STDY7619069_hg19-1_2_0'
in_files = glob.glob('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/{}/demuxlet_line_search_results/demuxlet_out.line_subset_*.best'.format(sample_id))

#pool_id = 'pool6'
#in_files = glob.glob('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/cellranger211_count_26558_5245STDY7602377_hg19-1_2_0/demuxlet_line_search_results/demuxlet_out.line_subset_*.best')
#in_files = glob.glob('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/cellranger211_count_26558_5245STDY7602379_hg19-1_2_0/demuxlet_line_search_results/demuxlet_out.line_subset_*.best')

#pool_id = 'pool4'
#in_files = glob.glob('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/cellranger211_count_25923_5245STDY7487301_hg19-1_2_0/demuxlet_line_search_results/demuxlet_out.line_subset_*.best')


list_of_dfs = []

for filename in in_files[:]:
    df = pd.read_csv(filename, sep='\t')
    df['SNG'] = df['BEST'].apply(lambda x: x.startswith('SNG'))
    df = df.query('SNG')
    list_of_dfs.append(df)

merged_df = pd.concat(list_of_dfs)


count_df = merged_df['BEST'].value_counts().to_frame()

count_df = count_df.query('BEST>5')

count_df['cell_line_id'] = [x.split('-')[-1] for x in count_df.index]

# extract just the 4-letter code for each donor
count_df['donor_id'] = count_df['cell_line_id'].apply(lambda x: x.split('_')[0])

count_df['in_experiment'] = count_df['donor_id'].apply(lambda x: x in pool_df.query('pool_id==@pool_id')['donor_id'].tolist())


pool_df = pool_df.query('pool_id==@pool_id')
pool_df['present_in_pool'] =  pool_df['donor_id'].apply(lambda x: x in count_df['donor_id'].tolist())

print(pool_df.query('not present_in_pool'))
print(count_df.query('not in_experiment'))
