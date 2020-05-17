import pandas as pd
import scanpy.api as sc
import glob
import os
import yaml
import argparse
from functools import reduce

''' Script to merge 10x samples together.

This incorporates donor ID mapping information for each cell (doublets are excluded).

Specification of which samples to include and exclude is done with a YAML config file.'''


raw_directory = '/hps/nobackup2/stegle/users/dseaton/hipsci/singlecell_neuroseq/data'
directory = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data'
sample_tracking_file = os.path.join(directory, 'metadata/sample_tracking.tsv')


parser = argparse.ArgumentParser()
parser.add_argument('--config_file',
                    help="Input config file", type=str)
parser.add_argument('--output_file',
                    help="Output results file", type=str)
args = parser.parse_args()

config_file = args.config_file
output_file = args.output_file

# configuration - to select which samples to include
with open(config_file, 'r') as f:
    config_dict = yaml.safe_load(f)


# memory efficiency step - exclude genes on the fly if we already know which genes to include
genelist_file = None

if genelist_file is not None:
    gene_list = sc.read(genelist_file, backed='r').var_names
else:
    gene_list = None

# loading 10x sample metadata
sample_df = pd.read_csv(sample_tracking_file, sep='\t')

# select only samples that match the conditions of the 'sample_selection'
# specified in the config file.
list_of_selections = []
for colname in config_dict['sample_selection'].keys():
    selection = sample_df[colname].isin(config_dict['sample_selection'][colname])
    list_of_selections.append(selection)
complete_selection = reduce(lambda x, y: x*y, list_of_selections)

sample_df = sample_df[complete_selection]

# point to file paths for expression counts and donor mapping/demultiplexing
sample_df['path_to_h5_file'] = sample_df['sample_id'].apply(lambda x: os.path.join(raw_directory, 'data_raw/{}/filtered_gene_bc_matrices_h5.h5'.format(x)))
sample_df['path_to_demuxlet_file'] = sample_df['sample_id'].apply(lambda x: os.path.join(directory, 'data_processed/demuxlet/{}/demuxlet_out.best'.format(x)))

# complain if some files don't exist
sample_df['file_exists'] = sample_df['path_to_h5_file'].apply(lambda x: os.path.exists(x))
if sample_df['file_exists'].sum()>0:
    print('h5 files missing for: {}'.format(sample_df.query('not file_exists')['sample_id'].tolist()))
# only keep files that exist
sample_df = sample_df.query('file_exists')

# exclude poor 10x samples
excluded_samples = ['cellranger211_count_25528_5245STDY7387188_hg19-1_2_0',
                    'cellranger211_count_28882_5245STDY7825050_hg19-1_2_0',
                    'cellranger211_count_28882_5245STDY7825051_hg19-1_2_0',
                    'cellranger211_count_29749_5245STDY7982704_hg19-1_2_0',
                    'cellranger211_count_29749_5245STDY7982705_hg19-1_2_0',
                    'cellranger211_count_29749_5245STDY7982706_hg19-1_2_0',
                    'cellranger211_count_29749_5245STDY7982707_hg19-1_2_0']                    
sample_df = sample_df.query('sample_id not in @excluded_samples')

datasets = []
count = 1
# iterate over samples, checking whether it should be included, then adding it to the dataset with donor ID info
for idx,row in sample_df.iterrows():
    h5_file = row['path_to_h5_file']
    demuxlet_file = row['path_to_demuxlet_file']
    print(count)
    if True:
        sample_id = row['code']
        print('sample: ' + sample_id)

        # read the donor mapping information
        donor_df = pd.read_csv(demuxlet_file, sep='\t', index_col=0)
        donor_df['TYPE'] = donor_df['BEST'].apply(lambda x: x.split('-')[0])
        # select only Single Cells - exclude doublets and ambiguous cells ("DBL","AMB")
        print('{} out of {} cells are singlets.'.format((donor_df['TYPE']=='SNG').sum(), donor_df.shape[0]))
        donor_df = donor_df.query('TYPE=="SNG"')
        if 'donor_filter' in config_dict.keys():
            # filter out specific donors
            pool_id = row['pool_id']
            if pool_id in config_dict['donor_filter'].keys():
                # filter selected donors from this sample
                donors_to_filter_out = config_dict['donor_filter'][pool_id]
                donor_df = donor_df[~donor_df['SNG.1ST'].isin(donors_to_filter_out)]
        selected_cells = donor_df.index

        # read the count info, and merge with donor mapping
        data = sc.read_10x_h5(h5_file, genome='hg19')
        data = data[selected_cells, :]
        data.var_names_make_unique()
        # label mitochondrial genes
        data.var['mito'] = [x.startswith('MT-') for x in data.var_names]
        # calculate qc metrics
        sc.pp.calculate_qc_metrics(data, qc_vars=['mito'], inplace=True)
        # apply subsampling to cells and counts
        if 'subsample_cells' in config_dict.keys():
            sc.pp.subsample(data, fraction=config_dict['subsample_cells'])
        if 'subsample_counts' in config_dict.keys():
            fraction = config_dict['subsample_counts']
            target_counts_per_cell = data.obs['total_counts'].apply(lambda x: int(x*fraction)).values
            sc.pp.downsample_counts(data, counts_per_cell=target_counts_per_cell)
        if 'cell_qc' in config_dict.keys():
            n_cells_start = data.n_obs
            for filter_query in config_dict['cell_qc']:
                data = data[data.obs.query(filter_query).index, :]
            print('cell QC applied: {} cells dropped'.format(n_cells_start-data.n_obs))
        # drop everything except gene ids and mito from data.var
        data.var = data.var[['gene_ids','mito']]
        if count>0:
            # drop everything else from data.var (duplicated with first dataset)
            data.var = data.var.iloc[:,:0]
        # apply gene list filtering if specified
        if gene_list is not None:
            data = data[:, gene_list]
        data.obs['sample_id'] = sample_id
        data.obs['donor_id'] = donor_df.loc[data.obs.index, 'SNG.1ST']
        print('{} obs, {} vars'.format(data.n_obs, data.n_vars))
        if data.n_obs>0:
            datasets.append(data)
            count+=1


merged_data = datasets[0].concatenate(*datasets[1:])
print('{} obs, {} vars'.format(merged_data.n_obs, merged_data.n_vars))

print('write results')
merged_data.write(output_file)
