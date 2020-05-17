import scanpy.api as sc
import tempfile
import os
import pandas as pd
import argparse

''' Create 'pseudobulk' expression dataframes as flat files, one for each celltype. '''

groupby_function = 'mean'
cluster = None
dataset_name = 'pool1_9_D52'
input_file = '../../data/data_processed/{}.scanpy.w_metadata.scanpy.h5'.format(dataset_name)
outfile_prefix = input_file.replace('.h5','')
subsample = None

parser = argparse.ArgumentParser()
parser.add_argument('--input', help="Input h5 file", type=str)
parser.add_argument('--groupby', help="List of factors to group cells by", type=str,
                    nargs='+',
                    default=['donor_id','pool_id','cluster_id'])
parser.add_argument('--output_prefix', help="Output prefix", type=str)
parser.add_argument('--function', help="Function to apply to groups of cells (mean or sum)", type=str,
                    default=groupby_function)
parser.add_argument('--selected_subset', help="Selected subset to filter cells to.", type=str,
                    default=cluster)
parser.add_argument('--subset_field', help="Field on which cell filter will be applied.", type=str,
                    default='cluster_id')
args = parser.parse_args()

outfile_prefix = args.output_prefix
input_file = args.input
groupby_categories = args.groupby
groupby_function = args.function
selected_subset = args.selected_subset
subset_field = args.subset_field

if groupby_function not in ['mean', 'sum']:
    raise ValueError

print('Grouping by: {}'.format(groupby_categories))

data = sc.read(input_file, backed='r+')
print('complete dataset:')
print(data)


if selected_subset is not None:
    selected_cells = data.obs[subset_field].apply(lambda x : str(x)==selected_subset)
    data = data[selected_cells, :]


# write subsampled dataset to a temporary file, and load
temporary_h5_file = tempfile.NamedTemporaryFile(suffix='.h5', dir=os.path.dirname(input_file))
data.write(temporary_h5_file.name, force_dense=False)
data = None
data = sc.read(temporary_h5_file.name)


# basic data pre-processing
min_cells = int(data.n_obs*0.01)
sc.pp.filter_genes(data, min_cells=min_cells)


if groupby_function == 'mean':
    sc.pp.normalize_per_cell(data)
    sc.pp.log1p(data)

print(data)

comb_count = 0


outfile_template = outfile_prefix + '.{subset_field}.{selected_subset}.groupedby.{groupedby}.{groupby_function}.tsv'
outfile = outfile_template.format(subset_field=subset_field,
                                  selected_subset=selected_subset,
                                  groupedby='-'.join(groupby_categories),
                                  groupby_function=groupby_function)


cell_metadata_df = data.obs[groupby_categories].reset_index().sort_values(by=groupby_categories)
# create multiindex manually instead of using set_index, so that we can accomodate groupby 1 category
cell_metadata_df.index = pd.MultiIndex.from_arrays([cell_metadata_df[x].values for x in groupby_categories], names=groupby_categories)

combinations = list(cell_metadata_df[~cell_metadata_df.index.duplicated()].index)
print('Number of combinations = {}'.format(len(combinations)))

#make a dataframe
selected_genes = list(data.var.index)
groupby_ds_index = groupby_categories + ['n_cells'] + selected_genes
list_of_groupby_dss = []

for combination in combinations:
    if comb_count%10 == 0:
        print(comb_count)
    selected_cells = cell_metadata_df.loc[combination, 'index']
    n_cells = len(selected_cells)
    if n_cells<10:
        continue
    sub_data = data[selected_cells, :]
    sub_data = sub_data[:, selected_genes]
    if groupby_function=='mean':
        values = sub_data.X.mean(axis=0).tolist()[0]
    elif groupby_function=='sum':
        values = sub_data.X.sum(axis=0).tolist()[0]

    groupby_ds = pd.Series(list(combination) + [n_cells] + values, index=groupby_ds_index)
    list_of_groupby_dss.append(groupby_ds)
    comb_count += 1
    
groupby_df = pd.DataFrame(list_of_groupby_dss)

groupby_df.to_csv(outfile, sep='\t', index=False)
