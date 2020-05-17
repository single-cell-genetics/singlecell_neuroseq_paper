import argparse
import pandas as pd

'''Count the number of cells in each category (i.e. collection of groupby
terms) specified by the user.'''

default_groupby = ['donor_id','pool_id','time_point']

parser = argparse.ArgumentParser()
parser.add_argument('--input',
                    help="Input obs_df tsv file", type=str)
parser.add_argument('--groupby',
                    type=str,
                    nargs='+',
                    help='List of factors to group by.',
                    default=default_groupby)
args = parser.parse_args()

input_file = args.input
groupby_cols = args.groupby
pivot_col = 'celltype'

groupby_label = '-'.join(groupby_cols)
outfile_prefix = input_file.replace('.tsv','.')
outfile_prefix += 'groupedby.' + groupby_label

df = pd.read_csv(input_file, sep='\t', index_col=0)

print(df.head())

df['n_cells'] = 1.0

#count the rows i.e. the number of cells
df2 = df.groupby(groupby_cols + [pivot_col]).count()[['n_cells']]

df3 = df.groupby(groupby_cols).count()[['n_cells']]

df4 = df2.apply(lambda x: x['n_cells']/df3.loc[x.name[:len(groupby_cols)],'n_cells'], axis=1)
df4 = df4.reset_index()
df4.columns = groupby_cols + [pivot_col, 'f_cells']

df2.to_csv(outfile_prefix+'.celltype_counts.tsv', sep='\t')
df3.to_csv(outfile_prefix+'.sample_counts.tsv', sep='\t')
df4.to_csv(outfile_prefix+'.celltype_fractions.tsv', sep='\t', index=False)

df5 = df4.pivot_table(index=groupby_cols,columns=[pivot_col], fill_value=0.0)['f_cells']
df5.to_csv(outfile_prefix+'.celltype_fractions_pivoted.tsv', sep='\t')
