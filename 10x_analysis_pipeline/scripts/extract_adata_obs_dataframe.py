import scanpy.api as sc
import pandas as pd
import argparse

'''Extract the obs dataframe (with cell metadata) from a scanpy h5 file.'''

parser = argparse.ArgumentParser()
parser.add_argument('--input',
                    help="Input scanpy h5 file with expression data", type=str)
parser.add_argument('--output', help="Output file", type=str)
args = parser.parse_args()

input_file = args.input
outfile = args.output

# load the dataset, extract the obs dataframe, and then free the memory
data = sc.read(input_file, backed='r+')
print('complete dataset:')
print(data)

obs_df = data.obs.copy()
data = None

obs_df.to_csv(outfile, sep='\t')
