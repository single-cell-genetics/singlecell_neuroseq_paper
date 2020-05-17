import pandas as pd
import sys

demuxlet_file = sys.argv[1]
out_file = sys.argv[2]

df = pd.read_csv(demuxlet_file, sep='\t')

df = df[['BARCODE','BEST']]

df = df[df['BEST'].apply(lambda x: x.startswith('SNG-'))]

df['donor_id'] = df['BEST'].apply(lambda x: x.replace('SNG-',''))
df['cell_id'] = df['BARCODE']

df = df[['cell_id','donor_id']]

df.to_csv(out_file, sep='\t', index=False)
