import numpy as np
import pandas as pd
import os
import shutil
import sklearn.decomposition
import qtl_config_utils
import sys

phenotype_input_file = sys.argv[1]
outdir = sys.argv[2]
if len(sys.argv)>=4:
    subset_field = sys.argv[3]
    selected_subset = sys.argv[4]
else:
    subset_field = None
    selected_subset = None


outdir = os.path.abspath(outdir)

# general set up
genotypes_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Plink/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.norm.renamed.recode.vcf.gz'
kinship_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered_Plink-f/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20170327.genotypes.norm.renamed.recode.filtered.rel'
annotation_file = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/scQTLs/annos/ensembl_gene_id_annos.tsv'
chunk_file = '/nfs/leia/research/stegle/mjbonder/ChunkFiles/Ensembl_75_Limix_Annotation_FC_Gene_step100.txt'

n_pcs = 15

base_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/'

# eQTL discovery parameters
#number of top expressed genes to test
n_genes = 50000


sample_mapping_file = os.path.join(outdir, 'samplemapping.tsv')
phenotype_file = os.path.join(outdir, 'phenotypes.tsv')
covariates_file = os.path.join(outdir, 'covariates.tsv')
noise_matrix_file = os.path.join(outdir, 'noise_matrix.tsv')
config_file = os.path.join(outdir, 'qtl_config.yaml')

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

#load map HGNC to ensembl
mapping_df = pd.read_csv('/nfs/leia/research/stegle/dseaton/genomes/hg19/annotation/geneid_mappings/hgnc_symbol2ensembl_gene_id.txt', sep='\t')

# process phenotype dataframe
groupedby_df = pd.read_csv(phenotype_input_file, sep='\t')

#subset columns
if selected_subset is not None:
    phenotype_df = groupedby_df.query('{}==@selected_subset'.format(subset_field))
else:
    phenotype_df = groupedby_df.copy()

# should only be cells from one celltype
assert(len(phenotype_df['celltype'].drop_duplicates())==1)

#create merged index
phenotype_df['index'] = phenotype_df['donor_id']
cell_count_ds = phenotype_df['n_cells']
cell_count_ds.index = phenotype_df['index'].tolist()

#take donor and pool cols for samplemapping df
samplemapping_df = phenotype_df[['donor_id','index']]
samplemapping_df.to_csv(sample_mapping_file, sep='\t', index=False, header=False)

#reorganise to just be expression data indexed by the merged index
cols_to_drop = list(set(phenotype_df.columns) & {'donor_id','celltype','n_cells','treatment','pool_id'})
phenotype_df = phenotype_df.drop(cols_to_drop, axis=1).set_index('index')
phenotype_df = phenotype_df.transpose()

#number of top expressed genes
selected_genes = list(phenotype_df.mean(axis=1).nlargest(n_genes).index)
phenotype_df = phenotype_df.loc[selected_genes, :]

phenotype_list = phenotype_df.index
mapping_df = mapping_df.query('hgnc_symbol in @phenotype_list')
mapping_df = mapping_df.drop_duplicates(subset=['hgnc_symbol'])
mapping_df = mapping_df.set_index('hgnc_symbol')

# limit only to hgnc symbols that map to ensembl gene IDs
phenotype_df = phenotype_df.loc[mapping_df.index,:]
phenotype_df.index = mapping_df['ensembl_gene_id']

phenotype_df.to_csv(phenotype_file, sep='\t')

pc_mat = sklearn.decomposition.PCA(n_components=n_pcs).fit_transform(phenotype_df.values.transpose())
pc_df = pd.DataFrame(data=pc_mat, index=phenotype_df.columns, columns=['PC{}'.format(x) for x in range(1,n_pcs+1)])

pc_df.to_csv(covariates_file, sep='\t')



noise_scaling_vector = [1/float(x) for x in cell_count_ds.tolist()]
noise_matrix = np.diag(noise_scaling_vector)
noise_matrix_df = pd.DataFrame(data=noise_matrix, index=cell_count_ds.index, columns=cell_count_ds.index)

noise_matrix_df.to_csv(noise_matrix_file, sep='\t')
kinship_file = noise_matrix_file


number_of_permutations = '0'
minor_allele_frequency = '0.05'
hwe = '0.000001'
call_rate = '1'
window_size = '500000'
block_size = '15000'


config_dict = {'af': annotation_file,
               'pf': phenotype_file,
               'cf': covariates_file,
               'kf': kinship_file,
               'smf': sample_mapping_file,
               'plink': genotypes_file,
               'maf': minor_allele_frequency,
               'hwe': hwe,
               'np': number_of_permutations,
               'cr': call_rate,
               'w': window_size,
               'bs': block_size,
               'chunk_file': chunk_file
}

# write config to a file
qtl_config_utils.write_config(config_dict, config_file)
