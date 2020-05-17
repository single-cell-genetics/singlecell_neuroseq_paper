import sys
import pandas as pd
import os
import shutil
import sklearn.decomposition
import qtl_config_utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--phenotype_file',
                    help="Input phenotype file", type=str)
parser.add_argument('--sample_counts_file',
                    help="Input sample counts file", type=str)
parser.add_argument('--to_sum',
                    help="Celltypes to be added together. Comma separated list of hyphen separated lists (e.g. 'DA-Sert,P_FPP-FPP')",
                    type=str, default=None)
parser.add_argument('--selected_variant_file',
                    help="Input selected variant file", type=str)
parser.add_argument('--mode',
                    help="Analysis mode - celltype_fractions or celltype_pcs", type=str, default='celltype_fractions')
parser.add_argument('--out_dir',
                    help="Output directory", type=str)
args = parser.parse_args()

phenotype_input_file = args.phenotype_file
sample_counts_input_file = args.sample_counts_file
qtl_file = args.selected_variant_file
analysis_type = args.mode
to_sum = args.to_sum
outdir = os.path.abspath(args.out_dir)
# add slash so that hipsci qtl pipeline postprocessing works properly
if not outdir.endswith('/'):
    outdir += '/'

# general set up
#genotypes_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered_Plink-f/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20170327.genotypes.norm.renamed.recode.filtered'
genotypes_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Plink/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.norm.renamed.recode.vcf.gz'
kinship_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered_Plink-f/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20170327.genotypes.norm.renamed.recode.filtered.rel'

base_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/'


qtl_df = pd.read_csv(qtl_file, sep='\t')
variant_df = qtl_df[['snp_id']]



annotation_file = os.path.join(outdir, 'annotation.tsv')
samplemapping_file = os.path.join(outdir, 'samplemapping.tsv')
phenotype_file = os.path.join(outdir, 'phenotypes.tsv')
covariates_file = os.path.join(outdir, 'covariates.tsv')
variantfilter_file = os.path.join(outdir, 'variant_filter.txt')
config_file = os.path.join(outdir, 'qtl_config.yaml')

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# write variant filter
variant_df.to_csv(variantfilter_file, sep='\t', index=False)

# process phenotype dataframe
groupedby_df = pd.read_csv(phenotype_input_file, sep='\t')
sample_cellcounts_df = pd.read_csv(sample_counts_input_file, sep='\t')


#subset
celltype_df = groupedby_df.iloc[:, :]

#create merged index
sample_cellcounts_df['index'] = sample_cellcounts_df['donor_id'] + '_' + sample_cellcounts_df['pool_id']
selected_samples = sample_cellcounts_df.query('n_cells>15')['index'].tolist()

celltype_df['index'] = celltype_df['donor_id'] + '_' + celltype_df['pool_id']
celltype_df = celltype_df.query('index in @selected_samples')


experimentmapping_df = celltype_df[['index','pool_id']].set_index('index')



#take donor and pool cols for samplemapping df
#samplemapping_df = celltype_df[['index','donor_id']]
samplemapping_df = celltype_df[['donor_id','index']]
samplemapping_df.to_csv(samplemapping_file, sep='\t', index=False, header=False)

#reorganise to just be expression data indexed by the merged index
celltype_df = celltype_df.drop(['donor_id','time_point','pool_id'], axis=1).set_index('index')

#represented cell types
selected_celltypes = celltype_df.mean(axis=0)>0.01
celltype_df = celltype_df.loc[:, selected_celltypes]

if analysis_type == 'celltype_fractions':
    phenotype_df = celltype_df.transpose()
    if to_sum is not None:
        to_sum = to_sum.split(',')
        for group in to_sum:
            indices = group.split('-')
            if not all([x in phenotype_df.index for x in indices]):
                print('Not able to apply sum of celltypes: {}'.format(group))
                continue
            phenotype_df.loc[group] = phenotype_df.loc[indices, :].sum(axis=0)

elif analysis_type == 'celltype_pcs':
    if to_sum is not None:
        raise(ValueError)
    n_components = 5
#    df = celltype_df.apply(lambda x: (x-x.mean())/x.std(), axis=0)
    model = sklearn.decomposition.PCA(n_components=n_components).fit(celltype_df.values.T)
    print(model.explained_variance_ratio_)
    pc_list = ['PC{}'.format(x) for x in range(1,n_components+1)]
    pca_df = pd.DataFrame(data=model.components_, index=pc_list, columns=celltype_df.index).transpose()
    phenotype_df = pca_df.transpose()
else:
    raise(ValueError)
    
phenotype_list = phenotype_df.index
phenotype_df.to_csv(phenotype_file, sep='\t')

### annotation file
feature_list = list(phenotype_df.index)
annotation_df = pd.DataFrame(columns=['feature_id','chromosome','start','end'])
annotation_df['feature_id'] = feature_list
annotation_df['chromosome'] = 1
annotation_df['start'] = 0
annotation_df['end'] = 0

annotation_df.to_csv(annotation_file, sep='\t', index=False)

# experiment details as covariates
dummy_df = pd.get_dummies(experimentmapping_df['pool_id'], drop_first=True)

dummy_df.to_csv(covariates_file, sep='\t')


num_perm = '1000'
minor_allele_frequency = '0.05'
hwe = '0.001'
call_rate = '1'
window_size = '0'
norm_method = 'standardize'

config_dict = {'af': annotation_file,
               'pf': phenotype_file,
               'cf': covariates_file,
               'kf': kinship_file,
               'smf': samplemapping_file,
               'plink': genotypes_file,
               'maf': minor_allele_frequency,
               'hwe': hwe,
               'np': num_perm,
               'cr': call_rate,
               'w': window_size
}

# write config to a file
qtl_config_utils.write_config(config_dict, config_file)



mjbonder_script = '/hps/nobackup/stegle/users/mjbonder/tools/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py'
dev_script = '/nfs/software/stegle/users/dseaton/dev/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py'

script = dev_script
script = mjbonder_script

test_shell_command = (
    " /nfs/software/stegle/users/mjbonder/conda-envs/limix_qtl/bin/python {script} "
    " --plink {gen} "
    " -af {af} "
    " -pf {pf} "
    " -od {od} "
    " --kinship_file {kf} "
    " --sample_mapping_file {smf} "
    " -vf {vf} "  # -vf for variant filter
#       " -gr {chunkFull} "
    " -np {np} "
#    " -wp "    #### write permutations
    " -maf {maf} "
    " -hwe {hwe} "
    " -t " #trans analysis
    " -gm {norm_method} "
    " -w 0 " # set window for SNPs around the gene
    " -cf {cf} "
    ).format(maf=minor_allele_frequency, hwe=hwe, np=num_perm, kf=kinship_file, smf=samplemapping_file,
             od=outdir, af=annotation_file, pf=phenotype_file, vf=variantfilter_file,
             gen=genotypes_file, script=script, norm_method=norm_method, cf=covariates_file)

merge_shell_command = (
    "/nfs/software/stegle/users/mjbonder/conda-envs/limix_qtl/bin/python /hps/nobackup/stegle/users/mjbonder/tools/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py  "
    " -id {} "
    " -od {}  -oo "
    "  ".format(outdir, outdir))

merge_shell_command_2 = (
    "/nfs/software/stegle/users/mjbonder/conda-envs/limix_qtl/bin/python /hps/nobackup/stegle/users/mjbonder/tools/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py  "
    " -id {} "
    " -od {}  -oo "
    " -sfo -tfb ".format(outdir, outdir))

print_command = (" head {}/qtl_results_all.txt | cut -f 1,2,3,4,5,6  ".format(outdir))

command = test_shell_command + ' ; ' + merge_shell_command + ' ; ' \
    + merge_shell_command_2 + ' ; ' + print_command



print()
print(command)
print()

with open(outdir+'/run_analysis.sh', 'w') as f:
    f.write(command)
