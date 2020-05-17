"""
Snakefile for single-cell dopaminergic neuron differentiation project for processing 10x samples

Author: Daniel Seaton
Affiliation: EMBL-EBI
Study: Single-cell dopaminergic neuron differentiation project
Date: Friday 11 May 2018
Run: snakemake -s Snakefile_10x.py --jobs 1000 --latency-wait 300 --cluster-config cluster.json --cluster 'bsub -q {cluster.queue} -J {cluster.name} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o {cluster.output} -e {cluster.error}' --use-conda

Based on snakefile by Davis McCarthy, EMBL-EBI

Raw data from Sanger first downloaded to:
../data/data_raw/

"""

import yaml
import glob
import os
from subprocess import run
import pandas as pd
import re

configfile: "config.yaml"
shell.prefix("set -euo pipefail;") 

permutation_mode = False
# turn temp file specifications on or off
def temp_context(infile):
    if permutation_mode:
        return temp(infile)
    else:
        return infile



## REFERENCE FILES
sample_metadata_file = config['sample_metadata_file']
pool_metadata_file = config['pool_metadata_file']
vcf_for_ase = ''

# metadata tables
sample_df = pd.read_csv(sample_metadata_file, sep='\t', dtype={'run':str})
sample_df = sample_df.set_index('sample_id', drop=False)
pool_df = pd.read_csv(pool_metadata_file, sep='\t', index_col=0)

# link sample names to cell lines
sample2lines_dict = dict()
for sample in list(sample_df.index):
    pool_id = sample_df.loc[sample, 'pool_id']
    line_list = pool_df.loc[pool_id, 'cell_line_id_long'].tolist()
    sample2lines_dict[sample] = line_list


# load config files. Note that this will silently drop files that are specified but that don't exist.
with open('./dataset_config_files.txt', 'r') as f:
    input_config_files = [x.strip() for x in f.readlines()]
config_files = [x for x in input_config_files if not x.startswith('#')]
config_files = [glob.glob(x) for x in config_files]
config_files = [y for x in config_files for y in x]
print('Running analysis for:')
for f in config_files:
    print(f)

dataset_config_dicts = dict()
# dataset names are the name of the yaml file
dataset_list = [os.path.basename(x).replace('.yaml','') for x in config_files]
# out prefix is everything up to the ".yaml" of the config file
out_prefix_list = [x.replace('.yaml','') for x in config_files]
# out dir is where the yaml file parent directory is
out_dir_list = [os.path.dirname(x) for x in config_files]
# load yaml files into dictionaries
for config_file in config_files:
    dataset = os.path.basename(config_file).replace('.yaml','')
    if os.path.basename(os.path.dirname(config_file)) != dataset:
        raise ValueError('yaml file and parent directory must have the same name, i.e.: {dir}/{dir}.yaml')
    with open(config_file, 'r') as f:
        config_dict = yaml.safe_load(f)
        dataset_config_dicts[dataset] = config_dict


raw_data_dir = '/hps/nobackup2/stegle/users/dseaton/hipsci/singlecell_neuroseq/data/data_raw/'


groupby_terms = ['donor_id-cluster_id','donor_id-cluster_id-time_point-treatment','donor_id-pool_id-cluster_id','pool_id','donor_id']
normalised_datasets = expand('{out_prefix}.scanpy.w_metadata.scanpy.h5', out_prefix=out_prefix_list)
normalised_dataset_obs_dfs = expand('{out_prefix}.scanpy.w_metadata.scanpy.obs_df.tsv', out_prefix=out_prefix_list)
merged_datasets = expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.h5', out_prefix=out_prefix_list)
merged_datasets += expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.obs_df.tsv', out_prefix=out_prefix_list)
merged_datasets += expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point.celltype_fractions_pivoted.tsv', out_prefix=out_prefix_list)
merged_datasets += expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point-treatment.celltype_fractions_pivoted.tsv', out_prefix=out_prefix_list)
#merged_datasets += expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.cluster_id.{selected_cluster}.groupedby.{groupby}.mean.tsv', out_prefix=out_prefix_list, selected_cluster=cluster_list, groupby=groupby_terms)


vardecomp_results = expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point.celltype_fractions_pivoted.var_decomp.tsv', out_prefix=out_prefix_list)



groupby_terms = ['donor_id-celltype','donor_id-celltype-treatment','donor_id-pool_id-celltype','donor_id-pool_id-celltype-treatment'][:]
celltype_list = ['DA','Sert','Epen1','Astro','P_FPP','FPP','NB'][-3:]
celltype_list = ['unknown']
merged_datasets += expand('{out_prefix}.scanpy.w_metadata.w_celltype.scanpy.celltype.{selected_celltype}.groupedby.{groupby}.mean.tsv', out_prefix=out_prefix_list, selected_celltype=celltype_list, groupby=groupby_terms)

#qtl_discovery_files = ['covariates.tsv', 'noise_matrix.tsv', 'phenotypes.tsv', 'qtl_config.yaml', 'samplemapping.tsv']
qtl_discovery_files = ['qtl_config.yaml']
cis_eqtl_discovery_files = expand('{out_dir}/qtl_analysis/eqtl_discovery/celltype_{celltype}/{file}', out_dir=out_dir_list, celltype=celltype_list, file=qtl_discovery_files)
cis_eqtl_sumstats_files = expand('{out_dir}/qtl_analysis/eqtl_discovery/celltype_{celltype}/summary_stats/qtl_config.yaml', out_dir=out_dir_list, celltype=celltype_list)
#cis_eqtl_discovery_files = expand('{out_dir}//qtl_analysis/eqtl_discovery/celltype_{celltype}/{treatment}/{file}', dataset=[x for x in dataset_list if x.endswith('D52')], celltype=celltype_list, file=qtl_discovery_files, treatment=['NONE','ROT'])
celltype_gwas_files = expand('{out_dir}/qtl_analysis/celltype_gwas/run_analysis.sh', out_dir=out_dir_list)


#ase_output = expand(raw_data_dir + '{sample}/possorted_genome_bam.cellsnp_ase.vcf.gz', sample=sample_list)
#ase_output = expand(raw_data_dir + 'ase/{sample}/cellsnp_ase.finished_dummy_file.txt', sample=sample_list[:120])


rule all:
    input:
#        ase_output,
        cis_eqtl_discovery_files,
#        cis_eqtl_sumstats_files,
#        celltype_gwas_files,
#        normalised_dataset_obs_dfs,
#        normalised_datasets,
#        merged_datasets

rule merge_scanpy:
    input:
        config_file = '{out_parent_dir}/{dataset}/{dataset}.yaml'
    output:
        output_file = temp_context('{out_parent_dir}/{dataset}/{dataset}.scanpy.h5')
    conda:
        'envs/scanpy_env.yaml'
    shell:
        '  python scripts/merge_scanpy.py --config_file {input.config_file} --output_file {output.output_file} '

rule pca_scanpy:
    input:
        input_file = '{out_parent_dir}/{dataset}/{dataset}.scanpy.h5'
    output:
        dataset_outfile = temp('{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.h5'),
        pca_outfile = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.PCA.tsv'
    conda:
        'envs/scanpy_env.yaml'
    shell:
        ' python scripts/dim_reduction_scanpy.py --input_file {input.input_file} '
        ' --dataset_outfile {output.dataset_outfile} --pca_outfile {output.pca_outfile} '

rule run_harmony:
    input:
        pca_file = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.PCA.tsv',
        metadata_file = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.obs_df.tsv'
    output:
        '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.tsv'
    params:
        harmony_batch_label  = lambda wildcards : dataset_config_dicts[wildcards.dataset]['harmony_batch_label'],
        # default theta to set if not described in config file
        harmony_theta = lambda wildcards: (dataset_config_dicts[wildcards.dataset]['harmony_theta']
                                           if 'harmony_theta' in dataset_config_dicts[wildcards.dataset].keys()
                                           else 2.0)
    shell:
#        '  source /nfs/software/stegle/system/Anaconda3-2018.12-Linux-x86_64/etc/profile.d/conda.sh ; '
#        '  source activate R_3p5 ; '
        '  /nfs/software/stegle/users/dseaton/conda-envs/R_3p5/bin/Rscript scripts/run_harmony.R '
        ' {input.pca_file} {output} '
        ' {input.metadata_file} '
        ' {params.harmony_batch_label} '
        ' {params.harmony_theta} '

 
rule cluster_cells_scanpy:
    input:
        h5_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.h5',
        pca_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.tsv'
    output:
        h5_output = temp('{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.clustered.h5'),
        cluster_output = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.cell_clustering_df.tsv'
    params:
        resolution  = lambda wildcards : dataset_config_dicts[wildcards.dataset]['resolution']
    conda:
        'envs/scanpy_env.yaml'
    shell:
        ' python scripts/cluster_cells_scanpy.py --h5_input {input.h5_input} --pca_input {input.pca_input} '
        ' --h5_output {output.h5_output} --cluster_output {output.cluster_output} --resolution {params.resolution} '


rule add_metadata_to_scanpy_h5:
    input:
        h5_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.h5',
        cluster_h5_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.clustered.h5',
        cluster_tsv_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.dimreduction.harmonyPCA.cell_clustering_df.tsv'
    output:
        h5_output = temp_context('{out_parent_dir}/{dataset}/{dataset}.scanpy.w_metadata.scanpy.h5')
    conda:
        'envs/scanpy_env.yaml'
    shell:
        ' python scripts/add_metadata_to_scanpy_h5.py --h5_input {input.h5_input} '
        ' --h5_output {output.h5_output} --cluster_tsv_input {input.cluster_tsv_input} '
        ' --cluster_h5_input {input.cluster_h5_input} '

rule add_celltype_to_scanpy_h5:
    input:
        h5_input = '{out_parent_dir}/{dataset}/{dataset}.scanpy.w_metadata.scanpy.h5',
        celltype_mapping_file = '{out_parent_dir}/{dataset}/{dataset}.celltype_mapping.tsv'
    output:
        h5_output = '{out_parent_dir}/{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.h5'
    conda:
        'envs/scanpy_env.yaml'
    shell:
        ' python scripts/add_celltype_to_scanpy_h5.py --h5_input {input.h5_input} '
        ' --h5_output {output.h5_output} --celltype_mapping {input.celltype_mapping_file} '


rule extract_obs_df:
    input:
        '{file_path}.h5'
    output:
        '{file_path}.obs_df.tsv'
    conda:
        'envs/scanpy_env.yaml'
    shell:
        '  python ./scripts/extract_adata_obs_dataframe.py '
        ' --input {input} '
        ' --output {output} '

rule count_celltypes:
    input:
        '{file_path}.obs_df.tsv'
    output:
        '{file_path}.obs_df.groupedby.{groupby}.celltype_fractions.tsv',
        '{file_path}.obs_df.groupedby.{groupby}.celltype_fractions_pivoted.tsv',
        '{file_path}.obs_df.groupedby.{groupby}.celltype_counts.tsv',
        '{file_path}.obs_df.groupedby.{groupby}.sample_counts.tsv'
    params:
        groupby_categories = lambda wildcards: wildcards.groupby.split('-')
    conda:
        'envs/scanpy_env.yaml'
    shell:
            "  python ./scripts/count_celltypes.py "
            " --input {input} "
            " --groupby {params.groupby_categories} "

rule create_pseudobulk_dataframes:
    input:
        exprs_file = '{file_path}.h5',
    output:
        '{file_path}.{subset_field}.{selected_subset}.groupedby.{groupby}.mean.tsv'
    params:
        output_prefix = '{file_path}',
        groupby_categories = lambda wildcards: wildcards.groupby.split('-')
    conda:
        'envs/scanpy_env.yaml'
    shell:
        " python ./scripts/create_pseudobulk_dataframes.py "
        " --input {input.exprs_file} "
        " --output_prefix {params.output_prefix} "
        " --groupby {params.groupby_categories} "
        " --function mean "
        " --selected_subset {wildcards.selected_subset} "
        " --subset_field {wildcards.subset_field} "


# rule summarise_donor_cluster_mapping:
#     input:
#         exprs_file = out_dir + '{sample}/expression_data.vargenefiltered.tsv',
#         metadata_file = out_dir + '{sample}/cell_metadata.tsv'
#     output:
#         fig_out_dir + '{sample}/heatmaps_of_cell_counts_across_cell_lines_and_clusters.png',
#         fig_out_dir + '{sample}/summary_of_cell_counts_across_cell_lines_and_clusters.png'
#     params:
#         out_dir = fig_out_dir
#     shell:
#         ' source activate py3 ; python ./scripts/summarise_donor_cluster_mapping.py '
#         ' {input.exprs_file} {input.metadata_file} {params.out_dir} '


# rule call_ase_cellsnp:
#     input:
#         barcodes=raw_data_dir + '{sample}/filtered_gene_bc_matrices/hg19/barcodes.tsv',
#         vcf=vcf_for_ase,
#         bam=raw_data_dir + '{sample}/possorted_genome_bam.bam'
#     output:
#         raw_data_dir + '{sample}/possorted_genome_bam.cellsnp_ase.vcf.gz'
#     conda:
#         'envs/cellsnp_env.yaml'
#     shell:
#         'cellSNP -s {input.bam} -o {output} '
#         ' -b {input.barcodes} '
#         ' -R {input.vcf} -p 20 '

# rule split_ase_by_donor:
#     input:
#         vcf = raw_data_dir + '{sample}/possorted_genome_bam.cellsnp_ase.vcf.gz',
#         mapping_file = out_dir + '{sample}/demuxlet_out.cell_to_donor_mapping.tsv'
#     output:
#         raw_data_dir + '{sample}/ase/cellsnp_ase.finished_dummy_file.txt'
#     params:
#         out_prefix = raw_data_dir + '{sample}/ase/cellsnp_ase'
#     conda:
#         'envs/cyvcf2_env.yaml'
#     shell:
#         'python /nfs/software/stegle/users/dseaton/dev/corrall/split_by_donor.py --out_file {output} --out_prefix {params.out_prefix} --vcf_file {input.vcf} --mapping_file {input.mapping_file}'


rule prepare_cis_eqtl_discovery_files:
    input:
        '{file_path}/{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.celltype.{celltype}.groupedby.donor_id-celltype.mean.tsv'
    output:
        '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype, [A-Za-z0-9\_]+}/qtl_config.yaml'
    params:
         out_dir = '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype}/'
    shell:
        ' python ./scripts/prepare_qtl_discovery_files.py {input} {params.out_dir} '


rule convert_qtl_config_for_sumstats_generation:
    input:
        '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype}/qtl_config.yaml'
    output:
        '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype, [A-Za-z0-9\_]+}/summary_stats/qtl_config.yaml'
    shell:
        ' python ./scripts/convert_qtl_config_for_sumstats_generation.py {input} {output} '


# rule prepare_cis_eqtl_discovery_files_treatment:
#     input:
#         os.path.join(out_dir, '{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.celltype.{celltype}.groupedby.donor_id-celltype-treatment.mean.tsv')
#     output:
#         '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype}/{treatment}/qtl_config.yaml')
#     params:
#         out_dir = '{file_path}/{dataset}/qtl_analysis/eqtl_discovery/celltype_{celltype}/{treatment}/')
#     shell:
#         ' python ./scripts/prepare_qtl_pipeline_files.py {input} {params.out_dir} treatment {wildcards.treatment}'

rule prepare_celltype_gwas_files:
    input:
        phenotype_file = '{file_path}/{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point.celltype_fractions_pivoted.tsv',
        sample_counts_file = '{file_path}/{dataset}/{dataset}.scanpy.w_metadata.w_celltype.scanpy.obs_df.groupedby.donor_id-pool_id-time_point.sample_counts.tsv'
    output:
        '{file_path}/{dataset}/qtl_analysis/celltype_gwas/run_analysis.sh'
    params:
        out_dir = '{file_path}/{dataset}/qtl_analysis/celltype_gwas/'
    shell:
        ' python ./scripts/prepare_celltype_gwas_files.py '
        ' --phenotype_file {input.phenotype_file} '
        ' --sample_counts_file {input.sample_counts_file} '
        ' --selected_variant_file /nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/qtl_analysis/eqtl_discovery_all_leads.fdr0.01.tsv '
        ' --out_dir {params.out_dir} '


rule cellfraction_variance_decomposition:
    input:
        infile = '{file_path}.obs_df.groupedby.{groupby}.celltype_fractions_pivoted.tsv',
    output:
        '{file_path}.obs_df.groupedby.{groupby}.celltype_fractions_pivoted.var_decomp.tsv',
    shell:
        ' PS1=${{PS1:-}} ; source /nfs/software/stegle/system/Anaconda3-2018.12-Linux-x86_64/etc/profile.d/conda.sh ; '
        ' conda activate limix_env ; '
        ' python ./scripts/cellfraction_variance_decomposition.py '
        ' --infile {input.infile} '
        ' --outfile {output} '
        ' --transform none '
        ' --effects donor_id pool_id '
