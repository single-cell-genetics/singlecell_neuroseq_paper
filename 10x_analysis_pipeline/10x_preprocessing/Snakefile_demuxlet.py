"""
Snakefile for single-cell dopaminergic neuron differentiation project for processing 10x samples

Author: Daniel Seaton
Affiliation: EMBL-EBI
Study: Single-cell dopaminergic neuron differentiation project
Date: Friday 11 May 2018
Run:
snakemake -s Snakefile_demuxlet.py --jobs 1000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -q {cluster.queue} -J {cluster.name} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o {cluster.output} -e {cluster.error}' --use-conda

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

## REFERENCE FILES
sample_metadata_file = config['sample_metadata_file']
pool_metadata_file = config['pool_metadata_file']

# metadata tables
sample_df = pd.read_csv(sample_metadata_file, sep='\t', dtype={'run':str})
sample_df = sample_df.set_index('sample_id', drop=False)
pool_df = pd.read_csv(pool_metadata_file, sep='\t', index_col=0)

pool_line_list_template = '../../data/metadata/pool2line_mapping_files/{pool_id}.tsv'

# link sample names to cell lines
sample2lines_dict = dict()
for sample in list(sample_df.index):
    pool_id = sample_df.loc[sample, 'pool_id']
    line_list = pool_df.loc[pool_id, 'cell_line_id_long'].tolist()
    sample2lines_dict[sample] = line_list

## variant files
neuroseq_samples_small_vcf = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/genotypes/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.allchr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.vcf' # complete version, with all donors that may be used in any experiment


## define commands
python_cmd = '/nfs/software/stegle/users/davis/conda-envs/py3/bin/python'
star_cmd = '/nfs/software/stegle/bin/STAR'
#gatk_cmd = 'java -jar /nfs/software/stegle/users/davis/GenomeAnalysisTK.jar'
gatk_cmd = '/nfs/software/stegle/users/dseaton/java/jdk1.8.0_112/bin/java -Xmx48g -Xms24g -Djava.io.tmpdir=/hps/nobackup/hipsci/scratch/singlecell_endodiff/tmp -jar /hps/nobackup/stegle/users/mjbonder/tools/GATK/GenomeAnalysisTK.jar'
picard_cmd='/nfs/software/stegle/users/dseaton/java/jdk1.8.0_112/bin/java -Xmx48g -Xms24g -Djava.io.tmpdir=/hps/nobackup/hipsci/scratch/singlecell_endodiff/tmp -jar /nfs/software/stegle/users/dseaton/picard/picard.jar'
#Rscript_cmd='/nfs/software/stegle/users/davis/conda-envs/py3/bin/Rscript'
Rscript_cmd = 'singularity exec /hps/nobackup/hipsci/scratch/biocplus.img Rscript'
R_cmd='/nfs/software/stegle/users/davis/conda-envs/py3/bin/R'
vcfsort_cmd='perl /nfs/software/stegle/users/dseaton/vcftools/src/perl/vcf-sort'


out_dir = '../../data/data_processed/'
fig_out_dir = '../figures/'

raw_data_dir = '/hps/nobackup2/stegle/users/dseaton/hipsci/singlecell_neuroseq/data/data_raw/'
processed_data_dir = '../../data/data_processed/'


## parameter objects and samples
with open('../list_of_excluded_neuroseq_10x_samples.txt', 'r') as f:
    excluded_sample_list = [x.strip() for x in f.readlines()]
#excluded_sample_list.append('cellranger211_count_27615_5245STDY7685554_hg19-1_2_0')
sample_list = [os.path.basename(x) for x in glob.glob(raw_data_dir + 'cellranger211*')]
sample_list = [x for x in sample_list if x not in excluded_sample_list]
sample_list = sorted(sample_list)[:]
#sample_list = sample_df.query('pool_id=="pool10" or pool_id=="pool11"')['sample_id'].tolist()

#additional_samples =  ['cellranger211_count_25789_5245STDY7458661_hg19-1_2_0', 'cellranger211_count_25790_5245STDY7458663_hg19-1_2_0']

#samples_to_run_donor_search = ['cellranger211_count_26746_5245STDY7619068_hg19-1_2_0','cellranger211_count_26746_5245STDY7619069_hg19-1_2_0']
samples_to_run_donor_search = ['cellranger211_count_31186_5245STDY8119869_hg19-1_2_0','cellranger211_count_31186_5245STDY8119870_hg19-1_2_0',
                               'cellranger211_count_31186_5245STDY8119871_hg19-1_2_0','cellranger211_count_31186_5245STDY8119872_hg19-1_2_0']

line_subsets = ['line_subset_{}'.format(x) for x in range(13)]
demuxlet_search_output = expand(out_dir + 'demuxlet/{sample}/demuxlet_line_search_results/demuxlet_out.{line_subset}.best', 
                                sample=samples_to_run_donor_search, line_subset=line_subsets)

#donor_id_demuxlet_output = expand(out_dir + '{sample}/demuxlet_out.best', sample=sample_list)
donor_id_demuxlet_output = expand(out_dir + 'demuxlet/{sample}/demuxlet_out.cell_to_donor_mapping.tsv', sample=sample_list)

rule all:
    input:
        donor_id_demuxlet_output, demuxlet_search_output

rule donor_id_demuxlet:
    input:
        bam=raw_data_dir+'{sample}/possorted_genome_bam.bam',
        barcodes=raw_data_dir+'{sample}/filtered_gene_bc_matrices/hg19/barcodes.tsv'
    output:
        out_dir + 'demuxlet/{sample}/demuxlet_out.best'
    params:
        outfile_prefix = out_dir + 'demuxlet/{sample}/demuxlet_out',
        line_list = lambda wildcards : pool_line_list_template.format(pool_id=sample_df.loc[wildcards.sample,'pool_id'])
    shell:
        ' /nfs/software/stegle/users/huangh/bin/bin/demuxlet --sam {input.bam} '
        ' --vcf {neuroseq_samples_small_vcf} --doublet-prior 0.05 '
        ' --group-list {input.barcodes} --out {params.outfile_prefix} '
        ' --sm-list {params.line_list} '

rule extract_cell_donor_mapping:
    input:
        out_dir + 'demuxlet/{sample}/demuxlet_out.best'
    output:
        out_dir + 'demuxlet/{sample}/demuxlet_out.cell_to_donor_mapping.tsv'
    shell:
        ' python ./scripts/extract_cell_donor_mapping.py {input} {output} '

rule donor_id_demuxlet_search:
    input:
        bam=raw_data_dir + '{sample}/possorted_genome_bam.bam',
        barcodes=raw_data_dir + '{sample}/filtered_gene_bc_matrices/hg19/barcodes.tsv',
        line_list = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/metadata/search_all_hipsci_lines_mapping_files/{line_subset}.txt'
    output:
        out_dir + 'demuxlet/{sample}/demuxlet_line_search_results/demuxlet_out.{line_subset}.best'
    params:
        outfile_prefix = out_dir + 'demuxlet/{sample}/demuxlet_line_search_results/demuxlet_out.{line_subset}'
    shell:
        ' /nfs/software/stegle/users/huangh/bin/bin/demuxlet --sam {input.bam} '
        ' --vcf {neuroseq_samples_small_vcf} --doublet-prior 0.05 '
        ' --group-list {input.barcodes} --out {params.outfile_prefix} '
        ' --sm-list {input.line_list} '
