#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(harmony)
library(dplyr)

pca_file = args[1]
out_file = args[2]
metadata_file = args[3]
merge_column = args[4]
theta = as.numeric(args[5])

#default theta is 2

#pca_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/pool1_pool2_D52/pool1_pool2_D52.scanpy.dimreduction.PCA.tsv'
#out_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/pool1_pool2_D52/pool1_pool2_D52.scanpy.dimreduction.harmonyPCA.tsv'
#metadata_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/pool1_pool2_D52/pool1_pool2_D52.scanpy.dimreduction.obs_df.tsv'


sample_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/metadata/sample_tracking.tsv'

pca_df = read.table(pca_file, sep='\t', header=TRUE, row.names=1)

metadata_df = read.table(metadata_file, sep='\t', header=TRUE)

sample_df = read.table(sample_file, sep='\t', header=TRUE)
sample_df$sample_id = sample_df$sanger_sample_id
sample_df = sample_df[,c('sample_id','pool_id','time_point','treatment')]

metadata_df = left_join(metadata_df, sample_df, by = "sample_id")
rownames(metadata_df) = metadata_df$index

head(metadata_df, 10)


harmony_embeddings <- HarmonyMatrix(pca_df, metadata_df, merge_column, theta=theta)

write.table(harmony_embeddings, out_file, sep='\t', quote=FALSE)