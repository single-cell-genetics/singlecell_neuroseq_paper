# Population-scale single-cell RNA-seq profiling across dopaminergic neuron differentiation

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq and bulk RNAseq data for our paper:

[Population-scale single-cell RNA-seq profiling across dopaminergic neuron differentiation](https://www.nature.com/articles/s41588-021-00801-6)

## Analysis scripts

The following folders contain scripts for data processing and analysis.
A short description can be found below:

* [Plotting Notebooks](../master/plotting_notebooks/) contains all jupyter notebooks to reproduce the individual figures (main and supplements).

* [10x Analysis Pipeline](../master/10x_analysis_pipeline/) contains a snakemake pipeline to process the individual scRNA-seq files generated by the 10x CellRanger software.

* [Neuronal differentiation outcome predictor](../master/differentiation_prediction_model/) contains a jupyter notebook and example data to predict differentiation outcome of iPSC lines based on iPSC gene expression.




