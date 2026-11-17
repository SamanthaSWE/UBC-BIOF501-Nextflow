# UBC-BIOF501-Nextflow
A basic nextflow pipeline evaluating expression data

## Project background
In an oncology clinic, a patient's tumor may be evaluated by transcriptome via RNA sequencing. A tumor's transcriptome provides insight to unique patterns of gene expression and disregulation, and may enable rationalizing non-standard targeted treatments. Characterizing the similarity of a patient's tumor transcriptome and experimental cell line transcriptome helps to assess the relevance of experimental findings to that patient's cancer.

Cancer Dependency Map (DepMap) is a large collection of molecular and experimental data primarily focused on human cancer cell lines. Like all cancers in patients, each cancer cell line has it's own molecular characterists. Many of these cell lines have been used in experiments including drug sensitivity screens, CRISPR knockout, and RNAi knockdown studies. These experiments have the potential to inform personalized cancer treatment. However, if an experimental cell line differs substantially from a patient's tumor, the insights from experiments with that cell line may be less applicable to that patient's cancer.

In this project, I do not have access to clinical tumor data. Instead, I will evaluate the transcriptomic similarity between cell lines using established computational methods that can serve as a foundation for future projects. Extending these methods to clinical tumor samples would require additional consideration, including accounting for differences in expression sequencing platforms, sample purity, and cellular heterogeneity.

##	An explanation of the data used

I have downloaded OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv which describes the expression of 19,218 genes across 2,130 cell lines. This dataset has already been pre-processed by DepMap to remove duplicates and QC failures, rename genes, TPM normalization, and log transforms raw values.

##	Rationale for the pipeline’s design

##	A Directed Acyclic Graph (DAG) and summary of the pipeline stages
##	An overview of core (e.g. Nextflow, conda) and pipeline dependencies (e.g. FastQC, STAR, DESeq2)
##	Step-by-step instructions for running the pipeline, including the installation guide for core dependencies
##	An overview of the directory structure
