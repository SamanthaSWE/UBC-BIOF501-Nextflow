# UBC-BIOF501-Nextflow
A basic nextflow pipeline evaluating expression data

## Project background
In an oncology clinic, a patient's tumor may be evaluated by transcriptome via RNA sequencing. A tumor's transcriptome provides insight to unique patterns of gene expression and disregulation, and may enable rationalizing non-standard targeted treatments. Characterizing the similarity of a patient's tumor transcriptome against experimental cell line transcriptomes may help to assess the relevance of experimental findings to that patient's cancer.

Cancer Dependency Map (DepMap) is a large collection of molecular and experimental data primarily focused on human cancer cell lines. Like all cancers in patients, each cancer cell line has it's own molecular characterists. Many of these cell lines have been used in experiments including drug sensitivity screens, CRISPR knockout, and RNAi knockdown studies. These experiments have the potential to inform personalized cancer treatment. However, if an experimental cell line differs substantially from a patient's tumor, the insights from experiments with that cell line may be less applicable to that patient's cancer.

In this project, I do not have access to clinical tumor data. Instead, I will evaluate the transcriptomic similarity between cell lines in DepMap. These computational methods that can serve as a foundation for additional development.


Caveats: a few major **simplifications** were taken in the pipeline:
1. Processing clinical tumor samples requires accounting for differences in expression sequencing platforms, sample purity, and cellular heterogeneity.
2. We only evaluate the accuracy of calculated outputs against known cell lineage. **However, calculating a sample's lineage similarity to cell lines doesn't mean a sample will respond the same to cell lines tested treatments.** 
    - To test treatment-relevant similarity analysis, additional research is needed. One strategy could modify the existing similarity calculations by first filtering to pathway-relevant RNA expression. Also, datasets including noncoding RNA may be beneficial.

##	Further explanation of DepMap data:
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv describes the expression of 19,218 protein-coding genes across 2,130 cell lines. This dataset has already been pre-processed by DepMap to remove duplicates and QC failures, rename genes, perform TPM normalization, and log transforms raw values.

As stated in the project background, this dataset is very suitable as many of these cell lines have been used in experiments including drug sensitivity screens, CRISPR knockout, and RNAi knockdown studies. The insights of these studies can be made relevant to a patient's tumor sample through similarity analysis.

##	Rationale for the pipeline’s design
TODO: with each step and main tool explained and justified
##	A Directed Acyclic Graph (DAG) and summary of the pipeline stages
flowchart TB
    subgraph " "
    v0["channel.fromPath"]
    v1["channel.fromPath"]
    v2["channel.fromPath"]
    v3["python_script"]
    v6["python_script"]
    v8["python_script"]
    v10["python_script"]
    end
    v4(["PreprocessExpression"])
    subgraph " "
    v5[" "]
    v12[" "]
    end
    v7(["PCAAnalysis"])
    v9(["SpearmanAnalysis"])
    v11(["SimilarityEvaluation"])
    v0 --> v4
    v1 --> v9
    v2 --> v11
    v3 --> v4
    v4 --> v7
    v4 --> v5
    v4 --> v9
    v6 --> v7
    v7 --> v11
    v8 --> v9
    v9 --> v11
    v10 --> v11
    v11 --> v12

##	An overview of core (e.g. Nextflow, conda) and pipeline dependencies (e.g. FastQC, STAR, DESeq2)

##	Step-by-step instructions for running the pipeline, including the installation guide for core dependencies
Install Prerequisites
1. [[Java 17 or later](https://docs.oracle.com/en/java/javase/27/install/overview-jdk-installation.html)] — required to run Nextflow
2. [[Nextflow](https://docs.seqera.io/nextflow/install)]— workflow execution engine
3. [[Docker Desktop](https://www.docker.com/)] — runs the pipeline's containerized Python scripts

Run Pipeline
1. Clone the repository:
    - ```git clone UBC-BIOF501-Nextflow```
    - ```cd UBC-BIOF501-Nextflow```
2. Start Docker Desktop - [[Docker Guide](https://docs.docker.com/go/guides/)]
3. Run nextflow .nf
    - ```nextflow run expression-similarity.nf```

### Developer Mode
Rebuild docker with command 
docker build --no-cache -t depmap-dev .
Ensure depmap-dev docker build

nextflow run expression-similarity.nf -with-timeline timeline.html

##	An overview of the directory structure
Root Directory (UBC-BIOF501-Nextflow) contains nextfl\
ow workflow definition and config, dockerfile, python requirements, and README.md
/data/ contains pre-requisite data (larger files pulled from git)
/python/ contains all python scripts 
/results/ store all intermediate and final results, including:
    - 

# TODO!!
- replace docker file with image pointing to git
- check timeline
- create DAG