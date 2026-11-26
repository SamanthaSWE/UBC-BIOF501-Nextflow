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

##	A Directed Acyclic Graph (DAG) and summary of the pipeline stages
from the root github directory, run:
- ```open pipeline_dag.html```
Preprocess expression: inputs only the original /Users/samanthaoxley/Git/UBC-BIOF501-Nextflow/data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv. Filters to 1-row per cell-line, dropped metadata, exclude incomplete genes.

Parallel Analysis processes: Principal Component Analysis (PCA)/Cosine Similarity and Spearman Correlation
PCA: Calculate 500 PCA components from all gene expression profiles and generate cosine similarity
Spearman: Input the Human Protein Atlas to filter genes to only known lineage-drivers, calc spearman correlation
Similarity Evaluation: visualization hub - generate UMAP for both pca cosine and spearman, also compare lineage agreement between them, and check unknown lineage for agreement

##	Rationale for the pipeline’s design
- Isolation of input processing allows for customization of filtering, and easily swapping input data. 
    - For example, remove depmap and insert other RNA expression dataset of choice. Can easily update the python code.
- Parallelization of PCA Cosine and Spearman Correlation: 
    - Can plug-and-play new analysis methods easily, the runtime only goes up as the longest analysis type (in this case: spearman)
- Final Visualization: 
    - lots of different visualization options can be included. 
    - This is the only place lineage mapping is included, as to not impact the similarity calculations.

##	An overview of core and pipeline dependencies
- Nextflow - manages each process into a specific workflow where outputs from one process can pass easily to the next process
- Docker - allows for python requirements to be installed within its own micro environment, and execution happens within environment as well
    - Python 3.13 image version, + custom requirements pinned to version for reproducibility

##	Step-by-step instructions for running the pipeline, including the installation guide for core dependencies
Install Prerequisites
1. [[Java 17 or later](https://docs.oracle.com/en/java/javase/27/install/overview-jdk-installation.html)] — required to run Nextflow
2. [[Nextflow](https://docs.seqera.io/nextflow/install)]— workflow execution engine
3. [[Docker Desktop](https://www.docker.com/)] — runs the pipeline's containerized Python scripts

Run Pipeline:
1. Clone the repository:
    - ```git clone UBC-BIOF501-Nextflow```
    - ```cd UBC-BIOF501-Nextflow```
2. Start Docker Desktop - [[Docker Guide](https://docs.docker.com/go/guides/)]
3. Run nextflow .nf
    - ```nextflow run expression-similarity.nf```

Compare Expected Results:
    - compare 4 outputs in /results/final_similarity/figures with /expected_results
    - NOTE: UPAM point LABELING is random and NOT EXPECTED to perfectly match, included for visual interpretation. The UMAP POINTS are reproducible and will NOT vary.

### Developer Mode
Rebuild docker with command 
docker build --no-cache -t depmap-dev .
Ensure depmap-dev docker build

nextflow run expression-similarity.nf -with-timeline timeline.html

##	An overview of the directory structure
**Root Directory** (UBC-BIOF501-Nextflow) contains nextfl\
ow workflow definition and config, dockerfile, python requirements, and README.md
**/data/** contains pre-requisite data (larger files pulled from git)
**/python/** contains all python scripts 
**/results/** store all intermediate and final results
**/expected-results/** store expected outputs. NOTE: UPAM point LABELING is random, for more dispersed lineages it is expected there is a different label of points. The UMAP points are reproducible and will NOT vary.