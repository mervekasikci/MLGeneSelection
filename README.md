# Machine Learning Based Gene Selection Algorithms for RNA-seq Gene Expression Data
In this GitHub repo, we provide biosigner algorithm [1], GMDH algorithm [2], Determan's optimal gene selection algorithm within support vector machines, random forest and elastic net generalized linear model [3] to perform gene selection for RNA-seq gene expression data. There are three folders in this repo named dataset, functions and gene selection. Descriptions of these folders and the files in them are below:

**Dataset:** We provide an example gene expression dataset which is available at https://portal.gdc.cancer.gov/projects/TCGA-KICH. Researcher use this example dataset or their own gene expression data.
**Functions:** We include the source codes required for the gene selection to this folder. 
Additionally, we add a source R code called "data_simulating.R" for researchers who want to generate their data via simulation. This code is a manipulated version of the sim.counts function in ssizeRNA R package [4].
**Geneselection:** We add codes for gene selection via machine learning algorithms to this folder. The 'geneselection_biosigner.R' code is used for gene selection with the biosigner algorithm, 'geneselection_gmdh.R' for the GMDH algorithm, and 'geneselection_omicsmarker.R' for Determan's algorithm within support vector machines, random forest and elastic net generalized linear models. The files in these folders not only enable gene selection with machine learning algorithms, but also enable the preprocessing steps to be applied to the data. The codes in this folder contain the following steps:

1. Data uploading: Researchers use example dataset (KICH_data.txt) or they can upload their raw RNA-seq gene expression dataset.
2. Filtering: Near-zero variance filtering is applied to remove genes with low variance.
3. Normalization: Median ratio normalization is performed to the filtered count data.
4. Transformation: Logarithmic transformation is applied to the normalized data to obtain a less skewed distribution.
5. Univariate analysis: Student's t-test is performed to reduce the number of genes and identify genes that differed significantly between two groups. As default, according to the results of the Student's t-test, genes are ranked from lowest to highest p-value, and the first 200 genes are selected.
6. Gene selection via machine learning algorithm: The machine learning model is developed according to the preferred algorithm among the biosigner, GMDH or Determan's optimal gene selection algorithms. The recommended gene list is obtained. Finally, the recommended gene list based on algorithm is obtained and can be saved in csv format.

References: 
[1] Rinaudo P, Boudah S, Junot C, Thevenot E (2016) biosigner: a new method for the discovery of significant molecular signatures from omics data. Frontiers in Molecular Biosciences 3:26.
[2] Dag O, Karabulut E, Alpar R (2019) Gmdh2: Binary classification via gmdh-type neural network algorithms - r package and web-based tool. International Journal of Computational Intelligence Systems 12(2):649–660.
[3] Determan C (2015) Optimal algorithm for metabolomics classification and feature selection varies by dataset. International Journal of Biology 7(1):100–115.
[4] Bi R, Liu P (2019, version 1.3.2) ssizeRNA: Sample Size Calculation for RNA-Seq Experimental Design, R Package. https://cran.r-project.org/web/packages/ssizeRNA/index.html.
