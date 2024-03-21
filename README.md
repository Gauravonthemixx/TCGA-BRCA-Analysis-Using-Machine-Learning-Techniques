# :rocket: Comprehensive Study of Clinical and Gene Expression TCGA Breast Cancer Data Using Machine Learning Techniques

## 📝 Description:
* We use Kaplan-Meier survival curves to investigate the effect of pathological tumor stage and the age of the patients on the overall survival rates of breast cancer.

* Predicting cancer staging in the patient using TCGA BRCA clinical and gene expression dataset using different models like Random Forest Survival, ANN and K-Nearest Neighbour to find out the most effective model
* We have used unsupervised learning methods like clustering on gene expression dataset to find potential patterns in the genomic data. 

* Clusters are used to improve supervised learning models and increase their predictive capabilities.

## ⌛ Dataset:

* To retrieve the data, we have used a R package called TCGAbiolinks. Package is freely and easily accessible in the Bioconductor Project at http://bioconductor.org/packages/TCGAbiolinks.

* Clinical and gene expression data is prepared using GDC query using TCGAbiolinks package in BCR.Biotab format.

### Preparing Clinical Supplement Datset in R:
```
# Extracting clinical data
gdcProject<- getGDCprojects()
TCGAbiolinks::getProjectSummary("TCGA-BRCA")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
) 

query_output<- getResults(query) 
GDCdownload(query)

clinical.BCRtab.all <- GDCprepare(query)

```

### Preparing Gene Expression Dataset in R:
```
# Extracting Gene Expression Data

GDCproject<-getGDCprojects()
TCGAbiolinks::getProjectSummary("TCGA-BRCA")

  query_TCGA <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    barcode = c("TCGA-*")) # 
    parameter enforced by GDCquery
GDCdownload(query_TCGA, method = "a
pi", files.per.chunk = 100)

getResults(query_TCGA)
tcga_gene<- GDCprepare(query_TCGA, s
ummarizedExperiment = TRUE)
gene_matrix<- assay(tcga_gene, 'unstranded')
```


