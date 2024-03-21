# :rocket: Comprehensive Study of Clinical and Gene Expression TCGA Breast Cancer Data Using Machine Learning Techniques

![hm1](https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/c384dc8b-ddfc-488c-b767-4a14290a5e7a)
![top (1)](https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/c69d9e5a-af28-4680-afcf-c06587cc7a31)
![kmeans(fix)](https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/ba123d93-98a5-490a-9bb9-fad3767d03b0)
<img width="708" alt="surv1" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/db5dc278-a01c-481d-9efa-6016dc257dbc">


## 📝 Description:
* We use Kaplan-Meier survival curves to investigate the effect of pathological tumor stage and the age of the patients on the overall survival rates of breast cancer.

* Predicting cancer staging in the patient using TCGA BRCA clinical and gene expression dataset using different models like Random Forest Survival, ANN and K-Nearest Neighbour to find out the most effective model
* We have used unsupervised learning methods like clustering on gene expression dataset to find potential patterns in the genomic data. 

* Clusters are used to improve supervised learning models and increase their predictive capabilities.

## ⌛ Dataset:

* To retrieve the data, we have used a R package called TCGAbiolinks. Package is freely and easily accessible in the Bioconductor Project at http://bioconductor.org/packages/TCGAbiolinks.

* Clinical and gene expression data is prepared using GDC query using TCGAbiolinks package in BCR.Biotab format.

### :pushpin: Preparing Clinical Supplement Datset in R:
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

### :pushpin: Preparing Gene Expression Dataset in R:
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

## 🖥️: Exploration of important variables in clinical dataset of breast cancer patients
<img width="617" alt="Screenshot 2024-03-14 at 21 41 07" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/effd74a5-69cf-4610-b2d4-336df4e59224">

<img width="687" alt="age(fix2)" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/6b7ea091-b684-41f2-979d-8fd9274e82ba">

* Patients around the age of 60 come to know about the tumors in their breast around this age in their life.

<img width="688" alt="race(fix)" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/cc7393bc-67ab-4479-a00b-d155135a9a50">

* The boxplot shows an analysis of the distribution of the current age across three racial groups; White, Black or African American and Asian.

* Invasive and Non Invasive cancer categorisation is achieved using analysis of Tumor Size, Node Stage and Metastasis Stage

<img width="592" alt="tumor(fix)" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/5e714c2b-b8d6-4d7a-a076-2df6ea8662b2">
<img width="565" alt="node(fix)" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/a9259f72-e28c-46f5-919d-a73583388df3">
<img width="531" alt="meta(fix)" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/5e169c25-e797-419b-8f06-0dfd1b1ce32b">

## 🩺 Survival Analysis of patients using clinical dataset:
<img width="708" alt="surv1" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/4353604e-8f7c-4a96-b40c-975723bf0485">
* While looking for impact of age on Overall Survival of patients, with the p-value of 0.62, the Kaplan-Meier curve shows no significant change across various age-groups.
<img width="1151" alt="surv2" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/78be4e93-bdec-482d-8d75-c6bcd7c4a3ae">
* While looking for impact of tumor stage on overall survival of patients, with the p-value of 0.00015, it shows the importance of the difference of the tumor stages and the survival chances. It shows the tumor stage to be a powerful prognosis at the time of the diagnosis for the overall survival of breast cancer patients.

## 📈 Result of Classification Model (Supervised and Unsupervised Learning Methods):
### Invasive and Non-Invasive cancer type classification:- 
* The Random Forest model performed extremely well in terms of precision and recall. The ANN model performed well in recognising invasive cancer, but it showed slightly less performance with non- invasive predictions. The KNN model displayed low performance for both the types of cancer as it struggled to predict the imbalanced class type.


<img width="605" alt="Screenshot 2024-03-18 at 16 58 06" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/09f15604-6adc-40b9-aadf-e6954a8ad080">
<img width="697" alt="Screenshot 2024-03-18 at 16 58 35" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/ee00011f-33ad-43d4-805c-21213e00e36d">

### Cluster analysis using PCA and K-Means:
* PCA is used to resolve high dimensional dataset into two principal components.
* The Elbow method is used to select the optimal number of clusters by choosing the k-value.
* The silhouette score reached its peak at 0.3929, showing similarity among the resultant clusters.

<img width="615" alt="Screenshot 2024-03-18 at 17 07 57" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/b82facad-beb0-40c3-a7d8-962e5dd12a92">

### Feature selection using Random Forest Method:-
* For each gene, an AUC score was determined, and the top 51 genes, each with an AUC score exceeding 0.7, were filtered for an in-depth examination, indicating their better predictive capacity. In our findings, particularly, the genes ‘SHE; and ‘SLC19A3’ stand out with their higher importance scores.
![top (1)](https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/c2f52be9-d3ae-4054-bd0f-6abf1d7f43cf)

### Improved classification models using engineered clusters of genes: 
* We performed cross-validation, the 5 -fold method is used to predict the performance of the machine learning models. The improved scores are shown below with the new improved ROC-curve.
  <img width="612" alt="Screenshot 2024-03-18 at 17 20 16" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/44dc3a62-cbd7-497b-a39e-191eb9e9a891">
<img width="585" alt="Screenshot 2024-03-18 at 17 19 39" src="https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/assets/91785440/b8f4dc1a-48fd-482b-b509-36850004fa68">

## Please Go Through [TCGA_BRCA- Machine learning models- HLD.docx](https://github.com/Gauravonthemixx/Text-Analytics-of-TED-Talks/files/14680993/TCGA_BRCA-.Machine.learning.models-.HLD.docx) for more info. 🔗



