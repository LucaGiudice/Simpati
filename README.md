# Simpati: pathway-based classifier <img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-Simpati/main/images/Simpati_logo.png" width=350 align="right" style="border:4px solid black;" />


[![R-CMD-check](https://github.com/jokergoo/cola/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/cola/actions)
[ ![bioc](https://bioconductor.org/shields/downloads/devel/cola.svg) ](http://bioconductor.org/packages/stats/bioc/cola)
[ ![bioc](http://bioconductor.org//shields/lastcommit/devel/bioc/cola.svg) ](http://bioconductor.org/checkResults/devel/bioc-LATEST/cola/)


## Features

1. It finds novel enriched pathways between two classes in comparison
2. It exploits the novel paradigm of patient similarity networks to detect pathways as no other method
3. It generates information and data to further manual analysis and explorations
4. It can be applied of multiple omics
5. It does not require technical parameters, only the patient's/samples's information
6. It allows to plot an enriched pathway in form of patient similarity network to have a beautiful informative image for the paper which allows to highlight the similarities and dissimilarities between patients
7. It performs an outlier dection and retrieves how much each patient is likely to not fit in its own class defined by its clinical information
8. It allows to explore the signature pathway of a patient class with a graphical interface

## Citation

Giudice Luca, et al., Simpati: pathway-based classifier exploits patient similarity networks for cancer stage prediction

## Install

*Simpati* is available on R, you can install it by:

```r
devtools::install_github(
  repo="LucaGiudice/Simpati",
  ref = "main",
  dependencies = "Depends",
  upgrade = "always",
  quiet = F
)
```

## Vignettes

There are the following vignettes:

1. [A Quick Start of Using Simpati Package and introduction to the results](https://github.com/LucaGiudice/Simpati/blob/main/vignettes/Classification_Mutations_introduction.Rmd)
2. [Advanced workflow of Simpati performed with all the operations explicit](https://github.com/LucaGiudice/Simpati/blob/main/vignettes/Classification_Mutations_advanced.Rmd)

## Classification and pathway detection

<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-Simpati/main/images/workflow.png" />

The cores of Simpati workflow:

Frame A shows the workflow
  1. Patient profiles colored by the class of belonging are described by biological features as genes. A gene-gene interaction network together with pathways are further input data required by the software. 
  2. All profiles are individually propagated over the network. The profile’s values are replaced by scores that reflect the gene’s starting information and connections. 
  3. Simpati proceeds by creating a patient similarity network for each pathway (psPSN). The pairwise similarity evaluates how much two patients have a similar pathway activity. 
  4. The psPSN is decomposed into three components. Two with the intra-similarities of the class specific representative patients, while one with the inter-similarities between them.  If the similarities of one class dominate over the other two components, then the psPSN is signature. 
  5. The latter is ultimately used to classify. An unknown patient is classified based on how much is like the other patients and on how much fits in the class specific representatives. 

Frame B shows the effect of the propagation 
1. The position of the gene and its level of expression leads a patient to act on the other genes and so on a pathway in a unique way. 

Frame C shows the patient similarity used in Simpati
1. It evaluates how much the genes between two patients are close and high in term of propagation values. Two patients that act on a pathway from the same gene positions and with the same expression values get the maximum similarity. 

Frame D illustrates the biological interpretation behind a signature psPSN. 
1. One class is cohesive because the disease is leading the patients to alter similarly the same pathway. One class is sparse because led by multiple factors.

### Usage

Few lines of code to perfrom *Simpati* analysis:

```r
library(Simpati)

#Get omic-specific patient profiles and their clinical data
geno=tcga_data$LIHC$`LIHC_Mutation-20160128`$assay_df;see(geno)
info=tcga_data$LIHC$`LIHC_Mutation-20160128`$clin_df;see(info)
info=info[,c("patientID","pathologic_stage")]

#Set name of the dataset
dataset_name="LIHC"
#Set the semantic type of the disease for the disgnet enrichment
disease_type=tcga_data$LIHC$semantic_type
#Set key words associated to the patient's disease
key_words=tcga_data$LIHC$key_words

#Gene interaction network
net=huri_net_l$net_adj

Simp_res=wrapper_human_mutations(geno,info,net,pathways_l,dataset_name,disease_type,key_words, n_cores=5,test_run=T,seed=0)
```

### Plots

Following plot shows you an example of pathway-specific patient similarity network between LIHC Late (L) stage cancer patients and Early (E)

<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-Simpati/main/images/BIOCARTA_MAPK_PATHWAY%20source-MSIGDB_C2%20source-BIOCARTA_MAPK_PATHWAY%20down-inv.png" /> <img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-Simpati/main/images/HALLMARK_HEDGEHOG_SIGNALING%20source-MSIGDB_C2%20source-HALLMARK_HEDGEHOG_SIGNALING%20up-inv.png" />

## License

MIT @ Giudice Luca
