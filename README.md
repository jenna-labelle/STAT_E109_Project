# STAT_E109_Project
## Final project for STAT E109: Predicting melanoma tumor grade using estimated cell composition

### Research question and motivation:

Cancer researchers frequently profile tumor tissue via “bulk” sequencing, which provides a cheap, quick, high-level look at the gene expression patterns for many samples. There are many publicly available cancer databases, such as the Cancer Genome Atlas [1], that provide this bulk gene expression data for large numbers of patient samples, as well as extensive clinical data such as age, patient outcome, and tumor grade. The major drawback to this bulk sequencing data is that all information on individual cell types within a tumor are lost. This information is critical, as recent work has highlighted the important anti-cancer role of immune cells within a tumor [2].

A more detailed approach is to sequence individual cells within a tumor- this provides extensive detail on the cell types present within a tumor but is much more expensive and time consuming. Integrating this “bulk” and “single cell” data, then, is of great interest to researchers. The goal of this project is to use the cell type composition information present in single cell samples to predict the composition in bulk samples, then use this information to predict tumor grade in bulk samples.

The plan for this project is two-fold: 
1) Using “single cell” melanoma data, generate a multivariate linear regression model to predict the proportion of different cell types in “bulk” melanoma samples based on expression of a few genes. 
2) Using these predicted cell type proportions for the “bulk” data, generate decision tree ensemble models to determine if these predicted cell type proportions can predict tumor grade.


### Data acquisition and processing:

Two main sources of data were used in this project: single cell sequencing data from Tirosh et al [3] and bulk sequencing data from The Cancer Genome Atlas [1]. Both data sources consist of large tables of genes x samples. From Tirosh et al, marker genes for 
6 cell types of interest were also obtained: T cells, B cells, macrophages, endothelial cells, CAFs, and melanoma cells. These marker genes, rather than the full set of ~20,000 genes, are used to predict cell composition in downstream analyses. 

In order to generate the linear regression models to predict these cell type proportions based on gene expression, a method referred to as “pseudo-bulking” was used. This essentially takes the single cell data and combines it into 50 random “pseudo samples” that are then comparable to actual bulk sequencing samples. 

Following this pseudo-bulking process, gene expression values for both the pseudo-bulked samples as well as the actual bulk samples were normalized via 3 standardly used methods: 1) Counts per million (each sample is normalized to have exactly 1 million gene counts) 2) log transformation 3) centering counts for each sample. For code to perform preprocessing, see Preprocessing.Rmd within the github repository [4].


### Predicting cell composition in pseudo-bulked samples: linear regression

The first goal of this project is to predict the cell composition in bulk sequencing samples. To accomplish this, the pseudo-bulked samples (where actual cell composition is known) were used to generate models for each cell type. Because both the outcome (cell proportion) and predictors (gene expression) are numeric variables, linear regression seemed to be a reasonable choice. However, this model assumes normality, which is not always the case for gene expression data. 

Six separate linear regression models are generated, one for each cell type. The input for each cell type consists of the known proportion of that cell type (outcome variables) and the expression values for 15 marker genes (predictor variables). 

Prior to generating the multivariate linear regression models, the data for each cell type was split into training (70%) and testing (30%). A variety of other testing/training splits were used (50/50, 60/40, and 80/20), but this training/testing split gave the most significant results, on average. To eliminate extraneous factors/genes, the AIC regression algorithm (direction = both) was used to keep only factors that are most significant. This method of regression is performed in a stepwise manner, where non-significant variables are iteratively added and removed until the final model is generated. For each cell type, 15 genes were input into the AIC regression algorithm and were subsequently reduced to:

  *	Bcells: 5 genes 
  * CAFs: 6 genes 
  * Endothelial cells: 11 genes 
  * Macrophages: 9 genes
  * Melanoma: 6 genes 
  * T cells: 6 genes

Multicollinearity may be an issue in gene expression data. To remove genes with high correlation, the VIF function was used for each cell types’ regression model. Overall, VIF values are relatively low, with all values below 20. A cutoff of 10 is generally a reasonable threshold; 6 genes across the 3 of the cell types showed a VIF > 10. These genes were removed, and additional iterative AIC regression models were generated for these 3 cell types. 

Most predictive factors are significant, with pvalues less than 0.05. While a few pvalues are higher than the 0.05 threshold, the AIC regression algorithm still weights these genes as predictive, and thus they were retained in the final model. Notably, the final B cell model contains only 2 predictive genes. However, based on biological knowledge of this cell type, this seems reasonable, as these markers are highly specific to B cells and may be sufficient to predict proportion of this cell type present.

Overall, all models are able to significantly predict cell type composition (p < 0.05). 

The models were then further validated using the testing datasets. For all cell types except for T cells, the model significantly predicted the proportion of each cell type in the pseudo-bulked samples (p<0.05). The T cell model is significant in the training dataset, suggesting a potential over-fitting issue. Because T.cells are the immune cell types that are able to directly kill tumor cells, they often correlate with a better prognosis. The fact that this model is unable to predict the proportion of T cells (and thus are not used downstream to predict T cell composition in the bulk sequencing data) is a major drawback of this analysis. Without T cell composition information, it may be more challenging to predict tumor grade. However, the models for the remaining 5 cell types are highly significant and can be used downstream to predict composition of these 5 cell types in bulk data.

For the 5 significant multivariate linear regression models (B.cells, CAFs, Endothelial cells, Macrophages, and melanoma), the proportion of these 5 cell types was predicted in bulk melanoma sequencing data (n=378).For CAFs and EndothelialCells, the majority of these predicted proportions are quite low, with most values less than 10%. This is what we may expect, however, as proportions of these cells are quite low in the single cell data as well (1.3% and 1.4%, respectively). 

### Predicting tumor grade based on predicted cell type composition: decision tree modeling

Bulk sequencing data was filtered to include only samples with tumor grade information available (n=116). Tumor grade was provided with high granularity, with each major class (Type1-IV) being split into additional subclasses. Because these sub-classes often contained only a few samples, tumor subtypes were simplified into their major classes. The data was then split into training (70%) and testing (30%) data sets (training settings: number=5, repeats=1, method=repeatedcv).

In order to predict tumor grade of bulk sequencing samples based on the predicted composition of cell types, a decision-tree model is a reasonable choice. Three ensemble methods were evaluated: bagging, random forest (final mtry value = 4), and boosting (xgbtree; settings: nrounds=500, max_depth=4, eta=0.28, gamma=2). For all models, the outcome variable is tumor grade and the predictors are the cell type composition predicted in the multivariate linear regression model. In addition to this, gender and age may be large confounders in the analysis and were included as additional covariates.  

For the boosting and random forest models, the most important cell type predictor is macrophages or Bcells. As hypothesized, then, immune cell proportion within a tumor tends to be the most important predictor of tumor stage.


The accuracyNull value (0.39) describes the baseline accuracy- that is, the accuracy if all observations were assigned to the maximum occurring group (in this case, Grade III). For the bagging model, accuracy is not improved further past this baseline accuracy (accuracy= 0.39; pvalue= 0.56). For both the random forest and boosting methods, accuracy is slightly improved beyond the baseline accuracy (accuracy: 0.54 and 0.44, respectively). However, the boosting model does not significantly predict tumor grade (p=0.31) while the random forest model is significant at the 95% confidence level (p=0.04). Thus, the random forest model would be selected for any further downstream analyses. Importantly, a random forest model built from age and gender alone cannot significantly predict tumor grade (p=0.432). 

However, while the random forest model is significant, the accuracy is still very low at only 53.7%. Plotting these predictions in the testing datasets highlights the low accuracy. The model is able to predict Stage III tumors with high accuracy, but Stages I, II, and IV have relatively low predictive accuracy. The code to generate all ensemble method models can be found in predictTumorGrade.Rmd within the github repository [4].


### Conclusions:

This project was composed of two main steps: predicting tumor cell type in pseudo-bulked data, and then using these predictions to predict tumor grade in bulk sequencing data. The first step was accomplished, as we were able to significantly predict tumor cell type from pseudo-bulked data for 5/6 cell types. Thus, for these 5 cell types, we reject the null hypothesis that gene expression does not predict cell composition. 

The second goal is to use this predicted cell type composition to predict tumor grade within bulk sequencing data. Here, we were able to generate a random forest model that significantly predicted tumor grade, and thus we reject the null hypothesis that cell type composition does not predict tumor grade. A caveat here is that, while significant, the accuracy of this model is quite low, and the model would obviously need to be improved before using it in any clinical or research setting.


### References:

[1] https://www.cancer.gov/tcga

[2] Gonzalez, H., Hagerling, C., & Werb, Z. (2018). Roles of the immune system in cancer: from tumor initiation to metastatic progression. Genes & development, 32(19-20), 1267–1284. https://doi.org/10.1101/gad.314617.118

[3] Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., 2nd, Treacy, D., Trombetta, J. J., Rotem, A., Rodman, C., Lian, C., Murphy, G., Fallahi-Sichani, M., Dutton-Regester, K., Lin, J. R., Cohen, O., Shah, P., Lu, D., Genshaft, A. S., Hughes, T. K., Ziegler, C. G., Kazer, S. W., … Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science (New York, N.Y.), 352(6282), 189–196. https://doi.org/10.1126/science.aad0501

[4] LaBelle (2022) STAT_E109_Project https://github.com/jenna-labelle/STAT_E109_Project


