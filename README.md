# Gene-expression data retrieval and pre-processing
The list of datasets has been provided in 'NeuroblastomaSI.xlsx'. 

Follow instructions from https://github.com/priyanka8993/EMT_score_calculation to retrieve and pre-process microarray data and https://github.com/sushimndl/EMT_Scoring_RNASeq for RNA-sequencing datasets.

# Random gene-swap analysis 
To generate histograms for PC1 variance for random comnbinations of NOR/MES genelists and housekeeping genes, use the PC1_Variance_Histogram.r script. 
Add pre-processed gene-expression matrix files as tab-delimited .txt files in the datasets folder.

To generate boxplots with random swaps of NOR/MES genelist with housekeeping genes (one at a time) use the PC1_Swap.r script

# Linear fit to PC1 means
Add mean values of PC1 variance for corresponding number of swaps to a GSEID.csv file (column1: Number of Swaps, Column2: Mean Variance) in the PC1_Means folder and use the Linear_fit.py script to visualize the fit and obtain R-squared value, mean squared error value, slope and intercept for the fit.

# PCA, K-Means and GSEA
To perform PCA and K-means clustering on gene-expression data, add pre-processed gene-expression matrix files as tab-delimited .txt files in the data folder. Use PCA.py script for PCA and K-Means and GSEA.py for GSEA on the two generated clusters. Use the signature.gmt file as input for GSEA. Create an Output folder to store GSEA outputs and PCA plots.
