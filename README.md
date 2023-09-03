# Random Swap analysis 
To generate histograms for PC1 variance for random comnbinations of NOR/MES genelists and Housekeeping genes, use the PC1_Variance_Histogram.r script 
Add pre-processed gene-expression matrix files as tab-delimited .txt files in the datasets folder.

To create boxplots with random swaps of NOR/MES genelist with housekeeping genes (one at a time) use the PC1_Swap.r script

# Linear Fit to PC1 Means
Add mean values of PC1 variance for corresponding number of swaps to a GSEID.csv file (column1: Number of Swaps, Column2: Mean Variance) in the PC1_Means folder and use the Linear_fit.py script to visualize the fit and obtain R-squared value, mean squared error value, slope and intercept for the fit.


