# Gene-Expression Data Retrieval and Pre-processing

The list of datasets used in this analysis is provided in supplementary table[`NeuroblastomaSI.xlsx`](https://www.tandfonline.com/action/downloadSupplement?doi=10.1080%2F15384047.2024.2301802&file=kcbt_a_2301802_sm6771.xlsx) of the [paper](https://doi.org/10.1080/15384047.2024.2301802).

To retrieve and pre-process gene-expression data:

- For **microarray datasets**, follow instructions from [EMT_score_calculation](https://github.com/priyanka8993/EMT_score_calculation).
- For **RNA-sequencing datasets**, refer to [EMT_Scoring_RNASeq](https://github.com/sushimndl/EMT_Scoring_RNASeq).

---

# Random Gene-Swap Analysis

To explore PC1 variance using random gene combinations:

- Use `PC1_Variance_Histogram.r` to generate histograms for PC1 variance across random combinations of NOR/MES gene lists and housekeeping genes.
- Add pre-processed gene-expression matrix files as tab-delimited `.txt` files in the `datasets/` folder.

- Use `PC1_Swap.r` to generate boxplots by swapping NOR/MES genes with housekeeping genes (one at a time).

---

# Linear Fit to PC1 Means

To visualize the relationship between gene swaps and PC1 variance:

1. Create a `GSEID.csv` file in the `PC1_Means/` folder with the following structure:  
   - **Column 1**: Number of Swaps  
   - **Column 2**: Mean Variance

2. Run `Linear_fit.py` to generate the linear fit and obtain:
   - R-squared value  
   - Mean squared error  
   - Slope and intercept

---

# PCA, K-Means Clustering, and GSEA

To perform dimensionality reduction and enrichment analysis:

- Add pre-processed gene-expression matrix files as tab-delimited `.txt` files to the `data/` folder.
- Run `PCA.py` for Principal Component Analysis and K-Means clustering.
- Use `GSEA.py` to perform Gene Set Enrichment Analysis on the resulting clusters.
- Provide a gene signature `signature.gmt` file as input for GSEA. Gene signatures can be obtained from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb).

- Create an `Output/` folder to store:
  - PCA plots  
  - GSEA results

---

## Citation

If you use this repository or its contents in your work, please cite the following publication:

**Mutually exclusive teams-like patterns of gene regulation characterize phenotypic heterogeneity along the noradrenergic-mesenchymal axis in neuroblastoma**  
*Cancer Biology & Therapy (2024)*  
DOI: [10.1080/15384047.2024.2301802](https://doi.org/10.1080/15384047.2024.2301802)

---


