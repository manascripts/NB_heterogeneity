import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.pyplot import figure
import gseapy as gp

df=pd.read_csv('./GSEA_PCA/Data/GSEID.txt', sep= " ",header=None).transpose()
df.columns=df.iloc[0]
df=df.drop([0])
df=df.fillna(0)

cutoff = 10; #cutoff is where you set the classification between clusters
cls=[]
for i in range(len(PC['PC1'])):
    if PC['PC1'][i] < cutoff:
        print('NOR')
        cls.append('NOR')
    else:
        print('MES')
        cls.append('MES')

df1=df.transpose()
gs_res = gp.gsea(data = df1,
                 gene_sets="./GSEA_PCA/signature.gmt",
                 cls= cls,
                 # set permutation_type to phenotype if samples >=15
                 permutation_type = 'phenotype',
                 permutation_num = 100,
                 outdir='../Output',
                 no_plot = False, # Skip plotting
                 method='signal_to_noise',
                 processes = 4, seed = 7,
                 format = 'png')
