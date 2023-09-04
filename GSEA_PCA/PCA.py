#%%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.pyplot import figure

#%%
os.chdir("./Datasets/GSEID/")
loc=os.getcwd()+"\\"
match=[]
for file in os.listdir(loc):
    if file.endswith("_gene-exp2.txt"):
        match.append(file)

#%%
def data_std(geneset,component):
    df=pd.read_csv(geneset,sep=" ",header=None).transpose()
    df.columns=df.iloc[0]
    df=df.drop([0])
    df=df.fillna(0)
    X_std = StandardScaler().fit_transform(df) #scaling to unit variance
    pca = PCA(n_components=component)
    principalComponents = pca.fit_transform(X_std)
    PC=pd.DataFrame(principalComponents,columns=["PC"+str(x+1) for x in range(component)])
    return PC,pca

#%%
def variance_per(pca,geneset):
    features = range(pca.n_components_)
    figure(figsize=(15, 7.5), dpi=100)
    plt.rcParams.update({'font.size': 30})
    plt.bar(features, pca.explained_variance_ratio_, color='black')
    plt.xlabel('PCA features')
    plt.ylabel('variance %')
    plt.xticks(features)

#%%
def kmean_PC(PC,geneset,pca):
    from matplotlib.pyplot import figure
    custom_colors = ['green', 'purple']
    figure(figsize=(15, 10), dpi=100)
    plt.rcParams.update({'font.size': 32})
    kmeans = KMeans(n_clusters=2, init='k-means++', max_iter=300, n_init=10, random_state=0)
    pred_y = kmeans.fit_predict(PC)
    pred = kmeans.predict(PC)
    plt.scatter(PC["PC1"],PC["PC2"], c = [custom_colors[label] for label in pred], s = 200)
    plt.xlabel("PC1"+' '+'('+str(round(xx.explained_variance_ratio_[0]*100, 2))+'%'+')', fontname="Arial", fontsize=45, fontweight="bold")
    plt.ylabel("PC2"+' '+'('+str(round(xx.explained_variance_ratio_[1]*100, 2))+'%'+')', fontname="Arial", fontsize=45, fontweight="bold")
    plt.title(geneset.split("_")[0], fontweight="bold", fontname="Arial", fontsize=48)
    plt.savefig('../Output/'+geneset.split("_")[0]+'_kmeans_PC1_PC2.png')
    plt.show()

#%%
for geneset in match:
    PC,xx=data_std(geneset,component=4)
    variance_per(xx,geneset)
    kmean_PC(PC,geneset,xx)
