#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error


# In[3]:


import os


# In[90]:


os.chdir("D:/ProjectWork/Neuroblastoma/PC1 paper/")


# In[200]:


df = pd.read_csv('./Neuro/PC1Means/GSE66586.csv')
df['PC1V'] = df['PC1V']/100


# In[201]:


df


# In[93]:


def hill(x, A, n):
    return A**n/(A**n+x**n)


# In[94]:


def exp(x, A, n):
    return A*np.exp(-n*x)


# In[95]:


def lin(x, a, b):
    return a*x+b


# In[202]:


def fit_fun(df,f):
    popt, pcov = curve_fit(f, df['Swaps'], df['PC1V'], p0=[1, 1], maxfev=5000)
    y_fit = f(df['Swaps'], *popt)
    MSE = mean_squared_error(df['PC1V'], y_fit)
    R2 = r2_score(df['PC1V'], y_fit)
    print(str(f.__name__)+str(MSE)+'    '+str(R2))
    plt.plot(df['Swaps'], df['PC1V'], 'o', label='data')
    plt.plot(df['Swaps'], y_fit, label='fit')
    plt.ylabel('Mean (PC1-Variance)')
    plt.xlabel('Swaps')
    plt.title('GSE66586', fontdict = {'size': 18})
    plt.text(18, 0.3, str(f.__name__)+'\nA: '+str(round(popt[0], 3))+'\nn: '+str(round(popt[1], 3))+'\nMSE: '+str(round(MSE, 6))+'\nR2: '+str(round(R2, 3)))
    plt.savefig('GSE66586'+'_'+'fit-'+str(f.__name__)+'.png')
    plt.show()


# In[203]:


fit_fun(df, hill)


# In[204]:


fit_fun(df, exp)


# In[205]:


fit_fun(df, lin)

