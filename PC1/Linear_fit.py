#%%
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
import os

#%%
os.chdir("./PC1/")
df = pd.read_csv('./Neuro/PC1_Means/GSEID.csv')

#%%
def fit_fun(df,f):
    popt, pcov = curve_fit(f, df['Swaps'], df['PC1V'], p0=[1, 1], maxfev=10000)
    y_fit = f(df['Swaps'], *popt)
    MSE = mean_squared_error(df['PC1V'], y_fit)
    R2 = r2_score(df['PC1V'], y_fit)
    print(str(f.__name__)+str(MSE)+'    '+str(R2))
    plt.plot(df['Swaps'], df['PC1V'], 'o', label='data')
    plt.plot(df['Swaps'], y_fit, label='fit')
    plt.ylabel('Mean (PC1-Variance)')
    plt.xlabel('Swaps')
    plt.title('GSE17714', fontdict = {'size': 18})
    plt.text(0.7, 0.22, str(f.__name__)+'\nA: '+str(round(popt[0], 3))+'\nn: '+str(round(popt[1], 3))+'\nMSE: '+str(round(MSE, 6))+'\nR2: '+str(round(R2, 3)))
    plt.savefig('GSE17714'+'_'+'fit-'+str(f.__name__)+'.png')
    plt.show()

#%%
def lin(x, a, b):
    return a*x+b

#%%
fit_fun(df, lin)

