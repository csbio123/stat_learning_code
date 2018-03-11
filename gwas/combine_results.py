
# coding: utf-8

# In[27]:


import pandas as pd 
import glob
import numpy as np
import os
from scipy.stats import chisquare
import sys
from scipy.stats import chi2_contingency


# In[28]:


path = sys.argv[1]  #"/Users/ti1/Downloads/test/" # sys.argv[0] 
target = sys.argv[2] #"/Users/ti1/Downloads/query_files/gwas_88.fasta.target.fasta_result.csv"

allFiles = glob.glob(path + "/*.csv")
frame = pd.DataFrame()
list_ = []
target = pd.read_csv(target)
target.iloc[-1, 1] = len((target["protein"].drop_duplicates()))
target.iloc[-1, 2] = len((target["gene"].drop_duplicates()))
target = target.iloc[-1, :]

for i, file_ in enumerate(allFiles):
    df = pd.read_csv(file_,index_col=None, header=0)
    df["subset"] = os.path.basename(file_)
    df.iloc[-1, 1] = len((df["protein"].drop_duplicates()))
    df.iloc[-1, 2] = len((df["gene"].drop_duplicates()))
    list_.append(df.iloc[-1, :])
combined = pd.concat(list_, axis=1).transpose().set_index("subset").sort_index()
combined = combined.apply(pd.to_numeric, errors='ignore')
chi_p_value = combined.copy()
chi_statistic = combined.copy()


# In[29]:


for i in range(0, combined.shape[0]):
    
    row = combined.iloc[i, ]
    for j in range(1, len(combined.columns)):
        
        obs = np.array([[row.iloc[j], row.iloc[0]], [row.iloc[j], row.iloc[0]]])
        xi_s = (chi2_contingency(obs, correction=False))
        chi_p_value.iloc[i, j] = xi_s[1]
        chi_statistic.iloc[i, j] = xi_s[0]


# In[32]:


chi_p_value.to_csv(path + "/chi_p_value.csv")
chi_statistic.to_csv(path + "/chi_statistic.csv")
combined.to_csv(path + "/summary.csv")
