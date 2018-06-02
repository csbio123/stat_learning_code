
# coding: utf-8

# In[89]:


import pandas as pd 
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from scipy import stats
import sys
from sklearn import metrics
import os
import seaborn as sns



def analyse_results(input_dir, output_path):

    results = [y for x in os.walk(input_dir) for y in glob.glob(os.path.join(x[0], '*predictions*'))]
    data = load_results(results)
    os.makedirs(output_path, exist_ok=True)

    data_mean = data.groupby(["type_1"]).mean().sort_values("prediction_rf_ape").round(2)
    data_indx = data.set_index(["type_1"]).round(2)

    data_indx.to_csv("{}/data_combined.csv".format(output_path), header=True)
    data_mean.to_csv("{}/data_combined_mean.csv".format(output_path), header=True)
    
    return data_indx, data_mean

def load_results(results):
    data = []
    for f in results:
        d = pd.read_csv(f, sep=" ")
        file = os.path.basename(f)
        f = f.replace(input_dir, "")
        f = f.replace(file, "")
        dir_splits = f.split("/")
        type_1 = dir_splits[1]
        type_2 = dir_splits[2]
        d["type_1"] = type_1
        #d["type_2"] = type_2
        data.append(d)
    data = pd.concat(data)
    data = data.round(2)
    return data

def select_extremes(data, data_mean, number):
    
    best = data_mean.index.values[:number]
    worst = data_mean.index.values[-number:]
    show = np.concatenate([best, worst])

    return data.loc[show].reset_index()

def violin_plot(data_show, output_dir):
    figure= plt.figure(figsize=[10, 7])
    sns.set_style("whitegrid")
    ax = sns.violinplot(data=data_show, x="type_1", y="prediction_rf_ape",  
                     fliersize=0, palette=sns.color_palette("Set1", n_colors=8, desat=.5))

    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.xaxis.grid(color='white', linestyle='dashed')
    plt.xlabel('Feature classes')
    plt.ylabel('Absolute Percentage Error [%]')
    plt.tight_layout()
    plt.savefig("{}/violin_plot.pdf".format(output_dir))
    
    
def boxplot(data_show, output_dir):
    figure= plt.figure(figsize=[10, 7])
    sns.set_style("whitegrid")
    ax = sns.boxplot(data=data_show, x="type_1", y="prediction_rf_ape",  
                     fliersize=0, palette=sns.color_palette("Set1", n_colors=8, desat=.5))

    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.xaxis.grid(color='white', linestyle='dashed')
    plt.xlabel('Feature classes')
    plt.ylabel('Absolute Percentage Error [%]')
    plt.tight_layout()
    plt.savefig("{}/boxplot_plot.pdf".format(output_dir))
    
    
# In[34]:
if __name__ == "__main__":
        input_dir = sys.argv[1]
        output_dir = sys.argv[2]


        data, data_mean = analyse_results(input_dir, output_dir)
        data_show = select_extremes(data, data_mean, 3)

        print(data_show)
        violin_plot(data_show, output_dir)
        boxplot(data_show, output_dir)

