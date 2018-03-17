# In[27]:
import pandas as pd
import glob
import numpy as np
import os
from scipy.stats import chi2_contingency


# In[28]:


def run(path, target, output):
    input_dir = os.path.basename(path)

    if not os.path.exists(output):
        os.makedirs(output)

    print("Input directory {} is being processed".format(path))
    print("Target {} is being processed".format(target))

    allFiles = glob.glob(path + "/*.csv")
    frame = pd.DataFrame()
    list_ = []
    target = pd.read_csv(target)
    target.iloc[-1, 1] = len((target["protein"].drop_duplicates()))
    target.iloc[-1, 2] = len((target["gene"].drop_duplicates()))
    target = target.iloc[-1, :]

    for i, file_ in enumerate(allFiles):
        df = pd.read_csv(file_, index_col=None, header=0)
        df["subset"] = os.path.basename(file_)
        df.iloc[-1, 1] = len((df["protein"].drop_duplicates()))
        df.iloc[-1, 2] = len((df["gene"].drop_duplicates()))
        list_.append(df.iloc[-1, :])
    combined = pd.concat(list_, axis=1).transpose().set_index("subset").sort_index()
    combined = combined.apply(pd.to_numeric, errors='ignore')

    print("Doing some stats...")

    chi_p_value = combined.copy()
    chi_statistic = combined.copy()

    # In[29]:

    for i in range(0, combined.shape[0]):

        row = combined.iloc[i,]
        for j in range(1, len(combined.columns)):
            obs = np.array([[row.iloc[j], row.iloc[0]], [row.iloc[j], row.iloc[0]]])
            xi_s = (chi2_contingency(obs, correction=False))
            chi_p_value.iloc[i, j] = xi_s[1]
            chi_statistic.iloc[i, j] = xi_s[0]

    # In[32]:
    chi_p_value_f = output + "/{}_chi_p_value.csv".format(input_dir)
    chi_statistic_f = output + "/{}_chi_statistic.csv".format(input_dir)
    summary_f = output + "/{}_summary.csv".format(input_dir)

    print("Writing Chi-squared p-value file to {}".format(chi_p_value_f))
    print("Writing Chi-statistic file to {}".format(chi_statistic_f))
    print("Writing Summary file to {}".format(summary_f))

    chi_p_value.to_csv(chi_p_value_f)
    chi_statistic.to_csv(chi_statistic_f)
    combined.to_csv(summary_f)

if __name__ == "__main__":

    path =   '/Users/ti1/Downloads/config/gwas/output_test'
    target = '/Users/ti1/Downloads/config/gwas/output_test/target_test.fasta_result.csv'
    output = '/Users/ti1/Downloads/config/gwas/summary_outputs'

    run(path, target, output)
