#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.utils import resample
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
from scipy.stats import rankdata
from scipy.stats import binomtest
import sys
import os
import re
import warnings
from itertools import pairwise  # Python 3.10+
from itertools import combinations, product


'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

Goal : To produce a RAUC/AUC plot and rank based analysis
'''


# Function to get color palettes
def get_color_palette(num_entries):
    """
    Determine the appropriate palette for the number of entries.

    Parameters:
    num_entries (int) : number of entries
    
    Returns:
    str: color palette
    """
    pastel1 = plt.get_cmap('Pastel1').colors
    tab20 = plt.get_cmap('tab20').colors
    if num_entries <= len(pastel1):
        return pastel1
    elif num_entries <= len(tab20):
        return tab20
    else:
        warnings.warn("Warning: Number of entries exceeds the available colors in tab20 palette.")
        return tab20




def plot_roc_auc_from_excel_summary(dfs, plotted_value, seed=42, positive_consq_add=('non-sense', 'splice','Positive Control Gene'),negative_genes=[], negative_consq_add=('Negative Control Gene','synonymous','No predicted Mutation'),positive_genes=[], out='output.png'):
    """
    Parameters:
        mageck_files: list of MaGeCK result file paths
        vep_file: path to VEP annotated file
        position: 'AA' or 'Nuc' for consequence positioning
        n_bootstrap: bootstrap iterations
        seed: random seed
        positive_consq: tuple of consequence categories for positive class
        negative_consq: tuple of consequence categories for negative class

    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.metrics import roc_curve, auc
    from sklearn.utils import resample
    import pandas as pd
    """

    # Parse VEP consequences
    labels = []
    datasets=[]
    for label, data_full in dfs.items():
        figB, axB= plt.subplots(1, 1,figsize=(15,15))
        for cat in sorted(data_full['Consequence_Detail'].unique()):
            r = data_full[data_full['Consequence_Detail'] == cat]["|".join([plotted_value, 'rank'])].sort_values()
            y = pd.Series(range(1, len(r)+1)) / len(r)
            axB.plot(r, y, label=cat)
        axB.set_xlabel('Rank')
        axB.set_ylabel('Cumulative Fraction')
        figB.legend()
        axB.set_title('Ranks by Category')
        figB.savefig("_".join([label,out,'rank_cumulative.pdf']),format="pdf")
        plt.close(figB)
        ###### RAUC
        df=data_full
        df['Truth']=None
        df.loc[[row['Consequence'] in negative_consq_add or row['proteins'] in negative_genes or row['Controls']=='negative_control'  for ind, row in df.iterrows()],'Truth']=0
        N_cons = sum([row['Consequence'] in negative_consq_add for ind, row in df.iterrows()])
        print(df.loc[[row['Consequence'] in negative_consq_add or row['proteins'] in negative_genes or row['Controls']=='negative_control'  for ind, row in df.iterrows()],'Consequence'].value_counts())
        df.loc[[row['Consequence'] in positive_consq_add or row['proteins'] in positive_genes or row['Controls']=='positive_control'  for ind, row in df.iterrows()],'Truth']=1
        print(df.loc[[row['Consequence'] in positive_consq_add or row['proteins'] in positive_genes or row['Controls']=='positive_control'  for ind, row in df.iterrows()],'Consequence'].value_counts())
        df = df.loc[df['Truth'].notna(), :]
        y_true = df['Truth']
        ranks=df["|".join([plotted_value, 'rank'])]
        inverted = max(ranks) + 1 - ranks ### Ranks are highest (1) to lowest (max) Rauc want scores with highest (max) lowest (min)
        y_pred = inverted
        datasets.append((y_pred, y_true))
        print(f'{label} : Total = {len(y_true)}, Positive = {sum(y_true==1)}, Negative = {sum(y_true==0)}')
        labels.append(label)  # filename as label
    # Plot
    plt.figure(figsize=(10, 8))
    np.random.seed(seed)
    palette = get_color_palette(len(datasets))

    for idx, (y_pred, y_true) in enumerate(datasets):
        y_pred, y_true = np.asarray(y_pred), np.asarray(y_true)
        pos_idx = y_true == 1
        neg_idx = y_true == 0

        n_min = min(pos_idx.sum(), neg_idx.sum())
        pos_sample = resample(y_pred[pos_idx], replace=True, n_samples=n_min)
        neg_sample = resample(y_pred[neg_idx], replace=True, n_samples=n_min)
        y_resampled = np.concatenate([np.ones(n_min), np.zeros(n_min)])
        y_pred_resampled = np.concatenate([pos_sample, neg_sample])

        fpr, tpr, _ = roc_curve(y_resampled, y_pred_resampled)
        roc_auc = auc(fpr, tpr)

        plt.plot(fpr, tpr, lw=2, label=f'{labels[idx]} (AUC = {roc_auc:.3f})', color=palette[idx % len(palette)])

    plt.plot([0, 1], [0, 1], linestyle='--', color='grey')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves')
    plt.legend(loc='lower right')
    plt.grid(False)
    plt.tight_layout()
    plt.gca().set_box_aspect(1) 
    plt.savefig(f'{out}.pdf',format="pdf")

def main():
    parser = argparse.ArgumentParser(description="Plot repeat correlations with experimental labels.")
    parser.add_argument('-I',"--input", required=True, dest='excel_file', help="MageCK files per genes")
    parser.add_argument('-R',"--Remove_sheet", required=False, dest='to_remove',  default =[], nargs='+', help="Sheets to remove")
    parser.add_argument('-V',"--value", required=True, dest='plotted_value', choices = {'pos', 'neg'}, help="value to be plotted (ex Rank|Pos)")
    parser.add_argument("--Gene_positive", required=False, dest='gene_positive', default =[], nargs='+', help='Which gene should be considered for Positive controls')
    parser.add_argument('-P',"--Positive_consequence", required=False, dest='positive_consquence', default= ['non-sense', 'splice'],nargs='+', help='Which annotation should be positive controls (default : ["non-sense", "splice"]')
    parser.add_argument("--Gene_negative", required=False, dest='gene_negative', default =[], nargs='+', help='Which gene should be considered for negative controls')
    parser.add_argument('-N',"--Negative_consequence", required=False, dest='Negative_consquence', default= ['No predicted Mutation'],nargs='+', help='Which annotation should be Negative controls (default : ["No predicted Mutation"]')
    parser.add_argument('-O',"--out", required=False, dest='out', default='RAUC.png', help="Output image path (e.g., output.png)")
    args = parser.parse_args()
    all_sheets = pd.read_excel(args.excel_file, sheet_name=None)
    for k in args.to_remove:
        all_sheets.pop(k, None)
    plot_roc_auc_from_excel_summary(all_sheets, args.plotted_value, seed=42, positive_consq_add=args.positive_consquence, negative_consq_add=args.Negative_consquence,  positive_genes= args.gene_positive , negative_genes=args.gene_negative, out=args.out)

if __name__ == "__main__":
    main()

