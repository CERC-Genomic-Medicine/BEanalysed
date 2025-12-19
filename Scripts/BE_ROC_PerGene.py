#!/usr/bin/env python3

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

### To define grid ####
def suite(n):
    summ=0
    for i in range(0,n):
        summ=summ+i
    return summ

def is_prime(n):
    """Check if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def find_factorial_grid_approx_prime(product):
    # If the product is a prime number greater than 3, approximate to the next value
    if is_prime(product) and product > 3:
        product += 1  # Approximate to the next number        
    
    a, b = 1, product
    
    # Loop to find the factors of product that are closest to each other
    for i in range(1, int(product ** 0.5) + 1):
        if product % i == 0:
            a, b = i, product // i
            
    return a, b


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


def plot_roc_auc_per_gene_from_excel_summary(dfs, plotted_value, genes, positive_consequences=[], negative_consequences=[], seed=42, out='output.png'):
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
    gene_data_dict = {}
    for gene in genes: 
        datasets={}
        for label, data_full in dfs.items():
            ###### RAUC
            df=data_full.loc[data_full['proteins']==gene, :].copy()
            df.loc[:,'Truth'] = np.nan
            print(f"{gene} : ")
            df.loc[[row['Consequence_Detail'] in negative_consequences for ind, row in df.iterrows()],'Truth']=0
            print(f' {gene} :  Negative Controls')
            print(df.loc[[row['Consequence_Detail'] in negative_consequences for ind, row in df.iterrows()],'Consequence_Detail'].value_counts())
            df.loc[[row['Consequence_Detail'] in positive_consequences for ind, row in df.iterrows()],'Truth']=1
            print(f' {gene} :  Positive Controls')      
            print(df.loc[[row['Consequence_Detail'] in positive_consequences for ind, row in df.iterrows()],'Consequence_Detail'].value_counts())
            df = df.loc[df['Truth'].notna(), :]
            y_true = df['Truth']
            ranks=df["|".join([plotted_value, 'rank'])]
            inverted = max(ranks) + 1 - ranks ### Ranks are highest (1) to lowest (max) Rauc want scores with highest (max) lowest (min)
            y_pred = inverted
            datasets[label]=pd.DataFrame({'y_true': y_true, 'y_pred': y_pred})
            print(f'{label} : Total = {len(y_true)}, Positive = {sum(y_true==1)}, Negative = {sum(y_true==0)}')
            labels.append(label)  # filename as label
        gene_data_dict[gene] = datasets
    # Plot
    plt.figure(figsize=(10, 8))
    np.random.seed(seed)
    palette = get_color_palette(max([len(datasets) for datasets in gene_data_dict.values()]))
    grid=find_factorial_grid_approx_prime(len(genes))
    fig, axs = plt.subplots(grid[0],grid[1], figsize=(grid[1]*8,grid[0]*8), dpi = 300)
    axs = np.atleast_1d(axs).flatten()
    for ind in range(0,len(genes)):
        keys=sorted(gene_data_dict[genes[ind]].keys())
        for label in keys:
            data=gene_data_dict[genes[ind]][label]
            print(f'{genes[ind]}, {label} :')
            print(data)
            y_pred, y_true = np.asarray(data['y_pred']), np.asarray(data['y_true'])
            # Debug info
            print(f'{genes[ind]}, {label}:')
            print(f'  y_true dtype: {y_true.dtype}, shape: {y_true.shape}')
            print(f'  y_true unique values: {np.unique(y_true)}')
            print(f'  Positives: {(y_true == 1).sum()}, Negatives: {(y_true == 0).sum()}')
            # Skip if invalid
            if len(y_true) == 0 or len(np.unique(y_true)) < 2:
                print(f'  Skipping - insufficient data')
                continue
            fpr, tpr, _ = roc_curve(y_true, y_pred)
            roc_auc = auc(fpr, tpr)
            axs[ind].plot(fpr, tpr, lw=2, label=f'{label} (AUC = {roc_auc:.3f})', color=palette[keys.index(label) % len(palette)])
        axs[ind].plot([0, 1], [0, 1], linestyle='--', color='grey')
        axs[ind].set_xlabel('False Positive Rate')
        axs[ind].set_ylabel('True Positive Rate')
        axs[ind].set_title(f'{genes[ind]}')
        axs[ind].legend(loc='lower right')
        axs[ind].grid(False)
    fig.suptitle('Per Gene RAUC analysis', fontsize=32)
    fig.gca().set_box_aspect(1) 
    fig.savefig(f'{out}.pdf',format="pdf")


def plot_roc_auc_per_gene_per_consequence_from_excel_summary(dfs, plotted_value, genes, positive_consequences=[], negative_consequences=[], seed=42, out='output.png'):
    """
    Plot ROC AUC analysis with one panel per gene × label combination.
    Each panel shows separate ROC curves for each positive consequence vs negative controls.
    
    Parameters:
        dfs: dict of DataFrames (sheet_name -> DataFrame)
        plotted_value: 'pos' or 'neg' for rank column
        genes: list of genes to analyze
        positive_consequences: list of positive consequence categories
        negative_consequences: list of negative consequence categories
        seed: random seed
        out: output filename
    """
    
    # Build data structure: gene -> label -> consequence -> DataFrame
    gene_data_dict = {}
    panel_keys = []  # List of (gene, label) tuples for panels
    
    for gene in genes:
        gene_data_dict[gene] = {}
        for label, data_full in dfs.items():
            df = data_full.loc[data_full['proteins'] == gene, :].copy()
            
            # Get ranks for all rows
            rank_col = "|".join([plotted_value, 'rank'])
            if rank_col not in df.columns or df.empty:
                continue
                
            # Identify negative controls (shared across all curves in this panel)
            neg_mask = df['Consequence_Detail'].isin(negative_consequences)
            neg_df = df.loc[neg_mask].copy()
            
            if neg_df.empty:
                print(f'Skipping {gene}/{label} - no negative controls')
                continue
            
            consequence_data = {}
            for pos_consq in positive_consequences:
                # Get positives for this specific consequence
                pos_mask = df['Consequence_Detail'] == pos_consq
                pos_df = df.loc[pos_mask].copy()
                
                if pos_df.empty:
                    continue
                
                # Combine positive and negative
                combined = pd.concat([pos_df, neg_df], ignore_index=True)
                combined['Truth'] = combined['Consequence_Detail'].apply(
                    lambda x: 1 if x == pos_consq else 0
                )
                
                ranks = combined[rank_col]
                inverted = ranks.max() + 1 - ranks
                
                consequence_data[pos_consq] = pd.DataFrame({
                    'y_true': combined['Truth'].values,
                    'y_pred': inverted.values
                })
                
                print(f'{gene}, {label}, {pos_consq}: Pos={pos_mask.sum()}, Neg={neg_mask.sum()}')
            
            if consequence_data:
                gene_data_dict[gene][label] = consequence_data
                panel_keys.append((gene, label))
    
    if not panel_keys:
        print("No valid data to plot")
        return
    
    # Create grid for gene × label panels
    n_panels = len(panel_keys)
    grid = find_factorial_grid_approx_prime(n_panels)
    fig, axs = plt.subplots(grid[0], grid[1], figsize=(grid[1]*6, grid[0]*6), dpi=300)
    axs = np.atleast_1d(axs).flatten()
    
    palette = get_color_palette(len(positive_consequences))
    
    for idx, (gene, label) in enumerate(panel_keys):
        ax = axs[idx]
        consequence_data = gene_data_dict[gene][label]
        
        for consq_idx, pos_consq in enumerate(positive_consequences):
            if pos_consq not in consequence_data:
                continue
                
            data = consequence_data[pos_consq]
            y_true = np.asarray(data['y_true']).astype(int)
            y_pred = np.asarray(data['y_pred']).astype(float)
            
            # Skip if insufficient data
            if len(y_true) == 0 or len(np.unique(y_true)) < 2:
                print(f'Skipping {gene}/{label}/{pos_consq} - insufficient data')
                continue
            
            fpr, tpr, _ = roc_curve(y_true, y_pred)
            roc_auc = auc(fpr, tpr)
            
            ax.plot(fpr, tpr, lw=2, 
                    label=f'{pos_consq} (AUC = {roc_auc:.3f})',
                    color=palette[consq_idx % len(palette)])
        
        ax.plot([0, 1], [0, 1], linestyle='--', color='grey')
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(f'{gene} - {label}')
        ax.legend(loc='lower right', fontsize=8)
        ax.set_aspect('equal')
    
    # Hide unused subplots
    for idx in range(len(panel_keys), len(axs)):
        axs[idx].set_visible(False)
    
    fig.suptitle('ROC AUC per Gene × Dataset (by Positive Consequence)', fontsize=16)
    fig.tight_layout()
    fig.savefig(f'{out}_per_consequences.pdf', format="pdf")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Plot repeat correlations with experimental labels.")
    parser.add_argument('-I',"--input", required=True, dest='excel_file', help="MageCK files per genes")
    parser.add_argument('-R',"--Remove_sheet", required=False, dest='to_remove',  default =[], nargs='+', help="Sheets to remove")
    parser.add_argument('-G',"--Genes", required=False, dest='Genes', nargs='+', help="List of genes to analyse")
    parser.add_argument('--Negative_Control_Consequences', metavar='str', dest='Control_N_Consequences', required=False, nargs='+', default=['No predicted Mutation'], choices=['splice','non-sense','missense','synonymous','non-coding','regulatory', 'No predicted Mutation','N/A'], type=str, help='list of negative control consequences')
    parser.add_argument('--Positive_Control_Consequences', metavar='str', dest='Control_P_Consequences', required=False, nargs='+', default=['splice','non-sense'], choices=['splice','non-sense','missense','synonymous','non-coding','regulatory', 'No predicted Mutation','N/A'], type=str, help='list of positive control consequences')
    parser.add_argument('-V',"--value", required=False, dest='plotted_value', choices = {'pos', 'neg'}, default='neg', help="value to be plotted (ex Rank|Pos)")
    parser.add_argument('-O',"--out", required=False, dest='out', default='RAUC_per_Gene', help="Output images path (default : RAUC_per_Gene) **Do not include extentions**")
    args = parser.parse_args()
    all_sheets = pd.read_excel(args.excel_file, sheet_name=None)
    for k in args.to_remove:
        all_sheets.pop(k, None)
    dffs = pd.concat(all_sheets.values())
    valid_genes = set(dffs['proteins'].unique().tolist())
    if args.Genes:
        list_genes=args.Genes
        # Check for duplicates
        if len(args.Genes) != len(set(args.Genes)):
            duplicates = [g for g in args.Genes if args.Genes.count(g) > 1]
            parser.error(f"Duplicate genes provided: {set(duplicates)}")
        # Check for invalid genes
        invalid = set(args.Genes) - valid_genes
        if invalid:
            parser.error(f"Invalid genes: {invalid}. Valid options: {valid_genes}")
    else : 
        list_genes = valid_genes
    plot_roc_auc_per_gene_from_excel_summary(all_sheets, args.plotted_value, list_genes, positive_consequences=args.Control_P_Consequences, negative_consequences=args.Control_N_Consequences, seed=42, out=args.out)
    plot_roc_auc_per_gene_per_consequence_from_excel_summary(all_sheets, args.plotted_value, list_genes, positive_consequences=args.Control_P_Consequences, negative_consequences=args.Control_N_Consequences, seed=42, out=args.out)
if __name__ == "__main__":
    main()


