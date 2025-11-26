#!/usr/bin/env python3

import argparse
import re
import os
import math
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
from itertools import pairwise  # Python 3.10+
from itertools import combinations, product


plt.rcParams.update({'font.size': 20})


def parse_treatment_levels(xvar):
	return [sorted(group.split('/'), key=len, reverse=True) for group in xvar.split(',')]

def extract_labels_from_filename(filename, levels):
	fname = os.path.basename(filename)
	labels = []
	for group in product(*levels):
		hits = [lvl for lvl in group if re.search(rf'{re.escape(lvl)}', fname)]
		sized = len(group)
		if len(hits) != sized:
			continue
		else :
			labels.extend(hits)
	if len(labels) > sized :
		raise ValueError(f"Filename '{fname}' must contain exactly one of {product(*levels)}, found: {labels}")
	return labels


def extract_repeat_id(sgrna):
	match = re.search(r'_r(\d+)$', sgrna)
	return int(match.group(1)) if match else None

def group_by_base_name(df):
	indexed = df.index.astype(str)  # Ensure string index
	df = df.copy()
	df['base'] = indexed.str.replace(r'_r\d+$', '', regex=True)
	df['repeat'] = indexed.map(extract_repeat_id)
	cols = ['base', 'repeat'] + [c for c in df.columns if c not in ['base', 'repeat']]
	return df[cols]



def make_dataset(file_list, var, xvar) :
	dfs = []
	for f in file_list:
		label = "_".join(extract_labels_from_filename(f, parse_treatment_levels(xvar)))
		df = pd.read_csv(f, sep="\t")
		df.index=df['sgrna']
		df[label]=df[var]
		df['label'] = label
		df = df[label]
		dfs.append(df)
	data = pd.concat(dfs, axis=1)
	data = group_by_base_name(data)
	return data

def plot_repeat_comparisons_grid_2(df, value_columns, out):
	repeats = df['repeat'].unique()
	repeats.sort()
	repeat_pairs = list(itertools.combinations(repeats, 2))
	fig, axes = plt.subplots(len(repeat_pairs), len(value_columns), figsize=(6 * len(value_columns), 6 * len(repeat_pairs)))
	if len(value_columns) == 1 and len(repeat_pairs) == 1:
		axes = [[axes]]
	elif len(value_columns) == 1:
		axes = [[ax] for ax in axes]
	elif len(repeat_pairs) == 1:
		axes = [axes]

	for col_index, col in enumerate(value_columns):
		for row_index, (r1, r2) in enumerate(repeat_pairs):
			print(f" {col} : {r1} vs {r2} ")
			ax = axes[row_index][col_index]
			d1 = df[df['repeat'] == r1][['base', col]].set_index('base')
			d2 = df[df['repeat'] == r2][['base', col]].set_index('base')
			merged = d1.join(d2, lsuffix=f'_r{r1}', rsuffix=f'_r{r2}', how='inner')
			merged = merged.dropna()
			print(merged.columns)

			if merged.empty:
				ax.text(0.5, 0.5, f"No overlap for r{r1} vs r{r2}", ha='center', va='center')
				ax.axis('off')
				continue

			x = merged[f'{col}_r{r1}']
			y = merged[f'{col}_r{r2}']
			ax.scatter(x, y, alpha=0.6, s=10)
			ax.set_title(f"{col}" if row_index == 0 else "")
			ax.set_xlabel(f"repeat {r1}")
			ax.set_ylabel(f"repeat {r2}")
			lims = [
		    min(ax.get_xlim()[0], ax.get_ylim()[0]),
		    max(ax.get_xlim()[1], ax.get_ylim()[1]),
			]
			ax.plot(lims, lims, 'k--', linewidth=1)  # black dashed 1:1 line
			# Compute and annotate correlatio
			if len(x) > 1 and len(y) > 1:
				r, _ = pearsonr(x, y)
				ax.text(0.95, 0.05, f"r = {r:.2f}", transform=ax.transAxes,
					ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
	for ax in axes.flatten():
		ax.set_aspect('equal', adjustable = 'datalim')
	plt.tight_layout()
	plt.savefig(out, format='pdf')
	plt.close()


def main():
	parser = argparse.ArgumentParser(description="Plot repeat correlations with experimental labels.")
	parser.add_argument('-F',"--files", nargs='+', help="Input TSV files")
	parser.add_argument("--var", required=False, default='LFC', help="Column to analyze (e.g., LFC, score)")
	parser.add_argument('-o',"--out", required=False, default='repeat_on_repeat.pdf', help="Output image path (e.g., output.png)")
	parser.add_argument("--xvar", required=True, help="Comma-separated experimental conditions (e.g., UNT/TREAT,KO/WT)")
	parser.add_argument('--pivot', action='store_true', help="Pivot the grid") 
	args = parser.parse_args()
	data  = make_dataset(args.files, args.var, args.xvar)
	print(data.columns)
	plot_repeat_comparisons_grid_2(data, data.columns[2:], args.out)

if __name__ == "__main__":
	main()
