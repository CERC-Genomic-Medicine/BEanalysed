#!/usr/bin/env python3
import argparse
import sys
import re
import warnings
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from itertools import combinations, product
from collections import defaultdict
from scipy.stats import pearsonr, rankdata, binomtest

parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-I', '--input', metavar='FILE', dest='excel_file', required=True, type=str, help='')
parser.add_argument('-R',"--Remove_sheet", required=False, dest='to_remove',  default =[], nargs='+', help="Sheets to remove")
parser.add_argument('-X', '--comparison', metavar='string', dest='comparison', required=True, type=str, help='')
parser.add_argument('-C', '--color', dest='color_type', required=False,default=None , type=str, choices={None,'Consequence'}, help='What are the colors and shapde of the points')
parser.add_argument('-A', '--alpha', dest='alpha_type', required=False,default=None , type=str, choices={None,'all_little_bit','Significance'}, help='What are the colors of the points')
parser.add_argument('-B', '--biological_line', dest='B_line', required=False,default=None , type=float, help='plotting the biolofical siginifcance line')
parser.add_argument('-S', '--Significance_threshold', dest='p_thresh', required=False,default=None , type=float, help='P-value significance threshold')
parser.add_argument('-e', '--elements_to_plot', dest='element_plot', required=False,default='all' , type=str, choices={"all",'coding_only'}, help='Which elements are to be plotted (all : no filter [default], coding_only : no non-coding and no controls)')
parser.add_argument('-F', '--Square_Format', dest='SQformat', action='store_true', help='represent as a square 1:1 x/y axis in the same range')
parser.add_argument('-P', '--protein', dest='given_protein',  default =[], nargs='+', help='Plot_by_Protein')
parser.add_argument('-V', '--Variable', metavar='string', dest='var', required=False, default='lfc', type=str, help='')
parser.add_argument('-O', '--Output', dest='output', required=False,default=None , type=str, help='Output Name, by default --comparison argument .png if not --protein else --comparison argument_protein.png')







plt.rcParams.update({'font.size': 35})

coding_seq = ['synonymous','missense','splice','non-sense']
consequence_mapping = {
	'synonymous':   ('green', 'D'),     # diamond
	'missense':     ('purple', 'o'),    # circle
	'non-sense':    ('red', 's'),       # square
	'splice':       ('gold', '^'),      # triangle
	'regulatory':   ('blue', 'X'),      # X marker
	'non-coding':   ('lightgray', 'v'), # downward triangle
	'No predicted Mutation': ('gray', 'P'),  # plus-filled
	None:           ('silver', ',')     # default fallback
}



def add_legend_2(ax, consequence_mapping=None, size_map=None, transparency=None, loc='upper left',):
	"""
	Add a legend to the plot for the given consequence mapping and P-value sizes.
	"""
	handles = []
	subtitle_fontsize = 'small'

	if consequence_mapping:
		consequence_legend_elements = [
			mlines.Line2D([], [], color=value[0], marker=value[1], linestyle='None', markersize=18, label=key)
			for key, value in consequence_mapping.items()
		]
		handles.extend(consequence_legend_elements)

	if size_map:
		size_legend_elements = [
			mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=size, label=f'{pvalue}')
			for pvalue, size in size_map.items()
		]
		handles.extend(size_legend_elements)

	if transparency:
		transparency_legend_elements = [
			mlines.Line2D([], [], color='purple', marker='o', linestyle='None', markersize=18, alpha=alpha, label=label)
			for alpha, label in transparency
		]
		handles.extend(transparency_legend_elements)

	if handles:
		ax.legend(handles=handles, loc=loc, bbox_to_anchor=(1.05, 1),
				  title='Legend', handlelength=1, fontsize=subtitle_fontsize,
				  labelspacing=1.25, frameon=False)

def plot_scatter(ax, data, variable, label_a, label_b, consequence_mapping=None, xlines=None, ylines=None):
	if consequence_mapping : 
		for Consequence, group in data.groupby('Consequence'):
			print(Consequence)
			ax.scatter(
				group["_".join([variable,'df1'])],
				group["_".join([variable,'df2'])],
				label=Consequence,
				marker=consequence_mapping.get(Consequence)[1],
				alpha=group['alpha'],
				c=consequence_mapping.get(Consequence)[0],
				linewidth=3,
				s=200
			)
	else :
		ax.scatter(
			data["_".join([variable,'df1'])],
			data["_".join([variable,'df2'])],
			alpha=data['alpha'],
			c='black',
			linewidth=0.5,
			s=200
		)
	r, _ = pearsonr(data["_".join([variable,'df1'])], data["_".join([variable,'df2'])])
	print(_)
	ax.text(0.95, 0.05, f"r = {r:.2f}", transform=ax.transAxes,
	ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
	if xlines:
		for line in xlines :
			if line:
				ax.axvline(x=line, color='red', linestyle='--', linewidth=1)
	if ylines:
		for line in ylines :
			if line:
				ax.axhline(y=line, color='red', linestyle='--', linewidth=1)
	ax.set_xlabel(f"{label_a}")
	ax.set_ylabel(f"{label_b}")



def find_strict_matched_pairs(encoded, items):
	x, y = encoded.split('/')
	item_map = {}
	print(items)

	for item in items:
		parts = item.split('_')
		x_count, y_count = parts.count(x), parts.count(y)
		if x_count + y_count != 1:
			raise ValueError(f"Each item must contain exactly one of {x} or {y}: {item}")
		if x_count > 1 or y_count > 1:
			raise ValueError(f"Too many encoded comparison ({x}/{y}) in: {item}")
		idx = parts.index(x) if x in parts else parts.index(y)
		key = tuple(None if p == x or p == y else p for p in parts)
		item_map.setdefault((key, idx), {})[x if x in parts else y] = item

	matched = []
	for (key, idx), group in item_map.items():
		if x in group and y in group:
			pa, pb = group[x].split('_'), group[y].split('_')
			if len(pa) != len(pb):
				raise ValueError(f"Length mismatch in: {group[x]} / {group[y]}")
			# Ensure only the encoded token differs
			for i in range(len(pa)):
				if i == idx:
					if pa[i] != x or pb[i] != y:
						raise ValueError(f"Encoded mismatch at position {idx}: {group[x]} / {group[y]}")
				else:
					if pa[i] != pb[i]:
						continue  # collect diff later
			matched.append((group[x], group[y], idx))
		else:
			raise ValueError(f"Incomplete pair: {group}")
	return matched

def build_matrix_dict(encoded, items):
	print(encoded)
	matched = find_strict_matched_pairs(encoded, items)
	len_items = [len(i.split("_")) for i in items ]
	items_split = [i.split("_") for i in items ]
	variable_list = [list(set([f[indexe] for f in items_split])) for indexe in range(0,len_items[0]) ]
	encoded_position=set([i[2] for i in matched])
	print(encoded_position)
	if len(encoded_position)>1:
		raise ValueError(f"Each sheets must contain the comparison variable at the same position ('_' delimited)")
	else:
		encoded_position=list(encoded_position)[0]
	index_rows_tokens = [next(
	(i for i in range(len(variable_list))
	 if len(set(variable_list[i])) > 1 and i != encoded_position),
	None  # if there is one one graph to plot
	)]
	matrix_dict = dict()
	row_tokens = [ variable_list[i] for i in index_rows_tokens]
	row_tokens = ["_".join(list(p)) for p in product(*row_tokens)]
	col_tokens = [element for i, element in enumerate(variable_list) if not i in index_rows_tokens and i != encoded_position]
	col_tokens = ["_".join(list(p)) for p in product(*col_tokens)]
	col_index = [i for i, element in enumerate(variable_list) if not i in index_rows_tokens and i != encoded_position]
	print(f'row tokens : {row_tokens}, index rows tokens : {index_rows_tokens}, col index : {col_index}, col tokens : {col_tokens}')
	for match in matched :
		list_items = match[0].split('_')
		row_tok = "_".join([i for i in [list_items[j] for j in index_rows_tokens]])
		col_tok = "_".join([i for i in [list_items[j] for j in col_index]])
		row=row_tokens.index(row_tok)
		col=col_tokens.index(col_tok)
		matrix_dict[(row,col)]=(match[0],match[1])
	return matrix_dict


if __name__ == '__main__':
	args = parser.parse_args()
	try :
		all_sheets = pd.read_excel(args.excel_file, sheet_name=None)
	except: 
		raise ValueError('Main input data (-i) does not have a supported format')
	if args.alpha_type == 'Significance' and args.p_thresh is None:
		parser.error("Argument --Significance_threshold is required when --alpha is 'Significance'")
	for k in args.to_remove:
		all_sheets.pop(k, None)
	var=args.var
	blines=dict()
	if args.B_line:
		for key in all_sheets.keys():
			data_full = all_sheets[key]
			distribution = data_full.loc[data_full['Controls']=='negative_control',args.var]
			blines[key]=(np.quantile(pd.to_numeric(distribution),float(args.B_line)),np.quantile(pd.to_numeric(distribution),1-float(args.B_line) ))
	pairs_matrix_dic = build_matrix_dict(args.comparison, list(all_sheets.keys()))
	sizex=max([i[0] for i in pairs_matrix_dic.keys()]) + 1
	sizey=max([i[1] for i in pairs_matrix_dic.keys()]) + 1
	print(f'size = ({sizex}, {sizey})')
	print(pairs_matrix_dic)
	#fig, axes = plt.subplots(sizex,sizey,figsize= (25,50),subplot_kw={'aspect': 'equal'})
	fig, axes = plt.subplots(sizex,sizey,figsize=(20 * sizey+6*((sizex*sizey)**0.5), 18 * sizex) ,subplot_kw={'aspect': 'equal'})
	#fig, axes = plt.subplots(sizex,sizey,figsize=(43,36) ,subplot_kw={'aspect': 'equal'})
	comparison_label = args.comparison.replace('/','_')
	for position in pairs_matrix_dic.keys():
		label_a = pairs_matrix_dic[position][0]
		label_b = pairs_matrix_dic[position][1]
		df_a = all_sheets[label_a].loc[:,['id',var,'Consequence','p_value','proteins']]
		df_b = all_sheets[label_b].loc[:,['id',var,'Consequence','p_value','proteins']]
		data_full = df_a.merge(df_b, on=['id','Consequence','proteins'], suffixes=('_df1', '_df2'),how='inner')
		if data_full.empty:
			raise ValueError(f'Comparison {label_a} and {label_b} do not share "id/Consequence" elements')
		if args.given_protein :
			data_full = data_full.loc[data_full['proteins'].isin(set(args.given_protein)), :]
		if data_full.empty:
			raise ValueError(f'{" ".join(args.given_protein)} are not in the merge dataframe of {label_b} and {label_a}')
		if args.alpha_type=='Significance' :
			data_full['Significant']=[(row.p_value_df2<=args.p_thresh) or (row.p_value_df1<=args.p_thresh) for index, row in data_full.iterrows() ]
			data_full['alpha'] = [1 if val else 0.4 for val in data_full['Significant']]
			transparency_labels = [
			fr"$p\text{{-}}value \leq {args.p_thresh}$" + "\n for either element of comparison",
			fr"$p\text{{-}}value \geq {args.p_thresh}$" + "\n for both element of comparison"
			]
			transparency_labels=zip([1,0.4],transparency_labels)
		elif args.alpha_type=='all_little_bit':
			data_full['alpha'] = 0.5
			transparency_labels = None
		else :
			data_full['alpha'] = 1
			transparency_labels = None
		if args.element_plot == 'coding_only' :
			data_full=data_full.loc[data_full['Consequence'].isin(set(coding_seq)),:]
			for i in set(consequence_mapping.keys()) - set(coding_seq):
				consequence_mapping.pop(i, None)
		if sizey==1 and sizex ==1 :
			axed=axes
		elif sizey==1 or sizex==1:
			axed = axes[max(position[0], position[1])-1]
		else :
			axed=axes[position[0],position[1]]
		plot_scatter(axed, data_full, var, label_a, label_b, consequence_mapping=consequence_mapping if args.color_type == 'Consequence' else None, xlines=blines[label_a] if args.B_line else None, ylines=blines[label_b] if args.B_line else None)
	if sizey==1 and sizex ==1 :
		axed=axes
	elif sizey>1 and sizex==1:
		axed=axes[sizey-1]
	elif sizey==1 and sizex>1:
		axed=axes[0]
	else :
		axed=axes[0,sizey-1]
	add_legend_2(axed, consequence_mapping if args.color_type == 'Consequence' else None, transparency=transparency_labels)
	plt.tight_layout(rect=[0, 0, 0.90, 0.90],w_pad=1.5)
	if args.SQformat:
		mins, maxs = [], []
		for ax in axes.flatten():
			xmin, xmax = ax.get_xlim()
			ymin, ymax = ax.get_ylim()
			mins.append(xmin); maxs.append(xmax)
			mins.append(ymin); maxs.append(ymax)

		# Set uniform square limits
		mined, maxed = min(mins), max(maxs)
		lims = (mined, maxed)
		print(lims)
		for ax in axes.flatten():
			ax.set_xlim(lims)
			ax.set_ylim(lims)
			ax.set_box_aspect(1)
			print(lims)
			xt = ax.get_xticks()
			ax.set_yticks(xt)
	else :
		xmins, xmaxs, ymins, ymaxs = [], [], [], []
		for ax in axes.flatten():
			xmin, xmax = ax.get_xlim()
			ymin, ymax = ax.get_ylim()
			xmins.append(xmin); xmaxs.append(xmax)
			ymins.append(ymin); ymaxs.append(ymax)

		# Set uniform square limits
		xmin, xmax = min(xmins), max(xmaxs)
		ymin, ymax = min(ymins), max(ymaxs)
		span = max(xmax - xmin, ymax - ymin)
		xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
		xlim = (xmid - span / 2, xmid + span / 2)
		ylim = (ymid - span / 2, ymid + span / 2)
		for ax in axes.flatten():
			ax.set_xlim(xlim)
			ax.set_ylim(ylim)
	#for ax in axes.flatten():
#		ax.set_aspect('equal',adjustable = 'box') 
	if args.output :
		fig.savefig(args.output + '.pdf',format="pdf")
	else :
		fig.savefig(f'{comparison_label}_{"_".join(args.given_protein)}_{args.var}.pdf',format="pdf")


				
				

