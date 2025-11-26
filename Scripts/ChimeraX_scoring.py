#!/usr/bin/env python3


'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

Goal :

To produce ChimeraX attribute file.


'''

import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import argparse
from scipy.stats import rankdata
from scipy.stats import binomtest
import seaborn as sns
import sys
import re
import warnings
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase


parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-i', '--input', metavar='FILE', dest='excel_file', required=True, type=str, help='Summary excel file')
parser.add_argument('-b', '--lift_over', metavar='FILE', dest='lift_over', required=True, type=str, help='Chain bed-like file for liftove from 3D model to actual position')
parser.add_argument('--Bio_threshold', dest='Biological_threshold', metavar='float',type=float, required=False, default=1, help='Biological threshold (two-sided)')
parser.add_argument('--duplicate_strategy', metavar='str', dest='duplicate_strategy', required=False,default='median' , type=str, choices={"max",'mean','median'}, help='which lfc score is to be reported (default : median)')
parser.add_argument('-p', '--Prob_Threshold', metavar='float',type=float, dest='p_thresh', required=False, default=1, help='P-Value Threshold representation (default no threshold)')
parser.add_argument('-R',"--Remove_sheet", required=False, dest='to_remove',  default =[], nargs='+', help="Sheets to remove")



def adjust_value(x, intervals):
	for ind, row in intervals.iterrows():
		if int(row['Start']) <= int(x) <= int(row['End']):
			return (int(x) - int(row['Start']) +1, row['Chain'])
	return (None,None)  # or np.nan if you prefer


if __name__ == '__main__':
	args = parser.parse_args()
	try :
		all_sheets = pd.read_excel(args.excel_file, sheet_name=None)
	except: 
		raise ValueError('Main input data (-i) does not have a supported format')
	lift_over = pd.read_csv(args.lift_over, sep='\t', header=0,dtype=str)
	val = []
	for k in args.to_remove:
		all_sheets.pop(k, None)
	if not set(['Protein','Chain','Start','End']).issubset(lift_over.columns):
		raise ValueError('Bed-like file does not contain the relevant columns')
	for label, data_full in all_sheets.items():
		data_full=data_full.loc[data_full['Controls'].isnull(),:]
		data_full=data_full.loc[[i=='missense' for i in data_full['Consequence']],:]
		proteins=set(data_full['proteins'])
		with open(f"{label}_Sig.defattr", 'w') as f, open(f"{label}_nonSig.defattr", 'w') as g:
			f.write("#\n#\n")
			f.write(f'attribute: {args.duplicate_strategy}\n')
			f.write("recipient: residues\n")
			g.write("#\n#\n")
			g.write(f'attribute: NS\n')
			g.write("recipient: residues\n")
			for protein in proteins :
				liftover=lift_over.loc[lift_over['Protein']==protein ,:]
				data = data_full.loc[data_full['proteins']==protein,:]
				for indexe, row in data.iterrows():
					try:
						_ = row['Mutations'].split(',')
					except AttributeError:
						print("Problematic value:", row['Mutations'], "of type", type(row['Mutations']), "row : ", row)
				data['Mutations']=[i.split(',') for i in data['Mutations']]
				data = data.explode('Mutations').dropna(subset=['Mutations'])
				data.loc[:,'Position'] = [i.split('_')[1] for i in data['Mutations']]
				data.loc[:,'Chain'] = data['Position'].apply(lambda x: adjust_value(x, liftover)[1])
				data.loc[:,'Position'] = data['Position'].apply(lambda x: adjust_value(x, liftover)[0])
				data = data[data['Position'].notnull()]
				data['Sig']=[row['Bio_p-value']<=args.Biological_threshold and row['p_value']<=args.p_thresh for index, row in data.iterrows()]
				data_NS = data.loc[~data['Sig'],:]
				data = data.loc[data['Sig'],:]
				if args.duplicate_strategy == 'max':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].max()
				if args.duplicate_strategy == 'median':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].median()
				if args.duplicate_strategy == 'mean':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].mean()
				for _, row in data.sort_values('Position').iterrows():
					value = float(row['lfc'])
					f.write(f"\t/{row['Chain']}:{int(row['Position'])}\t{value:.2f}\n")
					val.append(value)
				position_sig=set(["_".join([row.Chain, str(int(row.Position))]) for ind, row in data.iterrows()])
				position_NS =["_".join([row.Chain, str(int(row.Position))]) for ind, row in data_NS.iterrows()]
				position_NS = natsorted(list(set(position_NS) - position_sig))
				NS=[i.split('_') for i in position_NS]
				for listed in NS:
					g.write(f"\t/{listed[0]}:{listed[1]}\t1\n")

	print(f'max_value : {max(val)}, min_value {min(val)}')
