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
import os
import re
import shlex
import warnings
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase


parser = argparse.ArgumentParser(description='This script produce a lolipop plot (with domain/feature annotations) base on a MaGeCK, a Variant Effect Prediction file (VEP), and bed style file')
parser.add_argument('-i', '--input', metavar='FILE', dest='excel_file', required=True, type=str, help='Summary excel file')
parser.add_argument('-b', '--lift_over', metavar='FILE', dest='lift_over', required=False, default=None, type=str, help='Chain bed-like file for liftover from 3D model to actual position. Optional if --structure is given and the mapping can be derived from it.')
parser.add_argument('-s', '--structure', metavar='FILE', dest='structure', required=False, default=None, nargs='+', help='One or more PDB/mmCIF structure files. When provided, the script attempts to derive the lift_over table from DBREF (PDB) or _struct_ref / _struct_ref_seq (mmCIF) records, validated against the mutation positions present in the input. Falls back to --lift_over if derivation fails.')
parser.add_argument('--Bio_threshold', dest='Biological_threshold', metavar='float',type=float, required=False, default=1, help='Biological threshold (two-sided)')
parser.add_argument('--duplicate_strategy', metavar='str', dest='duplicate_strategy', required=False,default='max' , type=str, choices={"max",'mean','median'}, help='which lfc score is to be reported (default : max)')
parser.add_argument('-p', '--Prob_Threshold', metavar='float',type=float, dest='p_thresh', required=False, default=1, help='P-Value Threshold representation (default no threshold)')
parser.add_argument('-R',"--Remove_sheet", required=False, dest='to_remove',  default =[], nargs='+', help="Sheets to remove")



def adjust_value(x, intervals):
	for ind, row in intervals.iterrows():
		if int(row['Start']) <= int(x) <= int(row['End']):
			return (int(x) - int(row['Start']) +1, row['Chain'])
	return (None,None)  # or np.nan if you prefer


def _parse_pdb_dbref(path):
	"""Parse DBREF records from a legacy PDB file.
	Returns list of dicts: {Chain, db_name, db_accession, db_code, db_beg, db_end}."""
	records = []
	with open(path) as fh:
		for line in fh:
			if not line.startswith('DBREF '):
				continue
			try:
				chain   = line[12].strip()
				db_name = line[26:32].strip()
				db_acc  = line[33:41].strip()
				db_code = line[42:54].strip()
				db_beg  = int(line[55:60].strip())
				db_end  = int(line[62:67].strip())
			except (ValueError, IndexError):
				continue
			if chain and db_beg and db_end:
				records.append({
					'Chain': chain, 'db_name': db_name,
					'db_accession': db_acc, 'db_code': db_code,
					'db_beg': db_beg, 'db_end': db_end,
				})
	return records


def _iter_cif_loop(lines, start):
	"""Yield header list and data rows for a mmCIF loop_ beginning at lines[start]='loop_'."""
	headers = []
	j = start + 1
	while j < len(lines) and lines[j].lstrip().startswith('_'):
		headers.append(lines[j].strip())
		j += 1
	rows = []
	while j < len(lines):
		s = lines[j].strip()
		if not s or s.startswith('#') or s.startswith('loop_') or s.startswith('data_') or s.startswith('_'):
			break
		try:
			toks = shlex.split(s)
		except ValueError:
			toks = s.split()
		rows.append(toks)
		j += 1
	return headers, rows, j


def _parse_cif_struct_ref(path):
	"""Parse _struct_ref and _struct_ref_seq from an mmCIF file.
	Returns list of dicts in the same normalized shape as _parse_pdb_dbref."""
	with open(path) as fh:
		lines = fh.readlines()
	refs = {}         # ref_id -> {db_name, db_accession, db_code}
	seqs = []         # [{ref_id, Chain, db_beg, db_end}]
	i = 0
	while i < len(lines):
		s = lines[i].strip()
		if s == 'loop_':
			headers, rows, nxt = _iter_cif_loop(lines, i)
			if headers:
				cat = headers[0].split('.')[0]
				cols = [h.split('.', 1)[1] for h in headers]
				if cat == '_struct_ref':
					for row in rows:
						if len(row) < len(cols): continue
						entry = dict(zip(cols, row))
						refs[entry.get('id')] = entry
				elif cat == '_struct_ref_seq':
					for row in rows:
						if len(row) < len(cols): continue
						seqs.append(dict(zip(cols, row)))
			i = nxt
			continue
		if s.startswith('_struct_ref.') or s.startswith('_struct_ref_seq.'):
			cat = s.split('.', 1)[0]
			block = {}
			while i < len(lines) and lines[i].lstrip().startswith(cat + '.'):
				parts = shlex.split(lines[i].strip()) if lines[i].strip() else []
				if len(parts) >= 2:
					block[parts[0].split('.', 1)[1]] = parts[1]
				elif len(parts) == 1 and i + 1 < len(lines):
					block[parts[0].split('.', 1)[1]] = lines[i + 1].strip()
					i += 1
				i += 1
			if cat == '_struct_ref':
				refs[block.get('id')] = block
			else:
				seqs.append(block)
			continue
		i += 1
	records = []
	for se in seqs:
		rid = se.get('ref_id')
		ref = refs.get(rid, {})
		chain = se.get('pdbx_strand_id') or se.get('strand_id') or ''
		try:
			db_beg = int(se.get('db_align_beg'))
			db_end = int(se.get('db_align_end'))
		except (TypeError, ValueError):
			continue
		records.append({
			'Chain': chain.strip(),
			'db_name': (ref.get('db_name') or '').strip(),
			'db_accession': (ref.get('pdbx_db_accession') or '').strip(),
			'db_code': (ref.get('db_code') or ref.get('pdbx_db_code') or '').strip(),
			'db_beg': db_beg, 'db_end': db_end,
		})
	return records


def _parse_structure_file(path):
	ext = os.path.splitext(path)[1].lower()
	if ext in ('.cif', '.mmcif'):
		return _parse_cif_struct_ref(path)
	if ext in ('.pdb', '.ent'):
		return _parse_pdb_dbref(path)
	# Unknown extension: try mmCIF first, then PDB.
	try:
		recs = _parse_cif_struct_ref(path)
		if recs: return recs
	except Exception:
		pass
	try:
		return _parse_pdb_dbref(path)
	except Exception:
		return []


def _collect_missense_positions(all_sheets):
	"""Return {protein: set(int positions)} over missense, non-control rows across all sheets."""
	out = {}
	for _, df in all_sheets.items():
		if not {'Controls', 'Consequence', 'Mutations', 'proteins'}.issubset(df.columns):
			continue
		sub = df.loc[df['Controls'].isnull() & (df['Consequence'] == 'missense'),
		             ['proteins', 'Mutations']]
		for _, row in sub.iterrows():
			prot = row['proteins']
			muts = row['Mutations']
			if not isinstance(muts, str):
				continue
			for m in muts.split(','):
				parts = m.split('_')
				if len(parts) < 2:
					continue
				try:
					pos = int(parts[1])
				except ValueError:
					continue
				out.setdefault(prot, set()).add(pos)
	return out


def derive_liftover_from_structures(structure_paths, positions_by_protein):
	"""Try to build a lift_over DataFrame from one or more structure files.

	Matches the input's `proteins` values (case-insensitive) against each structure
	record's UniProt accession or code. Accepts a mapping for a protein only if
	every observed mutation position falls inside at least one [db_beg, db_end]
	interval associated with that protein. Returns (df_or_None, diagnostics)."""
	records = []
	for p in structure_paths:
		try:
			records.extend(_parse_structure_file(p))
		except Exception as exc:
			print(f"[lift_over] Warning: failed to parse {p}: {exc}", file=sys.stderr)
	if not records:
		return None, "no DBREF / _struct_ref_seq records found in structure file(s)"

	rows, missing, mismatched = [], [], []
	for prot, positions in positions_by_protein.items():
		if prot is None or (isinstance(prot, float) and np.isnan(prot)):
			continue
		key = str(prot).strip().upper()
		hits = [r for r in records
		        if key in (r['db_accession'].upper(), r['db_code'].upper())]
		if not hits:
			missing.append(prot)
			continue
		uncovered = {pos for pos in positions
		             if not any(r['db_beg'] <= pos <= r['db_end'] for r in hits)}
		if uncovered:
			mismatched.append((prot, len(uncovered), len(positions)))
			continue
		for r in hits:
			rows.append({'Protein': str(prot), 'Chain': r['Chain'],
			             'Start': str(r['db_beg']), 'End': str(r['db_end'])})
	if missing or mismatched:
		msg = []
		if missing:
			msg.append(f"no structure match for: {', '.join(map(str, missing))}")
		if mismatched:
			frag = ', '.join(f"{p} ({u}/{t} positions outside [Start,End])"
			                 for p, u, t in mismatched)
			msg.append(f"coverage mismatch: {frag}")
		return None, '; '.join(msg)
	if not rows:
		return None, "no proteins matched between input and structure"
	return pd.DataFrame(rows, columns=['Protein', 'Chain', 'Start', 'End']), None


if __name__ == '__main__':
	args = parser.parse_args()
	try :
		all_sheets = pd.read_excel(args.excel_file, sheet_name=None)
	except:
		raise ValueError('Main input data (-i) does not have a supported format')
	val = []
	for k in args.to_remove:
		all_sheets.pop(k, None)

	lift_over = None
	derivation_error = None
	if args.structure:
		positions_by_protein = _collect_missense_positions(all_sheets)
		if not positions_by_protein:
			derivation_error = "could not collect mutation positions from input"
		else:
			lift_over, derivation_error = derive_liftover_from_structures(
				args.structure, positions_by_protein)
		if lift_over is not None:
			print(f"[lift_over] Derived {len(lift_over)} mapping(s) from "
			      f"{len(args.structure)} structure file(s).", file=sys.stderr)
		else:
			print(f"[lift_over] Could not derive mapping from structure(s): "
			      f"{derivation_error}", file=sys.stderr)

	if lift_over is None:
		if not args.lift_over:
			raise ValueError(
				"lift_over could not be derived from structure file(s)"
				+ (f" ({derivation_error})" if derivation_error else "")
				+ ". Provide one explicitly via -b/--lift_over "
				"(tab-separated file with columns: Protein, Chain, Start, End)."
			)
		lift_over = pd.read_csv(args.lift_over, sep='\t', header=0, dtype=str)

	if not set(['Protein','Chain','Start','End']).issubset(lift_over.columns):
		raise ValueError('Bed-like file does not contain the relevant columns')


	for label, data_full in all_sheets.items():
		if args.Biological_threshold <1 and args.Biological_threshold>0 :
			data_negative_control = data_full.loc[data_full['Controls']=='negative_control',:]
			lower_biology=np.quantile(pd.to_numeric(data_negative_control['lfc']),float(args.Biological_threshold) )
			upper_biology=np.quantile(pd.to_numeric(data_negative_control['lfc']),1-float(args.Biological_threshold) )
			data_full['Bio_Sig'] = (data_full['lfc'] < lower_biology) | (data_full['lfc'] > upper_biology)
		elif args.Biological_threshold ==1 :    
			data_full['Bio_Sig'] = True
		else :
			raise ValueError('Biological threshold must be between 0 and 1')
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
				data = data_full.loc[data_full['proteins']==protein,:].copy()
				for indexe, row in data.iterrows():
					try:
						_ = row['Mutations'].split(',')
					except AttributeError:
						print("Problematic value:", row['Mutations'], "of type", type(row['Mutations']), "row : ", row)
				data['Mutations']=[i.split(',') for i in data['Mutations']]
				data = data.explode('Mutations').dropna(subset=['Mutations'])
				data['Position'] = pd.Series([i.split('_')[1] for i in data['Mutations']], index=data.index, dtype='object')
				data.loc[:,'Chain'] = data['Position'].apply(lambda x: adjust_value(x, liftover)[1])
				data.loc[:,'Position'] = data['Position'].apply(lambda x: adjust_value(x, liftover)[0])
				data = data[data['Position'].notnull()]
				data['Sig']=[row['Bio_Sig'] and row['p_value']<=args.p_thresh for index, row in data.iterrows()]
				data_NS = data.loc[~data['Sig'],:].copy()
				data = data.loc[data['Sig'],:].copy()
				if args.duplicate_strategy == 'max':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].max()
				if args.duplicate_strategy == 'median':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].median()
				if args.duplicate_strategy == 'mean':
				  data = data.groupby(['Position','Chain'], as_index=False)['lfc'].mean()
				for _, row in data.sort_values('Position').iterrows():
					value = float(row['lfc'])
					f.write(f"\t/{row['Chain']}:{int(row['Position'])}\t{value:.2f}\n")
				position_sig=set(["_".join([row.Chain, str(int(row.Position))]) for ind, row in data.iterrows()])
				position_NS =["_".join([row.Chain, str(int(row.Position))]) for ind, row in data_NS.iterrows()]
				position_NS = natsorted(list(set(position_NS) - position_sig))
				NS=[i.split('_') for i in position_NS]
				for listed in NS:
					g.write(f"\t/{listed[0]}:{listed[1]}\t1\n")
