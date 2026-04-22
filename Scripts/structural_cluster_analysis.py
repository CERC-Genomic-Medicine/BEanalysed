#!/usr/bin/env python3
"""
structural_cluster_analysis.py

Consume the outputs of ChimeraX_scoring.py (``*_Sig.defattr`` and
``*_nonSig.defattr``), the original summary Excel, and the 3D structure
(PDB or mmCIF) used to build those files, then:

  1. Cluster the significant residues in 3D (DBSCAN on Cα coordinates) and
     compute complementary clustering statistics (mean/median pairwise
     Cα-Cα distance, fraction of pairs within a contact threshold, total
     contact pairs, number of DBSCAN clusters, largest-cluster size).

  2. Run a bootstrap / Monte-Carlo randomization test: re-draw, ``n_boot``
     times, the same number of residues uniformly from the *targeted*
     background (significant ∪ non-significant residues written by
     ChimeraX_scoring), recompute the statistics, and report empirical
     one-sided p-values against the observed value.

Best-practice notes (and why):

  * The null is drawn from the residues that could have been hit by the
    library, not from all residues in the structure.  This is the
    convention used by HotSpot3D, HotMAPS, CLUMPS, and mutation3D — it
    removes bias from unequal guide-library coverage and from missing
    residues in the model.

  * Several statistics are reported because any single one is easy to game
    on small samples.  A cluster is convincing when distance statistics
    drop AND contact counts rise AND DBSCAN recovers a compact group.

  * Each (sheet, protein, chain) is analysed independently so that hits on
    different chains or datasets do not bleed into the same null.

  * Small-sample fix-up: the empirical p-value uses ``(k+1)/(n+1)`` so it
    can never be exactly 0, which keeps downstream log-transforms safe.

  * Cα distances are used to match the convention of the tools above and
    to keep runtime low; switching to closest-atom distances would change
    absolute numbers but rarely qualitative conclusions.

Usage:

  structural_cluster_analysis.py \\
      -d chimera_out_dir \\
      -s model.pdb \\
      -i summary.xlsx \\
      -o cluster_report \\
      --n_boot 10000 --contact_thresh 8 --dbscan_eps 10 --dbscan_min 3
"""

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Structure-file parsing (Cα coordinates only; dependency-free).
# ---------------------------------------------------------------------------

def read_ca_coords(path):
	"""Return {(chain, auth_seq): np.ndarray([x,y,z])} of Cα atoms
	from the first MODEL of a PDB or mmCIF file."""
	ext = os.path.splitext(path)[1].lower()
	if ext in ('.cif', '.mmcif'):
		return _read_ca_cif(path)
	return _read_ca_pdb(path)


def _read_ca_pdb(path):
	coords = {}
	in_first_model = True
	with open(path) as fh:
		for line in fh:
			if line.startswith('ENDMDL'):
				in_first_model = False
				continue
			if line.startswith('MODEL'):
				if coords:
					break
				in_first_model = True
				continue
			if not in_first_model:
				continue
			if not line.startswith(('ATOM  ', 'HETATM')):
				continue
			if line[12:16].strip() != 'CA':
				continue
			altloc = line[16]
			if altloc not in (' ', 'A'):
				continue
			chain = line[21].strip()
			try:
				resnum = int(line[22:26].strip())
				x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
			except ValueError:
				continue
			coords.setdefault((chain, resnum), np.array([x, y, z]))
	return coords


def _read_ca_cif(path):
	with open(path) as fh:
		lines = fh.readlines()
	coords = {}
	i = 0
	while i < len(lines):
		if lines[i].strip() == 'loop_':
			j = i + 1
			headers = []
			while j < len(lines) and lines[j].lstrip().startswith('_atom_site.'):
				headers.append(lines[j].strip().split('.', 1)[1])
				j += 1
			if headers:
				col = {h: idx for idx, h in enumerate(headers)}
				asym = 'auth_asym_id' if 'auth_asym_id' in col else 'label_asym_id'
				seq  = 'auth_seq_id'  if 'auth_seq_id'  in col else 'label_seq_id'
				required = ('label_atom_id', asym, seq, 'Cartn_x', 'Cartn_y', 'Cartn_z')
				if all(k in col for k in required):
					first_model = None
					while j < len(lines):
						s = lines[j].strip()
						if not s or s.startswith('#') or s.startswith('loop_') \
						   or s.startswith('_') or s.startswith('data_'):
							break
						parts = s.split()
						if len(parts) < len(headers):
							j += 1; continue
						if 'pdbx_PDB_model_num' in col:
							m = parts[col['pdbx_PDB_model_num']]
							if first_model is None: first_model = m
							elif m != first_model: break
						if parts[col['label_atom_id']].strip('"') != 'CA':
							j += 1; continue
						try:
							chain = parts[col[asym]]
							resnum = int(parts[col[seq]])
							x = float(parts[col['Cartn_x']])
							y = float(parts[col['Cartn_y']])
							z = float(parts[col['Cartn_z']])
						except ValueError:
							j += 1; continue
						coords.setdefault((chain, resnum), np.array([x, y, z]))
						j += 1
			i = j
			continue
		i += 1
	return coords


# ---------------------------------------------------------------------------
# .defattr parsing (ChimeraX_scoring outputs).
# ---------------------------------------------------------------------------

_DEFATTR_RE = re.compile(r'^\s*/([^:\s]+):(-?\d+)\s+(\S+)\s*$')


def parse_defattr(path):
	"""Return list of (chain, resnum, value) from a ChimeraX .defattr file."""
	out = []
	with open(path) as fh:
		for line in fh:
			m = _DEFATTR_RE.match(line)
			if not m: continue
			try:
				out.append((m.group(1), int(m.group(2)), float(m.group(3))))
			except ValueError:
				continue
	return out


def discover_defattr_pairs(paths):
	"""From files or directories, return {label: {'sig': path, 'ns': path}}."""
	files = []
	for p in paths:
		if os.path.isdir(p):
			files.extend(sorted(glob.glob(os.path.join(p, '*.defattr'))))
		elif os.path.isfile(p):
			files.append(p)
	pairs = {}
	for f in files:
		b = os.path.basename(f)
		if b.endswith('_Sig.defattr'):
			pairs.setdefault(b[:-len('_Sig.defattr')], {})['sig'] = f
		elif b.endswith('_nonSig.defattr'):
			pairs.setdefault(b[:-len('_nonSig.defattr')], {})['ns'] = f
	return pairs


# ---------------------------------------------------------------------------
# Clustering statistics.
# ---------------------------------------------------------------------------

def _pairwise_distances(points):
	if len(points) < 2:
		return np.array([])
	arr = np.asarray(points, dtype=float)
	diff = arr[:, None, :] - arr[None, :, :]
	d = np.sqrt((diff * diff).sum(-1))
	iu = np.triu_indices(len(arr), k=1)
	return d[iu]


def clustering_stats(points, contact_thresh):
	d = _pairwise_distances(points)
	if d.size == 0:
		return {'n': len(points), 'mean_d': np.nan, 'median_d': np.nan,
		        'frac_contact': np.nan, 'n_contact_pairs': 0}
	hits = d <= contact_thresh
	return {
		'n': len(points),
		'mean_d': float(d.mean()),
		'median_d': float(np.median(d)),
		'frac_contact': float(hits.mean()),
		'n_contact_pairs': int(hits.sum()),
	}


def dbscan(points, eps, min_samples):
	"""Return list of clusters (each a list of indices into `points`).
	Points not in any cluster are reported as noise and omitted."""
	n = len(points)
	if n == 0: return []
	arr = np.asarray(points, dtype=float)
	labels = np.full(n, -1, dtype=int)
	visited = np.zeros(n, dtype=bool)
	cid = -1

	def neighbors(i):
		d = np.sqrt(((arr - arr[i]) ** 2).sum(-1))
		return list(np.where(d <= eps)[0])

	for i in range(n):
		if visited[i]: continue
		visited[i] = True
		nbrs = neighbors(i)
		if len(nbrs) < min_samples:
			continue
		cid += 1
		labels[i] = cid
		queue = list(nbrs)
		seen = set(nbrs)
		k = 0
		while k < len(queue):
			j = queue[k]; k += 1
			if not visited[j]:
				visited[j] = True
				njb = neighbors(j)
				if len(njb) >= min_samples:
					for x in njb:
						if x not in seen:
							seen.add(x); queue.append(x)
			if labels[j] == -1:
				labels[j] = cid
	return [[int(i) for i in np.where(labels == c)[0]] for c in range(cid + 1)]


# ---------------------------------------------------------------------------
# Bootstrap / randomization.
# ---------------------------------------------------------------------------

def bootstrap_null(background_points, k, n_boot, contact_thresh, rng):
	"""Draw k residues without replacement from `background_points` n_boot
	times and return dict of np.arrays with the same statistic names as
	clustering_stats (excluding 'n')."""
	bg = np.asarray(background_points, dtype=float)
	stats = {'mean_d': [], 'median_d': [], 'frac_contact': [], 'n_contact_pairs': []}
	if len(bg) < k or k < 2:
		return {k_: np.array([]) for k_ in stats}
	idx = np.arange(len(bg))
	for _ in range(n_boot):
		sel = rng.choice(idx, size=k, replace=False)
		s = clustering_stats(bg[sel], contact_thresh)
		for key in stats:
			stats[key].append(s[key])
	return {k_: np.asarray(v, dtype=float) for k_, v in stats.items()}


def empirical_p(null, observed, tail):
	"""One-sided empirical p using the (k+1)/(n+1) correction."""
	null = np.asarray(null, dtype=float)
	null = null[~np.isnan(null)]
	if null.size == 0 or observed is None or (isinstance(observed, float) and np.isnan(observed)):
		return float('nan')
	if tail == 'less':
		k = int(np.sum(null <= observed))
	else:
		k = int(np.sum(null >= observed))
	return float((k + 1) / (null.size + 1))


# ---------------------------------------------------------------------------
# Summary-file annotation (optional).
# ---------------------------------------------------------------------------

def summary_annotations(summary_path, sheet_labels):
	"""Return {sheet: DataFrame(proteins, Consequence_Detail, Mutations, lfc, p_value)}.
	Empty dict on failure; the analysis does not depend on this."""
	if not summary_path:
		return {}
	try:
		book = pd.read_excel(summary_path, sheet_name=None)
	except Exception as exc:
		print(f"[annot] Could not read summary {summary_path}: {exc}", file=sys.stderr)
		return {}
	out = {}
	keep = ['proteins', 'Consequence', 'Consequence_Detail',
	        'Mutations', 'Mutations_full', 'lfc', 'p_value', 'fdr']
	for lbl in sheet_labels:
		df = book.get(lbl)
		if df is None: continue
		out[lbl] = df[[c for c in keep if c in df.columns]].copy()
	return out


# ---------------------------------------------------------------------------
# Pipeline.
# ---------------------------------------------------------------------------

def analyse_group(label, chain, sig_entries, ns_entries, coords,
                  contact_thresh, dbscan_eps, dbscan_min, n_boot, rng):
	"""Analyse one (sheet, chain) group.  Returns dict with stats + clusters
	+ bootstrap p-values, or None if there is not enough data."""
	sig_coords, sig_res, sig_val = [], [], []
	for ch, res, val in sig_entries:
		if ch != chain: continue
		c = coords.get((ch, res))
		if c is None: continue
		sig_coords.append(c); sig_res.append(res); sig_val.append(val)

	bg_coords, bg_res = [], []
	seen = set()
	for entries in (sig_entries, ns_entries):
		for ch, res, _ in entries:
			if ch != chain or (ch, res) in seen: continue
			c = coords.get((ch, res))
			if c is None: continue
			seen.add((ch, res))
			bg_coords.append(c); bg_res.append(res)

	if len(sig_coords) < 2 or len(bg_coords) < len(sig_coords) + 1:
		return None

	obs = clustering_stats(sig_coords, contact_thresh)
	clusters_idx = dbscan(sig_coords, dbscan_eps, dbscan_min)
	clusters = [{
		'members': sorted(sig_res[i] for i in members),
		'size': len(members),
		'mean_lfc': float(np.mean([sig_val[i] for i in members])),
		'centroid': [float(x) for x in np.mean([sig_coords[i] for i in members], axis=0)],
	} for members in clusters_idx]

	null = bootstrap_null(bg_coords, len(sig_coords), n_boot, contact_thresh, rng)

	# DBSCAN-based null (number of clusters, largest cluster size).
	null_n_clusters, null_largest = [], []
	if len(bg_coords) >= len(sig_coords):
		bg_arr = np.asarray(bg_coords)
		idx_all = np.arange(len(bg_arr))
		for _ in range(n_boot):
			sel = rng.choice(idx_all, size=len(sig_coords), replace=False)
			cs = dbscan(bg_arr[sel], dbscan_eps, dbscan_min)
			null_n_clusters.append(len(cs))
			null_largest.append(max((len(c) for c in cs), default=0))
	null_n_clusters = np.asarray(null_n_clusters, dtype=float)
	null_largest    = np.asarray(null_largest,    dtype=float)

	n_obs_clusters = len(clusters)
	largest_obs    = max((c['size'] for c in clusters), default=0)

	return {
		'sheet': label,
		'chain': chain,
		'n_sig': len(sig_coords),
		'n_background': len(bg_coords),
		'mean_d': obs['mean_d'],
		'p_mean_d': empirical_p(null['mean_d'], obs['mean_d'], 'less'),
		'median_d': obs['median_d'],
		'p_median_d': empirical_p(null['median_d'], obs['median_d'], 'less'),
		'frac_contact': obs['frac_contact'],
		'p_frac_contact': empirical_p(null['frac_contact'], obs['frac_contact'], 'greater'),
		'n_contact_pairs': obs['n_contact_pairs'],
		'p_n_contact_pairs': empirical_p(null['n_contact_pairs'], obs['n_contact_pairs'], 'greater'),
		'n_dbscan_clusters': n_obs_clusters,
		'p_n_dbscan_clusters': empirical_p(null_n_clusters, n_obs_clusters, 'greater'),
		'largest_cluster_size': largest_obs,
		'p_largest_cluster_size': empirical_p(null_largest, largest_obs, 'greater'),
		'clusters': clusters,
		'null': {k: v.tolist() for k, v in {
			'mean_d': null['mean_d'], 'median_d': null['median_d'],
			'frac_contact': null['frac_contact'],
			'n_contact_pairs': null['n_contact_pairs'],
			'n_dbscan_clusters': null_n_clusters,
			'largest_cluster_size': null_largest,
		}.items()},
	}


def write_cluster_defattr(out_path, label, clusters_by_chain):
	"""Write a ChimeraX .defattr mapping each clustered residue to its cluster id.
	Emit the id as a float so ChimeraX parses the attribute numerically — integer
	strings end up as str on the residue object, which breaks ``::cluster_id>0``
	selections."""
	with open(out_path, 'w') as fh:
		fh.write("#\n#\n")
		fh.write("attribute: cluster_id\n")
		fh.write("recipient: residues\n")
		counter = 0
		for chain in sorted(clusters_by_chain):
			for cluster in clusters_by_chain[chain]:
				counter += 1
				for res in cluster['members']:
					fh.write(f"\t/{chain}:{int(res)}\t{float(counter):.1f}\n")


def main():
	ap = argparse.ArgumentParser(
		description="3D clustering + bootstrap test on ChimeraX_scoring outputs.")
	ap.add_argument('-d', '--defattr', required=True, nargs='+',
		help='Directory (or individual files) containing *_Sig.defattr and '
		     '*_nonSig.defattr produced by ChimeraX_scoring.py.')
	ap.add_argument('-s', '--structure', required=True,
		help='PDB or mmCIF file of the 3D model used by ChimeraX_scoring.')
	ap.add_argument('-i', '--input', required=False, default=None,
		help='Original summary Excel (optional; used to enrich cluster reports '
		     'with mutation/consequence annotations).')
	ap.add_argument('-o', '--outdir', required=True,
		help='Output directory; created if missing.')
	ap.add_argument('--n_boot', type=int, default=10000,
		help='Bootstrap / randomization iterations (default 10000).')
	ap.add_argument('--contact_thresh', type=float, default=8.0,
		help='Cα-Cα contact threshold in Å (default 8).')
	ap.add_argument('--dbscan_eps', type=float, default=10.0,
		help='DBSCAN neighborhood radius in Å (default 10).')
	ap.add_argument('--dbscan_min', type=int, default=3,
		help='DBSCAN min_samples (default 3).')
	ap.add_argument('--seed', type=int, default=0,
		help='RNG seed for reproducibility (default 0).')
	args = ap.parse_args()

	os.makedirs(args.outdir, exist_ok=True)
	rng = np.random.default_rng(args.seed)

	coords = read_ca_coords(args.structure)
	if not coords:
		raise SystemExit(f"No Cα atoms parsed from {args.structure}")
	print(f"[structure] {len(coords)} Cα atoms across "
	      f"{len(set(c for c, _ in coords))} chain(s)", file=sys.stderr)

	pairs = discover_defattr_pairs(args.defattr)
	if not pairs:
		raise SystemExit("No *_Sig.defattr / *_nonSig.defattr files found.")
	print(f"[defattr] {len(pairs)} sheet(s): {', '.join(sorted(pairs))}", file=sys.stderr)

	annotations = summary_annotations(args.input, list(pairs))

	results = []
	clusters_by_label = {}
	for label, files in sorted(pairs.items()):
		sig = parse_defattr(files['sig']) if 'sig' in files else []
		ns  = parse_defattr(files['ns'])  if 'ns'  in files else []
		if not sig:
			print(f"[{label}] no significant residues; skipping.", file=sys.stderr)
			continue
		chains = sorted({c for c, _, _ in sig})
		clusters_by_label[label] = {}
		for chain in chains:
			res = analyse_group(label, chain, sig, ns, coords,
			                    args.contact_thresh, args.dbscan_eps,
			                    args.dbscan_min, args.n_boot, rng)
			if res is None:
				print(f"[{label}/{chain}] insufficient data; skipping.", file=sys.stderr)
				continue
			# Attach per-cluster mutation annotation from summary if available.
			annot = annotations.get(label)
			if annot is not None and 'Mutations' in annot.columns \
			   and 'proteins' in annot.columns:
				for cluster in res['clusters']:
					# We only know chain+residue_in_chain here; attach best-effort
					# via position offsets later. Skip unless user inspects JSON.
					cluster['mutations'] = []
			results.append(res)
			clusters_by_label[label][chain] = res['clusters']
		if clusters_by_label[label]:
			write_cluster_defattr(
				os.path.join(args.outdir, f"{label}_clusters.defattr"),
				label, clusters_by_label[label])

	# Compact TSV summary.
	tsv_cols = ['sheet', 'chain', 'n_sig', 'n_background',
	            'mean_d', 'p_mean_d', 'median_d', 'p_median_d',
	            'frac_contact', 'p_frac_contact',
	            'n_contact_pairs', 'p_n_contact_pairs',
	            'n_dbscan_clusters', 'p_n_dbscan_clusters',
	            'largest_cluster_size', 'p_largest_cluster_size']
	df = pd.DataFrame([{k: r[k] for k in tsv_cols} for r in results])
	df.to_csv(os.path.join(args.outdir, 'cluster_stats.tsv'),
	          sep='\t', index=False, float_format='%.4f')

	# Full JSON (stats + clusters + null distributions).
	with open(os.path.join(args.outdir, 'clusters.json'), 'w') as fh:
		json.dump(results, fh, indent=2, default=float)

	print(f"[done] wrote cluster_stats.tsv and clusters.json "
	      f"({len(results)} group(s)) to {args.outdir}", file=sys.stderr)


if __name__ == '__main__':
	main()
