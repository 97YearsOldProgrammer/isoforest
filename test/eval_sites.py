#!/usr/bin/env python3

"""
Splice site evaluation.

The key question for isoforest: does the HMM + Random Forest pipeline find
the right splice sites? The upstream isoforms project uses ALL GT/AG sites
(exhaustive). isoforest uses HMM posteriors + RF classification to narrow
the search space before scoring. This module measures how good that filtering is.

Evaluates:
  1. Site-level precision/recall/F1: did we find the real donors and acceptors?
  2. Hint filter comparison: which cutoff method (gap, topk, percentile, etc.)
     gives the best site selection?
  3. ROC/AUC: continuous evaluation of HMM posterior as a classifier
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))
import isoform

###################################
## EXTRACT REAL SITES FROM GFF   ##
###################################

def get_real_sites(gff, source=None):
	"""extract real donor/acceptor positions from annotation GFF

	Donors = intron start positions, Acceptors = intron end positions.
	Positions are 1-based (GFF convention).
	"""
	dons = set()
	accs = set()
	with open(gff) as fp:
		for line in fp:
			if line.startswith('#'): continue
			f = line.split('\t') if '\t' in line else line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			if f[6] != '+': continue
			if source is not None and f[1] != source: continue
			dons.add(int(f[3]))
			accs.add(int(f[4]))
	return dons, accs

def get_real_sites_0based(gff, source=None):
	"""same as get_real_sites but 0-based (for isoform.gtag_sites compatibility)"""
	dons, accs = get_real_sites(gff, source)
	return {d - 1 for d in dons}, {a - 1 for a in accs}

##############################
## PRECISION / RECALL / F1  ##
##############################

def precision_recall(pred, real):
	"""set-based precision, recall, F1 for splice site positions"""
	pred = set(pred)
	real = set(real)
	tp = len(pred & real)
	fp = len(pred - real)
	fn = len(real - pred)
	prec = tp / (tp + fp) if tp + fp > 0 else 0
	rec  = tp / (tp + fn) if tp + fn > 0 else 0
	f1   = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0
	return {
		'precision': prec,
		'recall':    rec,
		'f1':        f1,
		'tp': tp, 'fp': fp, 'fn': fn,
	}

def eval_sites(pred_dons, pred_accs, real_dons, real_accs):
	"""evaluate both donor and acceptor site predictions"""
	return {
		'donors':    precision_recall(pred_dons, real_dons),
		'acceptors': precision_recall(pred_accs, real_accs),
	}

##################################
## HINT FILTER METHOD COMPARISON ##
##################################

def eval_hint_methods(hmm_dons, hmm_accs, real_dons, real_accs, methods):
	"""compare multiple HMM hint filtering strategies

	hmm_dons, hmm_accs: list of (position, posterior_score) from HMM
	real_dons, real_accs: set of real positions
	methods: dict of {name: filter_function}
	    each filter_function takes [(pos, score), ...] and returns [pos, ...]

	Returns dict of {method_name: {donors: metrics, acceptors: metrics}}
	"""
	results = {}
	for name, fn in methods.items():
		pred_d = fn(hmm_dons)
		pred_a = fn(hmm_accs)
		results[name] = eval_sites(pred_d, pred_a, real_dons, real_accs)
	return results

#########################
## ROC / AUC           ##
#########################

def site_roc(scores, real_positions, n_points=100):
	"""compute ROC curve for HMM posteriors as splice site classifier

	scores: list of (position, posterior_score)
	real_positions: set of true positive positions
	n_points: number of threshold points to evaluate

	Returns (thresholds, tpr, fpr, auc)
	"""
	if not scores or not real_positions: return [], [], [], 0

	real = set(real_positions)
	vals = sorted([s for _, s in scores])
	lo = vals[0]
	hi = vals[-1]
	if lo == hi: return [lo], [1.0], [1.0], 0.5

	thresholds = np.linspace(lo, hi, n_points).tolist()

	tpr_list = []
	fpr_list = []

	total_pos = len(real)
	total_neg = len(scores) - total_pos

	for t in thresholds:
		pred = {pos for pos, val in scores if val >= t}
		tp = len(pred & real)
		fp = len(pred - real)
		tpr = tp / total_pos if total_pos > 0 else 0
		fpr = fp / total_neg if total_neg > 0 else 0
		tpr_list.append(tpr)
		fpr_list.append(fpr)

	# AUC via trapezoidal rule (fpr ascending)
	pairs = sorted(zip(fpr_list, tpr_list))
	auc = 0
	for i in range(1, len(pairs)):
		dx = pairs[i][0] - pairs[i-1][0]
		dy = (pairs[i][1] + pairs[i-1][1]) / 2
		auc += dx * dy

	return thresholds, tpr_list, fpr_list, auc

############################
## BASELINE COMPARISON    ##
############################

def eval_gtag_baseline(seq, real_dons, real_accs, flank=99, minex=25):
	"""evaluate the naive GT/AG baseline (all sites = predicted)

	This is what geniso does: use every GT as a donor, every AG as an acceptor.
	Measures how much noise the HMM/RF needs to filter out.
	"""
	all_dons, all_accs = isoform.gtag_sites(seq, flank, minex)
	return {
		'donors':       precision_recall(all_dons, real_dons),
		'acceptors':    precision_recall(all_accs, real_accs),
		'total_gt':     len(all_dons),
		'total_ag':     len(all_accs),
		'real_donors':  len(real_dons),
		'real_acceptors': len(real_accs),
		'noise_ratio_don': len(all_dons) / len(real_dons) if real_dons else 0,
		'noise_ratio_acc': len(all_accs) / len(real_accs) if real_accs else 0,
	}

##########
## CLI  ##
##########

if __name__ == '__main__':
	import argparse
	import json

	parser = argparse.ArgumentParser(description='evaluate splice site predictions')
	parser.add_argument('fasta', help='sequence fasta file')
	parser.add_argument('gff', help='real annotation GFF')
	parser.add_argument('--source', default=None, help='annotation source filter')
	parser.add_argument('--flank', type=int, default=99)
	arg = parser.parse_args()

	name, seq = next(isoform.read_fasta(arg.fasta))
	real_dons, real_accs = get_real_sites_0based(arg.gff, arg.source)

	# baseline: all GT/AG sites
	baseline = eval_gtag_baseline(seq, real_dons, real_accs, flank=arg.flank)

	print('GT/AG baseline:')
	print(f'  total GT sites:  {baseline["total_gt"]}')
	print(f'  total AG sites:  {baseline["total_ag"]}')
	print(f'  real donors:     {baseline["real_donors"]}')
	print(f'  real acceptors:  {baseline["real_acceptors"]}')
	print(f'  noise ratio don: {baseline["noise_ratio_don"]:.1f}x')
	print(f'  noise ratio acc: {baseline["noise_ratio_acc"]:.1f}x')
	print(f'  donor recall:    {baseline["donors"]["recall"]:.4f}')
	print(f'  acceptor recall: {baseline["acceptors"]["recall"]:.4f}')
	print(f'  donor precision: {baseline["donors"]["precision"]:.4f}')
	print(f'  acceptor precision: {baseline["acceptors"]["precision"]:.4f}')
