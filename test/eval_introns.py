#!/usr/bin/env python3

"""
Intron distribution evaluation.

The core question from the isoforms research: how well does the predicted
intron frequency distribution match the real one (from annotation/RNA-seq)?

Upstream project (isoforms) uses:
  - cmpiso: compare two GFFs via expdiff (Manhattan distance on intron probs)
  - optiso: genetic algorithm to minimize Manhattan distance via weight tuning
  - hintersection.py: histogram of intersection distances across all genes

This module provides:
  1. Per-gene intron distribution comparison (predicted vs real)
  2. Binary intron detection (did we find the right introns at all?)
  3. Batch evaluation across a smallgenes directory
  4. Summary statistics with numpy
"""

import copy
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))
import isoform

##############################
## INTRON EXTRACTION (CORE) ##
##############################

def get_introns(gff, source=None):
	"""extract normalized intron distribution from GFF

	When source is specified, only introns from that source are used.
	When all scores are 0 or '.', falls back to uniform distribution.
	This handles D. mel NCBI annotations (no expression scores).
	"""
	introns = {}
	total = 0
	with open(gff) as fp:
		for line in fp:
			if line.startswith('#'): continue
			f = line.split('\t') if '\t' in line else line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			if f[6] != '+': continue
			if source is not None and f[1] != source: continue
			beg = int(f[3])
			end = int(f[4])
			score = 0 if f[5] == '.' else float(f[5])
			if (beg, end) not in introns: introns[(beg, end)] = 0
			introns[(beg, end)] += score
			total += score

	if not introns: return {}

	# if no scores available, uniform distribution
	if total == 0:
		n = len(introns)
		return {k: 1.0 / n for k in introns}

	return {k: v / total for k, v in introns.items()}

def get_introns_by_source(gff):
	"""extract intron distributions from each source separately"""
	sources = set()
	with open(gff) as fp:
		for line in fp:
			if line.startswith('#'): continue
			f = line.split('\t') if '\t' in line else line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			sources.add(f[1])
	return {src: get_introns(gff, source=src) for src in sources}

#########################
## DISTANCE METRICS    ##
## (from isoforms lib) ##
#########################

def manhattan(p, q):
	d = 0
	for pi, qi in zip(p, q):
		d += abs(pi - qi)
	return d

def intersection(p, q):
	d = 0
	for pi, qi in zip(p, q):
		d += min(pi, qi)
	return 1 - d

def dkl(p, q):
	d = 0
	for pi, qi in zip(p, q):
		if pi > 0 and qi > 0:
			d += pi * math.log2(pi / qi)
	return d

def chebyshev(p, q):
	d = 0
	for pi, qi in zip(p, q):
		if abs(pi - qi) > d: d = abs(pi - qi)
	return d

#############################
## PER-GENE COMPARISON     ##
#############################

def compare_introns(pred, real):
	"""compare two intron distributions, returns all metrics + per-intron detail

	pred, real: dict of {(beg, end): probability}
	"""
	if not pred or not real: return None

	# union of all introns, fill missing with 0
	p = copy.deepcopy(pred)
	r = copy.deepcopy(real)
	for k in p:
		if k not in r: r[k] = 0
	for k in r:
		if k not in p: p[k] = 0

	keys = sorted(p.keys())
	pv = [p[k] for k in keys]
	rv = [r[k] for k in keys]

	details = [(k, p[k], r[k]) for k in keys]

	return {
		'manhattan':    manhattan(pv, rv),
		'intersection': intersection(pv, rv),
		'dkl':          dkl(rv, pv),   # KL(real || pred)
		'chebyshev':    chebyshev(pv, rv),
		'n_pred':       len(pred),
		'n_real':       len(real),
		'details':      details,
	}

def intron_f1(pred, real):
	"""binary intron detection — did we find the right introns at all?

	pred, real: dict of {(beg, end): probability} or set of (beg, end)
	Returns precision, recall, F1 on the intron SET (ignores probabilities).
	"""
	ps = set(pred.keys()) if isinstance(pred, dict) else set(pred)
	rs = set(real.keys()) if isinstance(real, dict) else set(real)
	tp = len(ps & rs)
	fp = len(ps - rs)
	fn = len(rs - ps)
	prec = tp / (tp + fp) if tp + fp > 0 else 0
	rec  = tp / (tp + fn) if tp + fn > 0 else 0
	f1   = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0
	return {'precision': prec, 'recall': rec, 'f1': f1, 'tp': tp, 'fp': fp, 'fn': fn}

##########################
## SINGLE GENE EVAL     ##
##########################

def eval_gene(pred_gff, real_gff, pred_source=None, real_source=None):
	"""full evaluation of one gene: distribution distance + binary detection"""
	pred = get_introns(pred_gff, source=pred_source)
	real = get_introns(real_gff, source=real_source)
	if not pred or not real: return None

	dist = compare_introns(pred, real)
	binary = intron_f1(pred, real)

	return {
		'manhattan':    dist['manhattan'],
		'intersection': dist['intersection'],
		'dkl':          dist['dkl'],
		'chebyshev':    dist['chebyshev'],
		'n_pred':       dist['n_pred'],
		'n_real':       dist['n_real'],
		'precision':    binary['precision'],
		'recall':       binary['recall'],
		'f1':           binary['f1'],
		'details':      dist['details'],
	}

###########################
## BATCH EVALUATION      ##
###########################

def eval_batch(pred_dir, real_dir, pred_source=None, real_source=None):
	"""evaluate all genes in pred_dir against matching GFFs in real_dir"""
	results = []
	for fn in sorted(os.listdir(pred_dir)):
		if not fn.endswith('.gff3'): continue
		real_gff = os.path.join(real_dir, fn)
		if not os.path.exists(real_gff): continue
		pred_gff = os.path.join(pred_dir, fn)
		r = eval_gene(pred_gff, real_gff, pred_source, real_source)
		if r is not None:
			r['gene'] = fn.replace('.gff3', '')
			results.append(r)
	return results

def summary(results):
	"""aggregate statistics across genes"""
	if not results: return {}

	def stats(arr):
		a = np.array(arr)
		return {
			'mean':   float(a.mean()),
			'median': float(np.median(a)),
			'std':    float(a.std()),
			'min':    float(a.min()),
			'max':    float(a.max()),
		}

	return {
		'n_genes':      len(results),
		'manhattan':    stats([r['manhattan'] for r in results]),
		'intersection': stats([r['intersection'] for r in results]),
		'f1':           stats([r['f1'] for r in results]),
		'precision':    stats([r['precision'] for r in results]),
		'recall':       stats([r['recall'] for r in results]),
	}

##########
## CLI  ##
##########

if __name__ == '__main__':
	import argparse
	import json

	parser = argparse.ArgumentParser(description='evaluate intron predictions')
	sub = parser.add_subparsers(dest='cmd')

	# single gene comparison
	p1 = sub.add_parser('gene', help='compare two GFF files')
	p1.add_argument('pred', help='predicted GFF')
	p1.add_argument('real', help='real annotation GFF')
	p1.add_argument('--pred_source', default=None)
	p1.add_argument('--real_source', default=None)

	# batch comparison
	p2 = sub.add_parser('batch', help='evaluate all genes in directories')
	p2.add_argument('pred_dir', help='directory with predicted GFFs')
	p2.add_argument('real_dir', help='directory with real annotation GFFs')
	p2.add_argument('--pred_source', default=None)
	p2.add_argument('--real_source', default=None)

	arg = parser.parse_args()

	if arg.cmd == 'gene':
		r = eval_gene(arg.pred, arg.real, arg.pred_source, arg.real_source)
		if r is None:
			print('no introns found', file=sys.stderr)
			sys.exit(1)
		print(f'manhattan:    {r["manhattan"]:.4f}')
		print(f'intersection: {r["intersection"]:.4f}')
		print(f'dkl:          {r["dkl"]:.4f}')
		print(f'chebyshev:    {r["chebyshev"]:.4f}')
		print(f'precision:    {r["precision"]:.4f}')
		print(f'recall:       {r["recall"]:.4f}')
		print(f'f1:           {r["f1"]:.4f}')
		print(f'n_pred={r["n_pred"]}  n_real={r["n_real"]}')
		for (b, e), pp, rr in r['details']:
			print(f'  {b}\t{e}\t{pp:.6f}\t{rr:.6f}')

	elif arg.cmd == 'batch':
		results = eval_batch(arg.pred_dir, arg.real_dir,
			arg.pred_source, arg.real_source)
		s = summary(results)
		print(json.dumps(s, indent=2))

	else:
		parser.print_help()
