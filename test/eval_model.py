#!/usr/bin/env python3

"""
Splice model quality evaluation.

Evaluates whether the built splice model (PWM, Markov, lengths) actually
captures the biology. The upstream isoforms project's ian/cmp-pwm.py
compares different PWM building methods. This module measures:

  1. PWM information content — do the donor/acceptor PWMs show strong consensus?
  2. Splice signal detection — do PWMs score real sites higher than random GT/AG?
  3. Markov model discrimination — do exon/intron models distinguish exonic vs intronic seq?
  4. Length model fit — does the length distribution match observed data?
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))
import isoform

########################
## PWM QUALITY        ##
########################

def pwm_info_content(pwm):
	"""information content per position of a PWM (bits, max 2 for DNA)"""
	ic = []
	for pos in pwm:
		h = 0
		for nt in pos:
			p = pos[nt]
			if p > 0: h -= p * math.log2(p)
		ic.append(2.0 - h)
	return ic

def pwm_total_info(pwm):
	"""total information content of PWM (sum across positions)"""
	return sum(pwm_info_content(pwm))

def eval_donor_pwm(don_pwm):
	"""evaluate donor PWM quality — GT consensus should be strong"""
	ic = pwm_info_content(don_pwm)
	return {
		'total_ic':      sum(ic),
		'per_position':  ic,
		'mean_ic':       np.mean(ic),
		'max_ic':        max(ic),
	}

def eval_acceptor_pwm(acc_pwm):
	"""evaluate acceptor PWM quality — AG consensus should be strong"""
	ic = pwm_info_content(acc_pwm)
	return {
		'total_ic':      sum(ic),
		'per_position':  ic,
		'mean_ic':       np.mean(ic),
		'max_ic':        max(ic),
	}

################################
## SPLICE SIGNAL DETECTION    ##
################################

def score_real_vs_random(pwm, seq, real_positions, flank=99, minex=25, site_type='don'):
	"""compare PWM scores at real splice sites vs random GT/AG sites

	Tests whether the model can discriminate real sites from noise.
	Returns score distributions for real and decoy sites.
	"""
	all_dons, all_accs = isoform.gtag_sites(seq, flank, minex)
	real = set(real_positions)

	if site_type == 'don':
		all_sites = all_dons
	else:
		all_sites = all_accs

	real_scores = []
	decoy_scores = []

	for pos in all_sites:
		# check bounds
		if site_type == 'don':
			start = pos
		else:
			start = pos - len(pwm) + 1
		if start < 0 or start + len(pwm) > len(seq): continue

		s = isoform.score_pwm(pwm, seq, start)
		if pos in real: real_scores.append(s)
		else:           decoy_scores.append(s)

	if not real_scores or not decoy_scores:
		return {'separable': False}

	real_arr = np.array(real_scores)
	decoy_arr = np.array(decoy_scores)

	return {
		'separable':   True,
		'real_mean':   float(real_arr.mean()),
		'real_std':    float(real_arr.std()),
		'decoy_mean':  float(decoy_arr.mean()),
		'decoy_std':   float(decoy_arr.std()),
		'separation':  float(real_arr.mean() - decoy_arr.mean()),
		'n_real':      len(real_scores),
		'n_decoy':     len(decoy_scores),
	}

##################################
## MARKOV MODEL DISCRIMINATION  ##
##################################

def eval_markov_discrimination(exs_model, ins_model, exon_seqs, intron_seqs,
		don_ctx=9, acc_ctx=9):
	"""test if exon/intron Markov models can discriminate sequence types

	Score exon sequences with both models, same for intron sequences.
	If models are good, exon model should score exons higher than introns.
	"""
	exon_exs = []   # exon seqs scored by exon model
	exon_ins = []   # exon seqs scored by intron model
	intron_exs = [] # intron seqs scored by exon model
	intron_ins = [] # intron seqs scored by intron model

	k = exs_model['k']

	for seq in exon_seqs:
		if len(seq) < k + 1: continue
		try:
			exon_exs.append(isoform.score_markov(exs_model, seq, 0, len(seq) - 1))
			exon_ins.append(isoform.score_markov(ins_model, seq, 0, len(seq) - 1))
		except KeyError:
			continue

	for seq in intron_seqs:
		interior_beg = don_ctx
		interior_end = len(seq) - acc_ctx - 1
		if interior_end - interior_beg < k + 1: continue
		try:
			intron_exs.append(isoform.score_markov(exs_model, seq,
				interior_beg, interior_end))
			intron_ins.append(isoform.score_markov(ins_model, seq,
				interior_beg, interior_end))
		except KeyError:
			continue

	if not exon_exs or not intron_exs:
		return {'testable': False}

	# per-base normalization
	ee = np.array(exon_exs)
	ei = np.array(exon_ins)
	ie = np.array(intron_exs)
	ii = np.array(intron_ins)

	return {
		'testable':             True,
		'exon_by_exon_model':   float(ee.mean()),
		'exon_by_intron_model': float(ei.mean()),
		'intron_by_exon_model': float(ie.mean()),
		'intron_by_intron_model': float(ii.mean()),
		'exon_discrimination':  float(ee.mean() - ei.mean()),
		'intron_discrimination': float(ii.mean() - ie.mean()),
		'n_exons':  len(exon_exs),
		'n_introns': len(intron_exs),
	}

#############################
## LENGTH MODEL FIT        ##
#############################

def eval_length_fit(length_model, observed_lengths):
	"""evaluate how well a length model fits observed data

	Computes mean log-likelihood of observed lengths under the model.
	Higher = better fit.
	"""
	scores = []
	for l in observed_lengths:
		scores.append(isoform.score_len(length_model, l))

	arr = np.array(scores)
	return {
		'mean_score':   float(arr.mean()),
		'std_score':    float(arr.std()),
		'min_score':    float(arr.min()),
		'max_score':    float(arr.max()),
		'n_observations': len(scores),
	}

###############################
## FULL MODEL EVALUATION     ##
###############################

def eval_splice_model(model_path):
	"""load and evaluate a splice model JSON"""
	model = isoform.read_splicemodel(model_path)

	don_eval = eval_donor_pwm(model['don'])
	acc_eval = eval_acceptor_pwm(model['acc'])

	return {
		'name':        model.get('name', ''),
		'donor_pwm':   don_eval,
		'acceptor_pwm': acc_eval,
		'don_length':  len(model['don']),
		'acc_length':  len(model['acc']),
		'exs_kmers':   len(model['exs']['mm']),
		'ins_kmers':   len(model['ins']['mm']),
		'exl_size':    model['exl']['size'],
		'inl_size':    model['inl']['size'],
		'inf':         model['inf'],
	}

##########
## CLI  ##
##########

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='evaluate splice model quality')
	parser.add_argument('model', help='splice model JSON file')
	arg = parser.parse_args()

	r = eval_splice_model(arg.model)
	print(f'model: {r["name"]}')
	print(f'donor PWM:    {r["don_length"]}bp, IC={r["donor_pwm"]["total_ic"]:.2f} bits')
	print(f'acceptor PWM: {r["acc_length"]}bp, IC={r["acceptor_pwm"]["total_ic"]:.2f} bits')
	print(f'exon Markov:  {r["exs_kmers"]} kmers')
	print(f'intron Markov: {r["ins_kmers"]} kmers')
	print(f'exon lengths: {r["exl_size"]} bins')
	print(f'intron lengths: {r["inl_size"]} bins')
	print(f'intron freq:  {r["inf"]:.3f}')
	print()
	print('donor IC per position:')
	for i, ic in enumerate(r['donor_pwm']['per_position']):
		bar = '#' * int(ic * 20)
		print(f'  {i}: {ic:.3f} {bar}')
	print()
	print('acceptor IC per position:')
	for i, ic in enumerate(r['acceptor_pwm']['per_position']):
		bar = '#' * int(ic * 20)
		print(f'  {i}: {ic:.3f} {bar}')
