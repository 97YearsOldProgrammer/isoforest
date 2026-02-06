#!/usr/bin/env python3

import argparse
import json
import math
import os
import sys

import numpy as np

from grimoire.sequence i	mport DNA
from grimoire.feature 	import Feature
from grimoire.genome 	import Reader

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))
import isoform
import randomf

################
## CONSTANTS ##
################

MIN_INTRON = 35
MIN_EXON   = 25
DON_CTX    = 9		# donor PWM context (3 upstream + GT + 4 downstream)
ACC_CTX    = 9		# acceptor PWM context (4 upstream + AG + 3 downstream)
MM_ORDER   = 1		# Markov model order

##############
## FILTERS ##
##############

def isolated(genes, log):
	"""marks BOTH overlapping partners, early break on sorted genes"""
	bad = set()
	for i in range(len(genes)):
		for j in range(i+1, len(genes)):
			if genes[j].beg > genes[i].end: break
			if genes[i].overlap(genes[j]):
				bad.add(i)
				bad.add(j)
	keep = []
	for i in range(len(genes)):
		if i in bad: log['gene-overlap'] += 1
		else:        keep.append(genes[i])
	return keep

def protein_coding(genes, log):
	keep = []
	for gene in genes:
		if gene.is_coding(): keep.append(gene)
		else:                log['non-coding'] += 1
	return keep

def has_introns(genes, log):
	keep = []
	for gene in genes:
		found = False
		for tx in gene.transcripts():
			if len(tx.introns) > 0:
				found = True
				break
		if found: keep.append(gene)
		else:     log['not-spliced'] += 1
	return keep

def canonical_introns(genes, log):
	keep = []
	for gene in genes:
		if gene.issues: log['non-canonical'] += 1
		else:           keep.append(gene)
	return keep

def not_too_long(genes, maxlen, log):
	keep = []
	for gene in genes:
		if gene.length > maxlen: log['too-long'] += 1
		else:                    keep.append(gene)
	return keep

def in_bounds(genes, padding, seqlen, log):
	keep = []
	for gene in genes:
		beg = gene.beg - padding - 1
		end = gene.end + padding
		if beg < 0 or end >= seqlen: log['out-of-bounds'] += 1
		else:                        keep.append(gene)
	return keep

#################
## EXTRACTION ##
#################

def count_introns(gene):
	best = 0
	for tx in gene.transcripts():
		n = len(tx.introns)
		if n > best: best = n
	return best

def build_smallgene(chrom, gene, padding, sources, maxiso):
	"""isoform filter + extraction in one pass. returns (dna, n_gt, n_ag) or None"""
	beg = gene.beg - padding - 1
	end = gene.end + padding
	seq = chrom.seq[beg:end+1]

	# isoform count on plus-strand-oriented sequence
	oriented = isoform.anti(seq) if gene.strand == '-' else seq
	dons, accs = isoform.gtag_sites(oriented, padding, MIN_EXON)
	ic = randomf.countiso(dons, accs, MIN_INTRON, MIN_EXON, limit=maxiso)
	if ic >= maxiso: return None

	# reuse same seq for DNA object
	desc = f'{chrom.name}:{gene.beg}-{gene.end} {gene.strand} {gene.id} ISO={ic}'
	dna = DNA(seq=seq, desc=desc)

	stuff = chrom.ftable.fetch(beg, end)
	for f in stuff:
		if f.strand != gene.strand: continue
		if f.beg < beg + padding:   continue
		if f.end > end - padding:   continue
		if f.source not in sources: continue
		dna.ftable.add_feature(Feature(dna, f.beg - beg, f.end - beg,
			f.strand, f.type, phase=f.phase,
			score=f.score, source=f.source, id=f.id, pid=f.pid))

	if gene.strand == '-': dna.revcomp()
	return dna, len(dons), len(accs)

def write_gene(outdir, name, dna):
	with open(f'{outdir}/{name}.fa', 'w') as fp:
		fp.write(dna.fasta())
	with open(f'{outdir}/{name}.gff3', 'w') as fp:
		for f in dna.ftable.features:
			fp.write(f.gff())
			fp.write('\n')

##########################
## SPLICE MODEL BUILDER ##
##########################

def collect_splice_data(gene):
	"""collect unique donor/acceptor/exon/intron sequences from gene"""
	intron_seen = set()
	exon_seen = set()
	dons, accs = [], []
	exon_seqs, intron_seqs = [], []
	exon_lens, intron_lens = [], []
	for tx in gene.transcripts():
		for intron in tx.introns:
			sig = (intron.beg, intron.end)
			if sig in intron_seen: continue
			intron_seen.add(sig)
			iseq = intron.seq_str()
			if len(iseq) >= max(DON_CTX, ACC_CTX):
				dons.append(iseq[:DON_CTX])
				accs.append(iseq[-ACC_CTX:])
			intron_seqs.append(iseq)
			intron_lens.append(intron.length)
		for exon in tx.exons:
			sig = (exon.beg, exon.end)
			if sig in exon_seen: continue
			exon_seen.add(sig)
			exon_seqs.append(exon.seq_str())
			exon_lens.append(exon.length)
	return dons, accs, exon_seqs, intron_seqs, exon_lens, intron_lens

def length_hist(lengths):
	"""normalized length distribution via numpy"""
	counts = np.bincount(lengths)
	return (counts / counts.sum()).tolist()

def pwm_to_logodds(pwm):
	"""convert probability PWM to log-odds [{nt: score}, ...]"""
	return [{nt: isoform.prob2score(pos[nt]) for nt in pos} for pos in pwm]

def markov_to_logodds(mm, order):
	"""convert probability Markov {ctx: {nt: p}} to log-odds {k, mm: {kmer: score}}"""
	scored = {}
	for ctx in mm:
		for nt in mm[ctx]:
			scored[ctx + nt] = isoform.prob2score(mm[ctx][nt])
	return {'k': order + 1, 'mm': scored}

def len_to_logodds(hist):
	"""convert probability histogram to log-odds length model {tail, size, val}"""
	size = len(hist)
	tail_p = hist[-1] if hist[-1] > 0 else 1e-10
	tail = isoform.find_tail(tail_p, size)
	expect = 1 / size
	val = []
	for p in hist:
		if p == 0: val.append(-100)
		else:      val.append(math.log2(p / expect))
	return {'tail': tail, 'size': size, 'val': val}

def write_model(outdir, gid, dons, accs, exon_seqs, intron_seqs,
		exon_lens, intron_lens, total_gt, total_ag):
	"""build splice model: raw component files + assembled JSON"""

	# raw probability models
	don_pwm = isoform.create_pwm(dons)
	acc_pwm = isoform.create_pwm(accs)
	exs_mm  = isoform.create_markov(exon_seqs, MM_ORDER, 0, 0)
	ins_mm  = isoform.create_markov(intron_seqs, MM_ORDER, DON_CTX, ACC_CTX)
	exl     = length_hist(exon_lens)
	inl     = length_hist(intron_lens)

	# write inspectable component files
	isoform.write_pwm(f'{outdir}/{gid}.don.pwm', don_pwm)
	isoform.write_pwm(f'{outdir}/{gid}.acc.pwm', acc_pwm)
	isoform.write_markov(f'{outdir}/{gid}.exon.mm', exs_mm)
	isoform.write_markov(f'{outdir}/{gid}.intron.mm', ins_mm)
	isoform.write_len(f'{outdir}/{gid}.exon.len', exl)
	isoform.write_len(f'{outdir}/{gid}.intron.len', inl)

	# intron frequency: log-odds of a GT being a real donor
	p_real = len(dons) / total_gt if total_gt > 0 else 0.01
	inf = math.log2(p_real / (1 - p_real)) if 0 < p_real < 1 else 0

	# assemble complete splice model JSON (all log-odds)
	model = {
		'name': gid,
		'don':  pwm_to_logodds(don_pwm),
		'acc':  pwm_to_logodds(acc_pwm),
		'exs':  markov_to_logodds(exs_mm, MM_ORDER),
		'ins':  markov_to_logodds(ins_mm, MM_ORDER),
		'exl':  len_to_logodds(exl),
		'inl':  len_to_logodds(inl),
		'inf':  inf,
	}
	with open(f'{outdir}/{gid}.splicemodel.json', 'w') as fp:
		json.dump(model, fp, indent=2)

	print(f'model: {len(dons)} donors, {len(accs)} acceptors', file=sys.stderr)
	print(f'model: {len(exon_seqs)} exons, {len(intron_seqs)} introns',
		file=sys.stderr)
	print(f'model: inf={inf:.3f} (P(real|GT)={p_real:.4f})', file=sys.stderr)
	print(f'model: wrote {gid}.splicemodel.json', file=sys.stderr)

##########
## MAIN ##
##########

parser = argparse.ArgumentParser(description='build smallgenes dataset')
parser.add_argument('fasta', type=str, metavar='<fasta>',
	help='path to genome fasta file')
parser.add_argument('gff', type=str, metavar='<gff>',
	help='path to genome gff3 file')
parser.add_argument('out', type=str, metavar='<out>',
	help='output directory')
parser.add_argument('--genome_id', required=False, type=str, default='dm',
	metavar='<str>', help='prefix for output files [%(default)s]')
parser.add_argument('--padding', required=False, type=int, default=99,
	metavar='<int>', help='flanking region on each end [%(default)d]')
parser.add_argument('--maxgene', required=False, type=int, default=1000,
	metavar='<int>', help='maximum gene length [%(default)d]')
parser.add_argument('--max_isoforms', required=False, type=float, default=10,
	metavar='<float>', help='maximum isoforms in millions [%(default)f]')
parser.add_argument('--sources', required=False, type=str, nargs='+',
	default=['Gnomon', 'BestRefSeq%2CGnomon', 'BestRefSeq'],
	metavar='<str>', help='annotation sources to keep [%(default)s]')
parser.add_argument('--build_model', action='store_true',
	help='build D. mel splice model (PWM, Markov, lengths, JSON)')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

if not os.path.exists(arg.out): os.mkdir(arg.out)

genome  = Reader(gff=arg.gff, fasta=arg.fasta)
maxiso  = int(arg.max_isoforms * 1e6)
sources = set(arg.sources)
idx     = 0

log = {
	'smallgenes':    0,
	'gene-overlap':  0,
	'non-coding':    0,
	'not-spliced':   0,
	'non-canonical': 0,
	'too-long':      0,
	'out-of-bounds': 0,
	'too-complex':   0,
}

# splice model accumulators
all_dons, all_accs             = [], []
all_exon_seqs, all_intron_seqs = [], []
all_exon_lens, all_intron_lens = [], []
total_gt, total_ag             = 0, 0

for chrom in genome:
	genes = chrom.ftable.build_genes()

	# cheap filters first
	genes = isolated(genes, log)
	genes = protein_coding(genes, log)
	genes = has_introns(genes, log)
	genes = canonical_introns(genes, log)
	genes = not_too_long(genes, arg.maxgene, log)
	genes = in_bounds(genes, arg.padding, len(chrom.seq), log)

	# expensive filter + extraction in single pass
	for gene in genes:
		result = build_smallgene(chrom, gene, arg.padding, sources, maxiso)
		if result is None:
			if arg.verbose: print('>', end='', file=sys.stderr, flush=True)
			log['too-complex'] += 1
			continue
		dna, n_gt, n_ag = result
		if arg.verbose: print('.', end='', file=sys.stderr, flush=True)

		idx += 1
		name = f'{arg.genome_id}.{count_introns(gene)}_{idx}'
		dna.name = name
		write_gene(arg.out, name, dna)
		log['smallgenes'] += 1

		if arg.build_model:
			d, a, es, is_, el, il = collect_splice_data(gene)
			all_dons.extend(d)
			all_accs.extend(a)
			all_exon_seqs.extend(es)
			all_intron_seqs.extend(is_)
			all_exon_lens.extend(el)
			all_intron_lens.extend(il)
			total_gt += n_gt
			total_ag += n_ag

	if arg.verbose:
		print(f'\n{chrom.name} {log}', file=sys.stderr, flush=True)

# assemble and write splice model
if arg.build_model and all_dons:
	write_model(arg.out, arg.genome_id, all_dons, all_accs,
		all_exon_seqs, all_intron_seqs, all_exon_lens, all_intron_lens,
		total_gt, total_ag)

print(log, file=sys.stderr)
