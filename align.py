#!/usr/bin/env python

"""
### Purpose

This program generates alignments based on the input A-site reads profile pairs that have
sufficient read coverage. By running the program once to generate alignments you can speed
up future analyses/running multiple sets of comparisons by reading in these alignments rather
than generating them each time.

"""
import os, sys
sys.path.append('/data/pegasus/nissley/trans_sched/update_v5/for_github/src')
import compare
from datetime import datetime

# usage statement
usage = './align.py [input reads.tab file] [outdir name]'
if len(sys.argv) != 3:
	print (usage)
	sys.exit()

# paths to input files
yeastSF             = 'inpfiles/Saccharomyces_cerevisiae_SUPERFAMILY_domains_SGD.txt'
scop                = 'inpfiles/dir.cla.scope.2.07-stable.txt'
yeast_CDSs          = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta'
yeast_aa            = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_trans_all_R64-2-1_20150113.fasta'
yeast_Asite_reads   = sys.argv[1]
outdir              = sys.argv[2]
min_DNA_ID          = 0.30
max_DNA_ID          = 0.80

# make output directory; if it already exists, delete it
compare.mkdir(outdir)

# make sure the A-site read file exists
if os.path.exists(yeast_Asite_reads) != True:
	print ('A-site reads file', yeast_Asite_reads, 'does not exist.')
	sys.exit()

# get dictionary of Domain class objects
yeast_S288C_domains = compare.getDomains(yeastSF, scop, yeast_CDSs, yeast_aa, verbose=True)

# generate read profiles for each Domain class object (as possible)
compare.getReadProfiles(yeast_S288C_domains, yeast_Asite_reads, normalize=False, verbose=False)

# process raw read profiles in single domain translation rate profiles
compare.getSingleDomainProfiles(yeast_S288C_domains, verbose=False)

# get pairs of domains within the same SUPERFAM family (using each pair once but using individual domains multiple times)
pairs = compare.getHomologs(yeast_S288C_domains, yeast_S288C_domains, verbose=False)

# perform AA and DNA alignments; add multi=False argument to run in series
compare.alignDomains(yeast_S288C_domains, yeast_S288C_domains, pairs, outdir)

# use this function with verbose = True to write the DNA sequence identity file to alignments/dna_alignments/dna_sequence_ID.dat
goodPairs = compare.findAcceptableDomainPairs(pairs, outdir, min=min_DNA_ID, max=max_DNA_ID, verbose=True)
