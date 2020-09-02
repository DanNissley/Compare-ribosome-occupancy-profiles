#!/usr/bin/env python

import os, sys
sys.path.append('/data/icarus/nissley/trans_sched/fullPipeline/update_v4/trans_sched/src')
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

# perform AA and DNA alignments
compare.alignDomains(yeast_S288C_domains, yeast_S288C_domains, pairs, outdir)
