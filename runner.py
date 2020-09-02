#!/usr/bin/env python

import os, sys
sys.path.append('/data/icarus/nissley/trans_sched/fullPipeline/update_v4/trans_sched/src/')
import compare
from datetime import datetime

# usage statement
usage = './runner.py [input reads.tab file] [outdir name] [k; 0 or 1] [random; True or False] [run/trajectory output label] > stream.out'
if len(sys.argv) != 6:
	print (usage)
	sys.exit()

# paths to input files
yeastSF             = 'inpfiles/Saccharomyces_cerevisiae_SUPERFAMILY_domains_SGD.txt'
scop                = 'inpfiles/dir.cla.scope.2.07-stable.txt'
yeast_CDSs          = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta'
yeast_aa            = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_trans_all_R64-2-1_20150113.fasta'
yeast_Asite_reads   = sys.argv[1]
outdir              = sys.argv[2]
k                   = int(sys.argv[3])
random              = compare.str2bool(sys.argv[4])
traj                = sys.argv[5]

# make output directory; if it already exists, delete it
#compare.mkdir(outdir)

# make sure the A-site read file exists
if os.path.exists(yeast_Asite_reads) != True:
	print ('A-site reads file', yeast_Asite_reads, 'does not exist.')
	sys.exit()

# make sure k was requested to have one of its allowed values
if k not in [0, 1]:
	print ('A value of k='+str(k)+' is not allowed; choose from k = [0, 1].')
	sys.exit()

# if random control analysis was requested
if random:

	# get dictionary of Domain class objects
	yeast_S288C_domains = compare.getDomains(yeastSF, scop, yeast_CDSs, yeast_aa, verbose=True)

	# generate read profiles for each Domain class object (as possible)
	compare.getReadProfiles(yeast_S288C_domains, yeast_Asite_reads, verbose=False)

	# process raw read profiles in single domain translation rate profiles
	compare.getSingleDomainProfiles(yeast_S288C_domains, verbose=False)

	# generate random pairs of domains
	pairs = compare.getRandomDomainPairs(yeast_S288C_domains, yeast_S288C_domains, 2000)

	# "align" domains; this is done naively based on the 5' end
	alignedDomains = compare.alignRandomReadProfiles(pairs, yeast_S288C_domains)

	# perform translation rate profile comparisons
	compare.runComparison(alignedDomains, yeast_S288C_domains, outdir, traj, k)

# if we want to run the actual comparisons between related domain profiles
else:

	# get dictionary of Domain class objects
	yeast_S288C_domains = compare.getDomains(yeastSF, scop, yeast_CDSs, yeast_aa, verbose=True)

	# generate read profiles for each Domain class object (as possible)
	compare.getReadProfiles(yeast_S288C_domains, yeast_Asite_reads, normalize=False, verbose=False)

	# process raw read profiles in single domain translation rate profiles
	compare.getSingleDomainProfiles(yeast_S288C_domains, verbose=False)

	# get pairs of domains within the same SUPERFAM family (using each pair once but using individual domains multiple times)
	pairs = compare.getHomologs(yeast_S288C_domains, yeast_S288C_domains)

	# perform AA and DNA alignments
	#compare.alignDomains(yeast_S288C_domains, yeast_S288C_domains, pairs, outdir)

	# select pairs of domains with acceptable DNA sequence identity
	goodPairs = compare.findAcceptableDomainPairs(pairs, outdir, min=0.30, max=0.80, verbose=False)

	# align translation rate profiles based on AA alignments, throwing out pairs with unacceptable gaps
	alignedDomains = compare.getAlignedReadProfiles(goodPairs, yeast_S288C_domains, outdir, maxGapSize=5, maxNumGaps=10,  verbose=False)

	# perform translation rate profile comparisons
	compare.runComparison(alignedDomains, yeast_S288C_domains, outdir, traj, k)
