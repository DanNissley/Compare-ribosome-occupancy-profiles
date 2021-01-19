#!/usr/bin/env python

"""
### PURPOSE

This program implements, using the functions in the src/compare.py module, comparisons between
either randomly selected pairs of domains or between pairs of evolutionarily related domains.

"""
import os, sys
sys.path.append('/data/pegasus/nissley/trans_sched/update_v5/for_github/src')
import compare
from datetime import datetime

### usage statement
usage = './runner.py [input reads.tab file] [out_dir name] [k; 0 or 1] [random; True or False] [run output label] > stream.out'
if len(sys.argv) != 6:
	print (usage)
	sys.exit()

### paths to input files needed for this analysis

# path to SUPERFAM domain definitions
yeastSF             = 'inpfiles/Saccharomyces_cerevisiae_SUPERFAMILY_domains_SGD.txt'

# path the Structural Classification of Proteins; gives class, superfamily, family, etc.
scop                = 'inpfiles/dir.cla.scope.2.07-stable.txt'

# path to the consensus S288C coding DNA sequences
yeast_CDSs          = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta'

# path to the consesus S288C translated AA sequences
yeast_aa            = 'inpfiles/S288C_reference_genome_R64-2-1_20150113/orf_trans_all_R64-2-1_20150113.fasta'

# the input A-site profiles to be used in this analysis
yeast_Asite_reads   = sys.argv[1]

# make sure the A-site read file exists
if os.path.exists(yeast_Asite_reads) != True:
	print ('A-site reads file', yeast_Asite_reads, 'does not exist.')
	sys.exit()

# directory to which data will be written
out_dir             = sys.argv[2]

# if it does not already exist, create it
if os.path.exists(out_dir) == False:
	os.system('mkdir '+out_dir)

# integer, either 0 or 1; determines whether the first (0) or second (1) related profile
# will be used in comparisons to unrelated profiles. Only has minor influence on the results,
# and does not change conclusions at all
k                   = int(sys.argv[3])

# make sure k has an allowed value
if k not in [0, 1]:
	print ('A value of k='+str(k)+' is not allowed; choose from k = [0, 1].')
	sys.exit()

# iff True, perform random domain comparisons as a control
random              = compare.str2bool(sys.argv[4])

# integer with which to label the output from this independent run
run_label           = sys.argv[5]

# only those domain pairs with min_DNA_ID <= (domain pair DNA ID) <= max_DNA_ID will be analyzed
min_DNA_ID          = 0.30
max_DNA_ID          = 0.80

# directory from which input sequence alignments between related domains will be read
align_dir           = 'alignments/'

# iff True, run %MinMax analysis for comparison to the ribosome occupancy analysis
minmax              = False

### Start analysis

# get dictionary of Domain class objects
yeast_S288C_domains = compare.getDomains(yeastSF, scop, yeast_CDSs, yeast_aa, verbose=True)

# generate read profiles for each Domain class object (as possible)
compare.getReadProfiles(yeast_S288C_domains, yeast_Asite_reads, verbose=False)

# process raw read profiles in single domain translation rate profiles
compare.getSingleDomainProfiles(yeast_S288C_domains, verbose=False)

# what you do next depends on whether or not random analysis was requested
if random:

	# generate random pairs of domains
	pairs = compare.getRandomDomainPairs(yeast_S288C_domains, yeast_S288C_domains, 2000)

	# "align" domains; this is done naively based on the 5' end
	alignedDomains = compare.alignRandomReadProfiles(pairs, yeast_S288C_domains)

# if we want to run the actual comparisons between related domain profiles
else:

	# get pairs of domains within the same SUPERFAM family (using each pair once but using individual domains multiple times)
	pairs = compare.getHomologs(yeast_S288C_domains, yeast_S288C_domains)

	# perform AA and DNA alignments - commented out so we can just use alignments in alignments/ and
	# skip that step for all future analyses; saves time and disk space!
	#compare.alignDomains(yeast_S288C_domains, yeast_S288C_domains, pairs, out_dir)

	# select pairs of domains with acceptable DNA sequence identity
	goodPairs = compare.findAcceptableDomainPairs(pairs, align_dir, min=min_DNA_ID, max=max_DNA_ID, verbose=False)

	# align translation rate profiles based on AA alignments, throwing out pairs with unacceptable gaps
	alignedDomains = compare.getAlignedReadProfiles(goodPairs, yeast_S288C_domains, align_dir, maxGapSize=5, maxNumGaps=10,  verbose=False)

# perform ribosome occupancy profile comparisons
compare.runComparison(alignedDomains, yeast_S288C_domains, out_dir, run_label, k)
