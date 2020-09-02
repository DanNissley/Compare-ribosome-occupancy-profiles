#!/usr/bin/env python

"""
Daniel A. Nissley, Ph.D.
University of Oxford
Department of Statistics
Deane Lab
Oxford Protein Informatics Group

Original version: 30-July-2019
Current version : 28-July-2020

Purpose: To compare ribosome occupancy profiles between related domains in yeast

Notes:

(1)  Look into the inconsistencies found in some SUPERFAM entries (i.e., finding '-' in a superfamily identifier, etc.)
(2)  Why does superfamily 103481 not have a SCOP class? - DONE (see not in getDomains function)
(3)  Why do some genes have an incorrect number of codons in the A-site read file? See getReadProfiles function
     RESOLVED - 28 July 2020. This was caused by an erroneous assumption that checking "if orf in key" would result in unique assignment of
     reads to the correct domain. However, because some domains have names with "-A" and "-B" subscripts, this leads to problems because they
     will not be uniquely identified. I have fixed this by changing the check in getDNAseq, getAAseq, and getReadProfiles to "if orf == key.split(';')[0]",
     which guarantees strict matching of ORF names between different input files. - DONE
(4)  Add back in controls that check number of positions with no reads, etc. before running alignments - DONE
(5)  Cannot use ORF ID alone as the key for the Domain object dictionaries as many genes have multiple domains; add the amino
     range as well so we can specify - DONE
(6)  Make sure we check for overlapping domain defintions when domains are in the same protein - DONE
(7)  Look into mismatches between domain definitions and sequence information. This is likely a problem with
     crossreferencing data between SGD and SUPERFAM. For example, the DNA sequence for YBL057C contains 209 codons
     in S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta, but SUPERFAM has a domain that
     spans YBL057C;93-214. There are <20 cases like this but it may indicate a subtle shift that influences other domains/genes
(8)  Another problem likely stemming from mismatches between SUPERFAM and SGD - SUPERFAM has the domain YDR098C;662-812
     but in S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta the protein only has 250 amino acids...
     These mismatches and those as in (7) are currently filtered out and not used for downstream analyses
(9)  Why are there so many domains in the SUPERFAM file that seem to be the same exact thing? For example, grep 662-812 and notice
     that despite having different ORF IDs they have the exact same residue numbers and closest PDB (1ASU).
(10) Found and corrected a minor big in previous version that led to the 61st amino acid being omitted from .fa files used for alignments; should not have any influence on results
(11) Currently doing DNA alignments using an extra codon on the 5' end; I could not decide whether I wanted to use the DNA sequence that determines
     the translation schedule for the domain or the DNA sequence that encodes for the amino acids in the domain, so I took the codons that do both (should not seriously influence results, just 3 nucleotides
     out of hundreds to thousands for a typical domain definition)
"""

### IMPORT MODULES

import os
import sys
import numpy as np
import multiprocessing as mp
from datetime import datetime
import scipy.integrate
import scipy.interpolate
import scipy.stats
import subprocess
from datetime import datetime
import random

### CLASSES

class Domain:

	# Domain objects contain information about a particular protein domain as defined by SUPERFAMILY

	# family             : str, the SUPERFAM family (recent common ancestor) to which this domain belongs
	# superfamily        : str, the SUPERFAM superfamily (distant common ancestor) to which this domain belongs
	# orfID              : str, the open reading frame ID for this domain (e.g., YAL002W, YBL080C, etc.)
	# resStart           : int, the first residue included in this domain in primary sequence
	# resEnd             : int, the last residue included in this domain in the primary sequence
	# SCOPclass          : str, the SCOP classification of this gene
	# readProfile        : np.array ((), dtype=np.float64), A-site reads per domain position
	# DNAseq             : str, the DNA sequence for the domain
	# AAseq              : str, the AA sequence for the domain

	# instantiation of new class objects based on domain parameters
	def __init__ (self, family, superfamily, orfID, resStart, resEnd, SCOPclass):
		self.family = family
		self.superfamily = superfamily
		self.orfID = orfID
		self.resStart = resStart
		self.resEnd = resEnd
		self.SCOPclass = SCOPclass
		self.readProfile = None
		self.DNAseq = None
		self.AAseq = None
		self.ORFlength = None
		self.singleDomainProfile = None

	# function to add A-site read profile for a given domain
	def addReadProfile (self, readProfile):
		self.readProfile = readProfile

	# function to add the DNA sequence for a given domain
	def addDNA (self, DNAseq):
		self.DNAseq = DNAseq

	# function to add amino acid sequence to facilitate quality controls
	def addAAseq (self, AAseq):
		self.AAseq = AAseq

	# function to add the ORF length - necessary to produce the trimmed translation rate profile
	def addORFlength(self, ORFlength):
		self.ORFlength = ORFlength

	# function to add single-domain profiles (see getSingleDomainProfiles for details and explanation)
	def addSingleDomainProfile(self, singleDomainProfile):
		self.singleDomainProfile = singleDomainProfile

### FUNCTIONS

# function to generate Domain class objects
def getDomains (f_superfam, f_scop, f_DNA_fasta, f_AA_fasta, min=100, verbose=False, ecoli=False):

	# This is the key initial function to call on a dataset - it brings together information on domain definitions,
	# structural classes, SUPERFAMILY annotations, and DNA sequences. It does not, however, parse information related
	# to the actual ribosome profiling data - this is handled by get readProfile.

	# f_superfam: str, path to SUPERFAM database file, e.g. Saccharomyces_cerevisiae_SUPERFAMILY_domains_SGD.txt
	# f_scop: str, path to SCOP database file, e.g., dir.cla.scope.2.07-stable.txt
	# f_DNA_fasta: str, path to file containing DNA sequences for this organism
	# f_AA_fasta: str, path to file containinf amino acid sequences for this organism
	# min: (optional, default=100) int, the minimum size of a domain in residues for us to include it the list of domains for future analysis
	# verbose: (optional, default=False) Boolean, determines whether or not additional non-essential information is printed
	# ecoli: (optional, default=False) Boolean, determines whether or not we get orfIDs based on the format in the E. coli SUPERFAM file, which differs from the yeast format

	# returns: a list of Domain class objects, one for each domain defined in the input SUPERFAM file

	domains = {}
	superfamDB = readFile(f_superfam)
	scopDB = readFile(f_scop)
	mapSUPERFAMILYtoSCOP = {} # dictionary; key -> superfamily, value -> SCOP class
	SCOPclass = ''
	superfamily = ''
	family = ''
	orfID = ''
	aaRange = ''
	N = 0

	if verbose:
		print ('Parsing domain defintions in', f_superfam)

	# populate mapSUPERFAMILYtoSCOP with superfamily/SCOP class pairs
	for m in scopDB:
		if m[0] == '#':
			continue
		m = m.split()
		SCOPclass = m[3].split('.')[0]
		superfamily = m[5].split(',')[2]
		if SCOPclass not in ['a','b','c','d','e','f','g','h','i','j','k','l'] and verbose:
			print ('SCOP class is not in the standard list of a, b, c, d, e, f, g, h, i, j, k, l for the line:')
			print (m)
			sys.exit()
		if 'sf' not in superfamily and verbose:
			print ('SUPERFAMILY incorrectly selected from SCOP database.')
			print (m)
			sys.exit()
		superfamily = superfamily.strip('sf=')
		mapSUPERFAMILYtoSCOP[superfamily] = SCOPclass

	# add in the SCOP class for superfamily 103481 manually; it appears as 'f' on the webserver (see https://scop.berkeley.edu/sunid=103481&ver=1.67) but not in the downloaded file
	mapSUPERFAMILYtoSCOP['103481'] = 'f'

	#for sf in mapSUPERFAMILYtoSCOP:
	#	print (sf, mapSUPERFAMILYtoSCOP[sf])

	# create a Domain class object for all domains listed in superfamDB
	superfamily = ''
	for k in superfamDB:
		k = k.strip('\n').split('\t')
		if k != [''] and k[0][0] not in ['#', 'G']:
			if ecoli:
				orfID = k[1].split('|')[3] # gi|16131605|ref|NP_418193.1|
			else:
				orfID = k[1]
				if orfID[0] in ['R', 'Q']: # ignore plasmids and mitochondrial genes
					continue
			aaRange = k[3]
			superfamily = k[5]
			family = k[8]
			if '-' in family: # some genes have weird families; ignore them
				continue
			if '-' in superfamily: # some genes have weird superfamilies; ignore them
				continue
			if ',' in aaRange: # also ignore non-contiguous domain definitions
				continue
			if superfamily not in mapSUPERFAMILYtoSCOP:
				if verbose:
					print ('SUPERFAMILY', superfamily, 'does not have a SCOP classification in', f_scop, 'and will be omitted.')
				continue
			# make sure the domain meets the size criterion set by the user
			if (int(aaRange.split('-')[1]) - int(aaRange.split('-')[0]) + 1) < min:
				continue
			# if all quality controls have been passed create a Domain object for this entry and append it to the list of domains
			domains[orfID+';'+aaRange.split('-')[0]+'-'+aaRange.split('-')[1]] = Domain(family, superfamily, orfID, int(aaRange.split('-')[0]), int(aaRange.split('-')[1]), mapSUPERFAMILYtoSCOP[superfamily])
			if verbose:
				print ((orfID+';'+aaRange.split('-')[0]+':'+aaRange.split('-')[1]).ljust(20)+'\t'+mapSUPERFAMILYtoSCOP[superfamily]+'\t'+superfamily+'\t'+family)
			N += 1

	# add DNA sequence for each Domain object
	getDNAseq(domains, f_DNA_fasta)
	# add amino acid sequence for each Domain object
	getAAseq(domains, f_AA_fasta, ecoli=ecoli)

	# quality control; make sure the amino acid and DNA sequences correspond exactly to one another
	if not qcDNAandAAseqs(domains):
		print ('Quality control checks for DNA and AA sequences failed - cannot proceed with incorrect sequence information!')
		sys.exit()
	if verbose:
		print ('Correctly parsed', N, 'domains in', f_superfam, '\n')

	# return the list of Domain class objects
	return domains

# function to get DNA FASTA sequences for each domain
def getDNAseq (domains, f_DNA_fasta):

	# domains: dictionary; keys are 'orfID;resStart-resEnd', values are the corresponding Domain objects
	# f_DNA_fasta: str; path to DNA FASTA file

	# returns: nothing, but the domains dictionary contents are modified to include an mRNA sequence as needed

	# N.B., this function takes an "extra" DNA codon to facilitate DNA/AA comparisons and translation rate profile generation.
	#       If a domain is residues 12-101, the codons that determine the dwell time at its positions are 13-102, while the codons
	#       that encode for this amino acid sequence are 12-101. Therefore, we will take codons 12-102 to facilitate both AA/DNA
	#       comparisons and translation rate profile generation. Note that when readProfile is generated for this hypothetical domain
	#       we consider only codons 13-102 (same is true of the DNA alignment step as well)

	DNAfasta = readFile(f_DNA_fasta)
	fullseq = ''
	domseq = ''
	n = 1

	for j in range (0, len(DNAfasta)):
		if '>' in DNAfasta[j].split()[0]:
			orf = DNAfasta[j].split()[0].strip('>')
			for key in domains:
				if orf == key.split(';')[0]:
					domain = domains[key]
					fullseq = ''
					domseq = ''
					n = 1
					# get the full DNA sequence for this ORF
					for k in range (j+1, j+100000):
						if DNAfasta[k].split() == []:
							continue
						#print (k, DNAfasta[k], DNAfasta[k].split())
						line = DNAfasta[k].split()[0]
						if '>' in line:
							break
						fullseq += line
					# make sure the sequence can be evenly split into codons
					if len(fullseq) % 3 != 0:
						print ('The DNA sequence extracted for '+orf+';'+str(domain.resStart)+'-'+str(domain.resEnd)+' is not evenly divisble into codons.')
						continue
					# get the specific portion corresponding to the domain of interest
					for i in range (0, len(fullseq), 3):
						if n >= (domain.resStart) and n <= (domain.resEnd+1): # for domain with residues 2-101, take codons 2-102
							domseq += (fullseq[i]+fullseq[i+1]+fullseq[i+2])
						if n > (domain.resEnd+1):
							break
						n += 1
					#print (key, 'DNAfull', fullseq)
					#print (key, 'DNAdom', domseq)
					# make sure the DNA sequence matches the length of the amino acid sequence plus 1 addition codon (see note in function header!)
					if len(domseq)/3 != (domain.resEnd-domain.resStart+2):
						print ('The DNA sequence extracted for '+orf+';'+str(domain.resStart)+'-'+str(domain.resEnd)+' does not match the length of the corresponding amino acid sequence + 1.')
						print ('Domain length = '+str(domain.resEnd-domain.resStart+1))
						print ('Number of DNA codons = '+str(len(domseq)/3))
						continue
					domains[key].addDNA(domseq)
	return

# function to get AA FASTA sequences for each domain
def getAAseq (domains, f_AA_fasta, ecoli=False):

	# domains: dictionary; keys are of the type 'orfID;resStart-resEnd', values are the corresponding Domain objects
	# f_AA_seq: str; path to AA FASTA file

	# N.B., this function does not apply an offset like getDNAseq does

	AAfasta = readFile(f_AA_fasta)
	fullseq = ''
	domseq = ''
	n = 1

	for j in range (0, len(AAfasta)):
		if '>' in AAfasta[j].split()[0]:
			if ecoli:
				orf = AAfasta[j].split()[0].strip('>').split('|')[3]
			else:
				orf = AAfasta[j].split()[0].strip('>')
			for key in domains:
				if orf == key.split(';')[0]:
					domain = domains[key]
					fullseq = ''
					domseq = ''
					n = 1
					# get the full AA sequence for this protein
					for k in range (j+1, j+100000):
						line = AAfasta[k].split()[0]
						if '>' in line:
							break
						fullseq += line
					# get the specific portion of this AA sequence corresponding to the current Domain object
					for i in range (0, len(fullseq)):
						if n >= domain.resStart and n <= domain.resEnd:
							domseq += fullseq[i]
						if n > domain.resEnd:
							break
						n += 1
					#print (key, 'AAfull', fullseq)
					#print (key, 'AAdom', domseq)
					# make sure the domain amino acid sequence matches the length of the domain
					if len(domseq) != (domain.resEnd-domain.resStart+1):
						print ('The AA sequence extracted for '+orf+';'+str(domain.resStart)+'-'+str(domain.resEnd)+' does not match the length of the domains.')
						print ('Domain length = '+str(domain.resEnd-domain.resStart+1))
						print ('Number of amino acids grabbed for sequence = '+str(len(domseq)))
						continue
					domains[key].addAAseq(domseq)
	return

# function to quality control DNA and AA sequences
def qcDNAandAAseqs (domains):

	# domains: dictionary; keys are of the type 'orfID;resStart-resEnd', values are the corresponding Domain objects

	# returns: Boolean; iff True, getDomains() function gives normal termination, otherwise an error is printed and the run dies

	# dictionary to map each DNA codon to the single-letter code for the amino acid for which it codes; '*' -> STOP codon
	mapDNAtoAA = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',}

	errors = {}
	error_count = 0
	for key in domains:
		domain = domains[key]
		aa = domain.AAseq
		dna = domain.DNAseq
		if aa is None or dna is None:
			continue
		# check their lengths
		if ((len(dna)/3)-1) != len(aa):
			print ('DNA and AA sequence lengths do not match for '+domain.orfID+';'+str(domain.resStart)+'-'+str(domain.resEnd))
			return False
		# make sure the sequences match one another
		j = 0
		for i in range (0, len(dna)-3, 3): # the '-3' is necessary to account for the +1 codon kept for DNA sequences
			if mapDNAtoAA[dna[i]+dna[i+1]+dna[i+2]] != aa[j]:
				if j == 0 and aa[j] == 'M': # this may be due to a non-AUG start codon; print an error but do not die
					print ('Possible non-ATG Start codon for', key, 'run will not be killed.')
				else:
					errors[error_count] = key+':'+dna[i]+dna[i+1]+dna[i+2]+'!='+aa[j]
					error_count += 1
			j += 1
	if len(errors) != 0:
		print (error_count, 'errors identified while quality controlling AA and DNA sequences.')
		for key in errors:
			print (key, errors[key])
		return True
	else:
		return True

# function to get raw A-site read profiles for each domain
def getReadProfiles (domains, f_asiteReads, normalize=True, frame=0, maxFracZeroReads=0.3, verbose=False):

	# domains: dictionary of Domain class objects; keys are 'orfID;resStart-resEnd' and values are corresponding Domain class objects
	# f_asiteReads: str, path to an A-site read file as output by Integer Programming method (Nabeel and O'Brien 2019
	# normalize: (optional, default=True) Boolean; default is True; if True, normalize readProfile by average number of reads in domain so average is 1.0 reads per codon position
	# frame: (optional, default=0) int; translation frame; 0 is the canonical frame and therefore the default, values of [0, 1, 21 are allowed values
	# maxFracZeroReads: (optional, default=0.3) float; the maximum proportion of positions within a domain allowed to have zero reads - if exceeded, no readProfile is added for the Domain object,
	#                   meaning that it is excluded from downstream analyses
	# verbose: (optional, default=False) Boolean, determines whether or not additional non-essential information is printed

	# returns: nothing, but modifies the Domain class objects in domains dictionary to include read profiles as possible
	#          based on the input A-site read profiles in f_asiteReads

	# N.B., this function applies an off-by-one adjustment to account for the fact that the translation time for a domain
	#       spanning residues 100-201 (for example) is dictated by the time required to decode codons 101-202, not 100-201.

	reads = readFile(f_asiteReads)
	readSeq = ''
	totalReads = 0
	codon = 1
	n = 0
	readProfile = None

	for i in range (0, len(reads)):
		for key in domains:
			# if we failed to get correct sequence information, do not bother parsing the read profile as we cannot use it anyways
			if domains[key].DNAseq is None or domains[key].AAseq is None:
				continue
			if reads[i].split()[0] == key.split(';')[0]:
				domain = domains[key]
				domains[key].addORFlength((int(reads[i].split()[1])-3)/3) # do not count the stop codon
				totalReads = np.float64(0.0)
				readSeq = reads[i].split()[2].split(',')
				readProfile = np.zeros((domain.resEnd-domain.resStart+1), dtype=np.float64)
				codon = 1
				n = 0
				# make sure the length of the raw read profile is evenly divisible by 3
				if len(readSeq) % 3 != 0:
					print ('A-site read sequence for '+domain.orfID+';'+str(domain.resStart)+'-'+str(domain.resEnd)+' is not evenly divisible by 3.')
					print ('This is potentially serious - check your A-site read file and input BAM file!')
					sys.exit()
				# get A-site reads for this domain
				for j in range (0, len(readSeq), 3):
					if codon >= (domain.resStart+1) and codon <= (domain.resEnd+1):
						totalReads += np.float64(readSeq[j+frame])
						readProfile[n] = np.float64(readSeq[j+frame])
						n += 1
					elif codon > (domain.resEnd+1):
						break
					elif codon < (domain.resStart+1):
						pass
					else:
						print ('Failed to correctly parse read sequence for '+domain.orfID+';'+str(domain.resStart)+'-'+str(domain.resEnd))
					codon += 1
				# make sure we counted the correct number of codons for this domain
				if n != (domain.resEnd-domain.resStart+1):
					if verbose:
						print ('Incorrect number of codons counted for '+domain.orfID+';'+str(domain.resStart)+'-'+str(domain.resEnd))
						print (n, 'codons were counted for this domain')
					continue
				# if no reads were found for this domain skip to the next one
				if totalReads == np.float64(0.0):
					continue
				# check the number of domain positions that have no reads, and skip this domain if the proportion is too high
				numZero = 0.0
				for m in range (0, len(readProfile)):
					if readProfile[m] == np.float64(0.0):
						numZero += 1.0
				if numZero/float(len(readProfile)) > maxFracZeroReads:
					if verbose:
						print (key, 'has insufficient read coverage and will be excluded.')
					continue
				if verbose:
					temp = open('new_output/rawReadProfiles/'+domains[key].orfID+'_'+str(domains[key].resStart)+'-'+str(domains[key].resEnd)+'.txt', 'w')
					for z in range (0, len(readProfile)):
						temp.write('%.0f' %(domains[key].resStart+z)+'\t'+'%.9f' %readProfile[z]+'\n')
				# normalize if requested by the user, adjusting values so the overall average for the domain is 1.0 reads per codon
				if normalize:
					readProfile = readProfile/(totalReads/np.float64(len(readProfile)))
				# update this Domain object with its readProfile
				domain.addReadProfile(readProfile)
	return

# function to generate "single-domain translation rate profiles" for domains
def getSingleDomainProfiles(domains, verbose=False):

	# domains: dictionary of Domain class objects; keys are of the form 'YEL066W;32-173', values are Domain class objects
	# verbose: (optional, default=False) Boolean; whether or not to print non-essential information to screen

	# returns: nothing, but modifies the Domain class objects in domains by adding single domain read profiles as possible

	for key in domains:
		domain = domains[key]
		if domain.AAseq is None or domain.DNAseq is None or domain.readProfile is None:
			continue
		z = domain.resStart+1 # codon whose reads are placed at first position in readProfile
		x = 1 # counter for relative domain position
		pos = np.array([], dtype=int) # list that will contain positions to be deleted from the profile
		relativeDomainPosition = np.zeros((len(domain.readProfile)), dtype=np.float64)

		if verbose:
			print (domains[key].orfID+'_'+str(domains[key].resStart)+'-'+str(domains[key].resEnd)+' ORFlength='+str(domains[key].ORFlength))

		# find positions corresponding to the first 40 or last 20 positions in the gene
		for i in range (0, len(domain.readProfile)):
			if z <= 40: # we want to remove these positions due to biases from initiation
				pos = np.append(pos, i)
			elif z > (domain.ORFlength-20): # and we want to remove these due to biases from termination
				pos = np.append(pos, i)
			else:
				pass
			relativeDomainPosition[i] = x
			z += 1
			x += 1

		if verbose:
			print (domain.orfID, domain.resStart, domain.resEnd)
			for i in range (0, len(domain.readProfile)):
			       print (domain.resStart+i+1, domain.readProfile[i], relativeDomainPosition[i])

			print ('Positions to delete:', pos)

		# make a copy of the read profile
		readProfileCopy = np.zeros((len(domain.readProfile)), dtype=np.float64)
		for i in range (0, len(readProfileCopy)):
			readProfileCopy[i] = domain.readProfile[i]

		# delete domain positions corresponding to the first 40 or last 20 codons from both readProfileCopy and
		# from relativeDomainPosition
		readProfileCopy = np.delete(readProfileCopy, pos)
		relativeDomainPosition = np.delete(relativeDomainPosition, pos)

		if verbose:
			print ('After deleting first 40 and last 20 positions:')
			for i in range (0, len(relativeDomainPosition)):
			       print (i, relativeDomainPosition[i], readProfileCopy[i])

		pos = np.array([], dtype=int) # set pos to be an empty array once more

		# make an array of indices of readProfileCopy
		indices = np.zeros((len(readProfileCopy)), dtype=int)
		for i in range (0, len(indices)):
			indices[i] = i

		# find all front positions equal to zero before the first non-zero position in readProfileCopy
		for i in range (0, len(readProfileCopy)):
			if readProfileCopy[i] == np.float64(0.0):
				pos = np.append(pos, i)
			elif readProfileCopy[i] > np.float64(0.0):
				break
			else:
				print ('Error - reads should never be negative.')
				sys.exit()

		# find all back positions equal to zero before the first non-zero position (counting from the end)
		for i in range (1, len(readProfileCopy)):
			if readProfileCopy[-i] == np.float64(0.0):
				pos = np.append(pos, -i)
			elif readProfileCopy[-i] > np.float64(0.0):
				break
			else:
				print ('Error - reads should never be negative.')
				sys.exit()

		# if we have found something that needs deleted
		if pos.size > 0:
			readProfileCopy = np.delete(readProfileCopy, indices[pos]) # select pos positions from indices array and delete - this is necessary because np.delete ignores negative obj indices
			relativeDomainPosition = np.delete(relativeDomainPosition, indices[pos])

		if verbose:
			print ('After trimming leading and trailing zeros:')
			for i in range (0, len(relativeDomainPosition)):
			       print (i, relativeDomainPosition[i], readProfileCopy[i])

		# in this new, trimmed readProfileCopy, make sure that we have a reasonable profile length
		if len(readProfileCopy) < 64: # we want the smoothed profile to have at minimum 50 positions
			print ('Insufficient length in trimmed singleDomainProfile for:', domain.orfID, domain.resStart, domain.resEnd)
			continue

		# now, spline the trimmed read profile so no positions have a zero reads
		splinedProfile = splineProfile(readProfileCopy)

		if splinedProfile is None:
			print ('Failed to produce splined singleDomainProfile - likely due to insufficient non-zero positions - for:', domain.orfID, domain.resStart, domain.resEnd)
			continue

		if verbose:
			print ('After splining:')
			for i in range (0, len(relativeDomainPosition)):
			       print (i, relativeDomainPosition[i], readProfileCopy[i], splinedProfile[i])

		# then smooth it - also, need to be careful to keep track of the bin centers for domain alignment later
		smoothedSplinedProfile = running_mean(splinedProfile, 15)
		relativeDomainPosition = running_mean(relativeDomainPosition, 15)

		if verbose:
			print ('After smoothing:', len(relativeDomainPosition))
			for i in range (0, len(relativeDomainPosition)):
			       print (i, relativeDomainPosition[i], smoothedSplinedProfile[i])

		if verbose:
			temp = open('new_output/rawSingleDomainProfiles/'+domains[key].orfID+'_'+str(domains[key].resStart)+'-'+str(domains[key].resEnd)+'.txt', 'w')
			for i in range (0, len(relativeDomainPosition)):
				temp.write('%.2f' %relativeDomainPosition[i]+'\t'+'%.9f' %smoothedSplinedProfile[i]+'\n')
			temp.close()

		# place the domain indices and reads into the same 2-column array
		if len(smoothedSplinedProfile) != len(relativeDomainPosition):
			print ('ERROR: normedSmoothedSplinedProfile and relativeDomainPosition must have the same length!')
			sys.exit()

		finalProfile = np.zeros((len(relativeDomainPosition),2), dtype=np.float64)
		finalProfile[:,0] = relativeDomainPosition
		finalProfile[:,1] = smoothedSplinedProfile

		domains[key].addSingleDomainProfile(finalProfile)
		print ('Correctly generated singleDomainProfile for', key)

	return

# function to extract pairs of evolutionarily related domains at either the family or superfamily level
def getHomologs (set1, set2, maxSizeDiff=25, level='family', requireReadProfile=True, verbose=False):

	# set1: a dictionary of Domain class objects
	# set2: a second dictionary of Domain class objects; may be identical to or different from set1 depending on what you are trying to do
	# maxSizeDiff: (optional, default=25) int; the maximum size different allowed between a pair of domains
	# level: (optional, default=family) str, the "level" at which homologs will be determined; choose from ['superfamily', 'family']
	# requireReadProfile: (optional, default=True) Boolean, iff True make sure that related domains have read profiles (i.e., domain.readProfile != None)
	# verbose: (optional, default=False) Boolean, determines whether or not additional non-essential information is printed

	# returns: a list of lists, each list being a pair of domain identifies of the type 'orfID;resStart-resEnd' (one from domainList1 and the other from domainList2) that are related at the requested level

	# N.B., this function requires that a domain have both a DNAseq and a readProfile in order for it to be included in downstream analyses

	if level not in ['superfamily', 'family']:
		print (level, 'is not an allowed value for the parameter \'level\' - choose from superfamily and family.')
		sys.exit()

	pairs = []
	Npairs = 0
	unique = []

	for d1 in set1:
		# make sure this Domain object has a readProfile, DNAseq, AAseq, and singleDomainProfile
		if set1[d1].readProfile is None and requireReadProfile:
			continue
		if set1[d1].DNAseq is None or set1[d1].AAseq is None:
			continue
		if set1[d1].singleDomainProfile is None:
			continue
		for d2 in set2:
			if d1 == d2:
				continue
			if set2[d2].readProfile is None and requireReadProfile:
				continue
			if set2[d2].DNAseq is None or set2[d2].AAseq is None:
				continue
			if set2[d2].singleDomainProfile is None:
				continue
			if np.abs((set1[d1].resEnd-set1[d1].resStart+1) - (set2[d2].resEnd-set2[d2].resStart+1)) > maxSizeDiff:
				continue
			if d1.split(';')[0] == d2.split(';')[0]: # check for overlapping domain definitions on the same gene
				r1temp = range(set1[d1].resStart, set1[d1].resEnd+1)
				r2temp = range(set2[d2].resStart, set2[d2].resEnd+1)
				r1 = set(r1temp)
				if len(r1.intersection(r2temp)) != 0:
					if verbose:
						print ('The domains '+set1[d1].orfID+';'+str(set1[d1].resStart)+'-'+str(set1[d1].resEnd)+' and '+\
        	                                       set2[d2].orfID+';'+str(set2[d2].resStart)+'-'+str(set2[d2].resEnd)+' are in the same gene and overlap - will be ignored.')
					continue
			if level == 'superfamily':
				if set1[d1].superfamily == set2[d2].superfamily:
					if [d1, d2, ] not in pairs and [d2, d1, ] not in pairs:
						pairs.append( [ d1, d2, ] )
						if set1[d1].superfamily not in unique:
							unique.append(set1[d1].superfamily)
			elif level == 'family':
				if set1[d1].family == set2[d2].family:
					if [d1, d2, ] not in pairs and [d2, d1, ] not in pairs:
						pairs.append( [ d1, d2, ] )
						if set1[d1].family not in unique:
							unique.append(set1[d1].family)
			else:
				print (level, 'is not an allowed value for the parameter \'level\' - choose from superfamily and family.')
				sys.exit()
	if verbose:
		if level == 'family':
			print (len(pairs), 'pairs of related domains within', len(unique), 'unique families identified.')
		elif level == 'superfamily':
			print (len(pairs), 'pairs of related domains within', len(unique), 'unique superfamilies identified.')
	return pairs

# function to get random pairs of domains; the same size selection criteria are used as for selecting pairs of homologs, but there is no requirement
# that pairs of domains be in the same SUPERFAM family (though it can happen by random chance)
def getRandomDomainPairs(set1, set2, Npairs, maxAttempts=5000000, maxSizeDiff=25, requireReadProfile=True, verbose=False):

	# set1: a dictionary of Domain class objects
	# set2: a second dictionary of Domain class objects; may be identical to or different from set1 depending on what you are trying to do
	# Npairs: int; the number of random pairs of domains to select; random selections will be made until this number of reasonable pairs are found or until the number of attemps == maxAttempts
	# maxAttemps: (optional, default=10000) int; the number of times to try different random selections
	# maxSizeDiff: (optional, default=25) int; the maximum size different allowed between a pair of domains
	# verbose: (optional, default=False) Boolean; determines whether or not additional non-essential information is printed

	# returns: a list of lists, each list being a pair of domain identifies of the type 'orfID;resStart-resEnd' (one from domainList1 and the other from domainList2)

	# N.B., this function does not prevent pairs of related domains from being chosen, it just performs totally random selection of pairs of domains (I think this makes more
	#       sense than intentionally excluding them but it would be easy to add an extra condition to control for this, something like "if d1.family == d2.family:\n continue" would do the trick)

	pairs = []
	attempts = 0

	d1_list = [] # make lists of dictionary keys for sets 1 and 2
	for d1 in set1:
		d1_list.append(d1)
	d2_list = []
	for d2 in set2:
		d2_list.append(d2)

	run = True

	while run == True:
		r1 = np.random.randint(0, len(d1_list)) # get random integer
		d1 = set1[d1_list[r1]] # use it to get random domain
		r2 = np.random.randint(0, len(d2_list)) # get a second random integer
		d2 = set2[d2_list[r2]] # and a use it to get a second random domain

		# make sure the pair passes certain criteria
		if d1.readProfile is None or d2.readProfile is None: # must have readProfile
			continue
		if d1.singleDomainProfile is None or d2.singleDomainProfile is None: # and singleDomainProfile
			continue
		if d1.DNAseq is None or d2.DNAseq is None: # and DNA sequence
			continue
		if d1.AAseq is None or d2.AAseq is None: # and AA sequence
			continue
		if np.abs(( d1.resEnd-d1.resStart+1) - (d2.resEnd-d2.resStart+1)) > maxSizeDiff: # and have roughly the same size within tolerance = maxSizeDiff
			continue
		pairs.append([d1_list[r1], d2_list[r2]]) # passing all these checks, accept the random pair and add a list to the list of good pairs

		# exit conditions
		if len(pairs) == Npairs:
			run = False
		if attempts == maxAttempts:
			run = False
		attempts += 1

	# only return the random list of pairs if a sufficient number were found
	if len(pairs) != Npairs:
		print ('Failed to select', Npairs, 'random pairs of domains.')
		sys.exit()
	else:
		return pairs

# function to find the list of domains between two datasets that have the same domain definition and acceptable read profiles
# this is analogous to getHomologs() but, rather than enable comparisons between translation rate profiles for related pairs of domains, it enables
# comparisons between translation rate profiles in the same domain within one dataset (simple control) or between two datasets from the same organism
def getMatching(set1, set2, verbose=False):

	# set1: a dictionary of Domain class objects; keys have the form 'YEL066W;32-173'
        # set2: a second dictionary of Domain class objects

	# returns: list of lists in the same format as the list of lists output by getHomologs

	# N.B., run the result of this function through the same set of analyses as is performed for comparisons between related domains

	pairs = []
	for key in set1:
		if set1[key].singleDomainProfile is not None and set2[key].singleDomainProfile is not None: # just make sure we have a singleDomainProfile for both domains!
			pairs.append([ key, key, ])
	return pairs

# function to perform DNA sequence alignment between pairs of domains
def alignDomains (set1, set2, pairs, outdir, multi=True, nprocs=10, verbose=False):

	# set1: dictionary of Domain class objects with keys of the type 'orfID;resStart-resEnd' as generated by getDomains()
	# set2: a second dictionary of Domain class objects with keys of the type 'orfID;resStart-resEnd' as generate by getDomains()
	# pairs: list of lists, where each inner list has two elements: [ ..., [domain from set1, domain from set2], ... ]; for example
	#       [ ..., ['YEL066W;32-173', 'YPR051W;13-159'], ... ]
	# multi: (optional, default=True) Boolean, iff True use multiprocessing to speed up alignments
	# nprocs: (optional, default=10) int, the number of processors to use for multiprocessing. If multi==False, has no effect
	# verbose: (optional, default=False) Boolean, determines whether or not additional non-essential information is printed

	# returns: nothing, but calls runMUSCLE to generate both DNA and AA alignments

	# N.B., DNA alignments are placed in 'dna_alignments/' and amino acid alignments are placed in 'aa_alignments/'
	#       DNA sequence identity is used for selecting related domains, while amino acid alignments are used for generating
	#       aligned translation rate profiles

	if outdir[-1] != '/':
		outdir += '/'

	# do DNA and AA alignmentsalignments
	mkdir(outdir+'dna_alignments/')
	mkdir(outdir+'aa_alignments/')
	cmds = []
	for pair in pairs:
		d1, d2 = pair[0:2]
		cmds.append( [set1[d1], set2[d2], outdir,] )
	if multi:
		pool = mp.Pool(nprocs)
		results = [pool.apply_async( runMUSCLE, cmd ) for cmd in cmds]
		for result in results:
			result.get()
	else:
		for cmd in cmds:
			runMUSCLE(cmd[0], cmd[1], cmd[2])
			sys.exit()
	return

# function to write a FASTA format DNA sequence file given a pair of Domain class objects
def runMUSCLE (domain1, domain2, outdir):

	# domain1: first Domain class object
	# domain2: second Domain class object

	# returns: nothing, but writes input FASTA and output alignment files to dna_alignments/ and aa_alignments/
	if outdir[-1] != '/':
		outdir += '/'

	# generate the input DNA FASTA file and run its alignment with MUSCLE
	ofile = open(outdir+'dna_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.fa', 'w')
	ofile.write('>'+domain1.orfID+';'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'\n')
	ofile.write(domain1.DNAseq+'\n')
	ofile.write('>'+domain2.orfID+';'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'\n')
	ofile.write(domain2.DNAseq+'\n')
	ofile.close()
	subprocess.run('./src/muscle3.8.31_i86linux64 -seqtype auto -msf '+\
                       '-in '+outdir+'dna_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.fa'+\
                       ' -out '+outdir+'dna_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.afa',
                       capture_output=True, shell=True)

	# generate the input AA FASTA file and run alignment with MUSCLE
	ofile = open(outdir+'aa_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.fa', 'w')
	ofile.write('>'+domain1.orfID+';'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'\n')
	ofile.write(domain1.AAseq+'\n')
	ofile.write('>'+domain2.orfID+';'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'\n')
	ofile.write(domain2.AAseq+'\n')
	ofile.close()
	subprocess.run('./src/muscle3.8.31_i86linux64 -seqtype auto -msf '+\
                       '-in '+outdir+'aa_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.fa'+\
                       ' -out '+outdir+'aa_alignments/'+domain1.orfID+'_'+str(domain1.resStart)+'-'+str(domain1.resEnd)+'_'+domain2.orfID+'_'+str(domain2.resStart)+'-'+str(domain2.resEnd)+'.afa',
                       capture_output=True, shell=True)
	return

# function to get DNA sequence identity in order to determine which pairs of domains are usable for analysis
def findAcceptableDomainPairs (pairs, outdir, min=0.30, max=0.80, verbose=False):

	# pairs: a list of lists; each sublist contains two domain identifiers of the type 'orfID;resStart-resEnd'
	# min: (optional, default=0.30) float; minimum sequence identity to be kept
	# max: (optional, default=0.80) float; maximum sequence identity to be kept
	# verbose: (optional, default=False) Boolean; determines whether or not additional non-essential information is printed

	# returns: a list of the pairs that have a DNA identity greater than or equal to ID_limit

	good_pairs = []
	if outdir[-1] != '/':
		outdir += '/'
	if verbose:
		ofile = open(outdir+'dna_alignments/dna_sequence_ID.dat', 'w')

	for pair in pairs:
		d1, d2 = pair[0:2]
		f_alignment = outdir+'dna_alignments/'+d1.split(';')[0]+'_'+d1.split(';')[1].split('-')[0]+'-'+d1.split(';')[1].split('-')[1]+'_'+\
                                              d2.split(';')[0]+'_'+d2.split(';')[1].split('-')[0]+'-'+d2.split(';')[1].split('-')[1]+'.afa'
		if os.path.exists(f_alignment):
			alignment = readFile(f_alignment)
			seq1 = ''
			seq2 = ''
			found = False
			x = 0
			same = 0
			for i in range (0, len(alignment)):
				line = alignment[i].split()
				if line == []:
					continue
				if '//' == line[0]:
					found = True
					continue
				if found:
					if x == 0:
						for j in range (1, len(line)):
							seq1 += line[j]
					elif x == 1:
						for j in range (1, len(line)):
							seq2 += line[j]
					else:
						print ('Failed to parse '+f_alignment+' correctly; the counter x should only ever be 0 or 1!')
						sys.exit()
					if x == 0:
						x = 1
					elif x == 1:
						x = 0
					else:
						print ('Failed to parse '+f_alignment+' correctly; the counter x should only ever be 0 or 1!')
						sys.exit()
			if len(seq1) != len(seq2):
				print ('The length of the two aligned sequences plus gaps should always be equal!')
				sys.exit()
			for k in range (0, len(seq1)):
				if seq1[k] == '.' or seq2[k] == '.': # do not count gaps
					continue
				if seq1[k] == seq2[k]:
					same += 1
				if verbose:
					print (d1, d2, k+1, seq1[k], seq2[k], same)
			if verbose:
				ofile.write((d1+' : '+d2)+', '+'%.3f' %(float(same)/float(len(seq1))*100.0)+'%'+'\n')
			if float(same)/float(len(seq1)) >= min and float(same)/float(len(seq1)) <= max:
				good_pairs.append(pair)
	if verbose:
		ofile.close()

	return good_pairs

# function to generate aligned translation rate profiles based on amino acid MUSCLE alignments
def getAlignedReadProfiles (pairs, domains1, outdir, domains2=None, maxGapSize=5, maxNumGaps=10,  verbose=False):

	# pairs: list of lists; each sublist is of the form '[ ..., ['YEL066W;32-173', 'YPR051W;13-159'], ... ]'
	#        this list should be generated by passing the list generated by getHomologs through findAcceptableDomainPairs after performing DNA and AA alignments
	# domains1: dictionary of Domain class objects (needed to get access to read profiles)
	# domains2: (optional, default= None) a second dictionary of Domain class objects; if provided, assume that d2 in pairs has a read profile stored in domains2 dictionary
	#           When this argument is left to a default of None, assume we are doing a comparison within one ribosome profiling experiment
	# maxGapSize: (optional, default=5) int; the maximum gap size allowed for a usable alignment
	# maxNumGaps: (optional, default=10) int; the maximum number of gaps (of size 1 or greater) allowed for a usable alignment

	# returns: dictionary with a pair of domains as the key, e.g. 'YEL066W;32-173:YPR051W;13-159'
	#          value is an np.array representing the aligned and processed translation rate profiles
	if outdir[-1] != '/':
		outdir += '/'
	alignedProfiles = {}
	for pair in pairs:

		d1, d2 = pair[0:2]

		f_alignment = outdir+'aa_alignments/'+d1.split(';')[0]+'_'+d1.split(';')[1].split('-')[0]+'-'+d1.split(';')[1].split('-')[1]+'_'+\
                                               d2.split(';')[0]+'_'+d2.split(';')[1].split('-')[0]+'-'+d2.split(';')[1].split('-')[1]+'.afa'
		alignment = readFile(f_alignment)

		# these variables need to be reset for each new pair of domains
		alignmentLength = 0
		inGap = False
		inEndGap = False
		numGaps = 0 # total number of gaps
		numEndGaps = 0 # number of end gaps; == 0, 1, or 2
		gapSizes = [] # list of gap sizes
		gapSize = 0
		maxGapSizeFound = 0 # size of largest gap
		alignmentMatrix = None
		gapVector = None
		found = False
		row = 0
		seq1count = 0
		seq2count = 0
		c = 0
		d1ReadMatrix = None # reads per aligned position for domain 1
		d2ReadMatrix = None # reads per aligned position for domain 2
		d1count = 0
		d2count = 0
		section = ''

		print ('Processing alignment for:', d1, d2)
		for line in alignment:
			if 'Len' in line:
				alignmentLength = int(line.split()[3])
				#print (f_alignment, alignmentLength)
				alignmentMatrix = np.zeros((alignmentLength, 2), dtype=str)
				gapVector = np.ones((alignmentLength), dtype=int)
			if line == '\n':
				continue
			if line == '//\n':
				found = True
				continue
			if found == True:
				line = line.split()
				for i in range (1, len(line)):
					section = line[i].strip()
					if row == 0:
						for s in section:
							alignmentMatrix[seq1count, row] = s
							seq1count += 1
					elif row == 1:
						for s in section:
							alignmentMatrix[seq2count, row] = s
							seq2count += 1
					else:
						print ('Failed to correctly determine the row in alignment file', f_alignment)
						sys.exit()
				if row == 0:
					row = 1
				elif row == 1:
					row = 0
				else:
					print ('Failed to correctly modify the value of \'row\' while generating alignmentMatrix for', f_alignment)
					sys.exit()
		# generate the gapVector for this alignment
		for i in range (0, len(alignmentMatrix)):
			if alignmentMatrix[i,0] == '.':
				gapVector[i] = 0
			if alignmentMatrix[i,1] == '.':
				gapVector[i] = 0

		# determine gap sizes, number of intrasequence gaps, and number of end gaps
		for j in range (0, len(gapVector)):
			# if you have found a new gap
			if gapVector[j] == 0 and inGap == False:
				inGap = True
				numGaps += 1
				gapSize += 1
				if j == 0 or j == (len(gapVector)-1):
					inEndGap    = True
					numEndGaps += 1
			# if you have found another element in a gap for the current gap
			elif gapVector[j] == 0 and inGap == True:
				gapSize += 1
				if j == (len(gapVector)-1):
					numEndGaps += 1
			# if you have found the first non-gap element (i.e., the end of the gap)
			elif gapVector[j] == 1 and inGap == True:
				inGap = False
				if inEndGap == True:
					inEndGap = False
				else:
					gapSizes.append(gapSize)
				gapSize = 0
			else:
				pass
		# determine the maximum gap size if there were any non-end gaps
		if gapSizes != []:
			maxGapSizeFound = np.max(gapSizes)
		# make sure we have a physically reasonable number of end gaps
		if numEndGaps not in [0, 1, 2]:
			print ('An alignment can have at max 2 end gaps - a value of', numEndGaps, 'is not realistic.')
			sys.exit()
		# determine whether or not we have an acceptable alignment based on the number and size of identified gaps
		if maxGapSizeFound > maxGapSize: # check max gap size
			if verbose:
				print ('Failed to accept alignment; max gap size exceeded for'+ d1 + ':'+ d2 + '; '+str(maxGapSizeFound)+' > '+str(maxGapSize))
			continue
		if (numGaps - numEndGaps) > maxNumGaps: # check total number of non-end gaps
			if verbose:
				print ('Failed to accept alignment; total number of allowed gaps exceeded for'+ d1 + ':'+ d2 + '; '+str(numGaps - numEndGaps)+' > '+str(maxNumGaps))
			continue

		if verbose:
			print (d1, d2, 'gaps:', gapSizes, '; the largest gap is', maxGapSizeFound, 'positions')

		if verbose:
			for z in range (0, len(alignmentMatrix)):
				print (d1, d2, z+1, alignmentMatrix[z,:])

		# If we have reached this point, the alignment is "accepted" and we can proceed to generate the aligned read profiles
		d1ReadMatrix = np.zeros((alignmentLength), dtype=np.float64)
		d2ReadMatrix = np.zeros((alignmentLength), dtype=np.float64)

		# We need to differentiate gaps from positions that simply have no reads - this is critical for processAlignedProfiles() function
		# if a position is in a gap assign it a value of -1

		if domains2 is None:
			readProfile1 = domains1[d1].readProfile
			readProfile2 = domains1[d2].readProfile
		else:
			readProfile1 = domains1[d1].readProfile
			readProfile2 = domains2[d2].readProfile

		for i in range (0, len(alignmentMatrix)):
			if alignmentMatrix[i,0] != '.':
				d1ReadMatrix[i] = readProfile1[d1count]
				d1count += 1
			else:
				d1ReadMatrix[i] = -1.0
			if alignmentMatrix[i,1] != '.':
				d2ReadMatrix[i] = readProfile2[d2count]
				d2count += 1
			else:
				d2ReadMatrix[i] = -1.0
		if verbose:
			for z in range (0, len(d1ReadMatrix)):
				print (d1, d2, z+1, d1ReadMatrix[z], d2ReadMatrix[z])

		if domains2 is not None:
			processedProfile1, processedProfile2 = processAlignedReadProfiles(d1, d1ReadMatrix, d2, d2ReadMatrix, domains1, domains2=domains2, verbose=verbose)
		else:
			processedProfile1, processedProfile2 = processAlignedReadProfiles(d1, d1ReadMatrix, d2, d2ReadMatrix, domains1, verbose=verbose)

		if processedProfile1 is None or processedProfile2 is None:
			if verbose:
				print ('Failed to correctly process aligned rate profiles for', d1, d2)
			continue
		alignedProfiles[d1+':'+d2] = [processedProfile1, processedProfile2, ]

		print ('Aligned rate profiles for', d1, d2, 'correctly processed.')

	return alignedProfiles

# function to perform naive domain profile alignments for randomly selected profiles
def alignRandomReadProfiles(pairs, domains1, domains2=None):

	# pairs: list of lists; each sublist is of the form '[ ..., ['YEL066W;32-173', 'YPR051W;13-159'], ... ]'
	#        this list does NOT need to be passed through alignment steps as the aligned profiles will be generated based only on the singleDomainProfiles
	# domains1: dictionary of Domain class objects (needed to get access to singleDomain read profiles)
	# domains2: (optional, default= None) a second dictionary of Domain class objects; if provided, assume that d2 in pairs has a read profile stored in domains2 dictionary
	#           When this argument is left to a default of None, assume we are doing a comparison within one ribosome profiling experiment

	# returns: dictionary with a pair of domains as the key, e.g. 'YEL066W;32-173:YPR051W;13-159'
	#          value is an np.array representing the aligned and processed singleDomainProfiles. The length of the aligned profiles is taken to be the length of the
	#          smaller of the two profiles in each pair

	alignedProfiles = {}
	for pair in pairs:
		d1, d2 = pair[0:2]
		if domains2 != None:
			length = min([len(domains1[d1].singleDomainProfile), len(domains2[d2].singleDomainProfile)])
			temp = alignSingleDomainProfiles(domains1[d1].singleDomainProfile, domains2[d2].singleDomainProfile, length)
			if temp[0] is None or temp[1] is None:
				continue
		else:
			length = min([len(domains1[d1].singleDomainProfile), len(domains1[d2].singleDomainProfile)])
			temp = alignSingleDomainProfiles(domains1[d1].singleDomainProfile, domains1[d2].singleDomainProfile, length)
			if temp[0] is None or temp[1] is None:
				continue
		alignedProfiles[d1+':'+d2] = [ normProfile(temp[0][:,1]), normProfile(temp[1][:,1]), ]
	return alignedProfiles

# function to process aligned read profiles into their final form for comparison
def processAlignedReadProfiles (d1, d1ReadMatrix, d2, d2ReadMatrix, domains1, domains2=None, verbose=False):

	# d1: str; identifier for the first domain of the type 'YEL066W;32-173'
	# d1ReadMatrix: np.array((), dtype=np.float64); reads per aligned position for d1 in the d1/d2 alignment
	# d2: str; identifier for the second domain of the type 'YPR051W;13-159'
	# d2ReadMatrix: np.array((), dtype=np.float64); reads per aligned position for d2 in the d1/d2 alignment
	# domains1: dictionary; values are domain identifiers of the type YEL066W;32-173, values are corresponding Domain class objects
	# domains2: optional second dictionary; if provided, comparisons will be made between information in domains1 and domains2
	# verbose: (optional, default=False) Boolean; if True, print some extra information to file and to screen

	# returns: the processed aligned translation rate profiles

	#print (d1, domains1[d1].ORFlength, d2, domains1[d2].ORFlength)

	if domains2 is None:
		z1 = domains1[d1].resStart+1 # counter for domain position in domain 1
		z2 = domains1[d2].resStart+1 # counter for domain position in domain 2
	else:
		z1 = domains1[d1].resStart+1
		z2 = domains2[d2].resStart+1
	pos1 = np.array([], dtype=int) # array of positions to delete from domain 1 array
	pos2 = np.array([], dtype=int) # array of positions to delete from domain 2 array

	# first, identify end gap positions at the front for d1 and d2
	for i in range (0, len(d1ReadMatrix)):
		if d1ReadMatrix[i] != np.float64(-1.0):
			break
		else:
			pos1 = np.append(pos1, i)
	for i in range (0, len(d2ReadMatrix)):
		if d2ReadMatrix[i] != np.float64(-1.0):
			break
		else:
			pos2 = np.append(pos2, i)

	if verbose:
		print (d1, 'front end gaps:', pos1)
		print (d2, 'front end gaps:', pos2)

	# and also end gaps at the back for d1 and d2
	for i in range (1, len(d1ReadMatrix)):
		if d1ReadMatrix[-i] != np.float64(-1.0):
			break
		else:
			pos1 = np.append(pos1, -i)
	for i in range (1, len(d2ReadMatrix)):
		if d2ReadMatrix[-i] != np.float64(-1.0):
			break
		else:
			pos2 = np.append(pos2, -i)

	if verbose:
		print (d1, 'back end gaps:', pos1)
		print (d2, 'back end gaps:', pos2)

	# now, identify positions in the matrices corresponding to gene positions prone to biases due to initiation and termination
	if domains2 is not None:
		check = domains2[d2].ORFlength
	else:
		check = domains1[d2].ORFlength
	for i in range (0, len(d1ReadMatrix)):
		if d1ReadMatrix[i] != np.float64(-1): # check domain 1
			if z1 <= 40:
				pos1 = np.append(pos1, i)
			elif z1 > (domains1[d1].ORFlength-20):
				pos1 = np.append(pos1, i)
			else:
				pass
			z1 += 1
		if d2ReadMatrix[i] != np.float64(-1): # check domain 2
			if z2 <= 40:
				pos2 = np.append(pos2, i)
			elif z2 > (check - 20):
				pos2 = np.append(pos2, i)
			else:
				pass
			z2 += 1

	if verbose:
		print (d1, 'all positions to pass to zero:', pos1)
		print (d2, 'all positions to pass to zero:', pos2)

	for pos in pos1:
		d1ReadMatrix[pos] = 0.0
	for pos in pos2:
		d2ReadMatrix[pos] = 0.0
	for i in range (0, len(d1ReadMatrix)):
		if d1ReadMatrix[i] == np.float64(-1):
			d1ReadMatrix[i] = 0.0
		if d2ReadMatrix[i] == np.float64(-1):
			d2ReadMatrix[i] = 0.0

	readProfile1, readProfile2 = trimReadProfiles(d1ReadMatrix, d2ReadMatrix)

	if readProfile1 is None or readProfile2 is None:
		return None, None

	if len(readProfile1) < 64: # note that trimReadProfiles already enforces len(readProfile1) == len(readProfile2)
		return None, None

	# spline to cover regions with 0 reads
	splinedProfile1 = splineProfile(readProfile1)
	splinedProfile2 = splineProfile(readProfile2)

	if splinedProfile1 is None or splinedProfile2 is None:
		return None, None

	if verbose:
		for i in range (0, len(splinedProfile1)):
			print ('splined:', d1, d2, splinedProfile1[i], splinedProfile2[i])

	# smooth with a 15-position running average
	smoothedSplinedProfile1 = running_mean(splinedProfile1, 15)
	smoothedSplinedProfile2 = running_mean(splinedProfile2, 15)

	if verbose:
		for i in range (0, len(smoothedSplinedProfile1)):
			print ('smoothed:', d1, d2, smoothedSplinedProfile1[i], smoothedSplinedProfile2[i])

	# normalize to have an area under the curve of 1.0
	normedSmoothedSplinedProfile1 = normProfile(smoothedSplinedProfile1)
	normedSmoothedSplinedProfile2 = normProfile(smoothedSplinedProfile2)

	if verbose:
		for i in range (0, len(normedSmoothedSplinedProfile1)):
			print ('normalized:', d1, d2, normedSmoothedSplinedProfile1[i], normedSmoothedSplinedProfile2[i])

	if verbose:
		temp = open('new_output/alignedParalogProfiles/'+d1.split(';')[0]+'_'+d1.split(';')[1]+'_'+d2.split(';')[0]+'_'+d2.split(';')[1]+'.txt', 'w')
		for i in range (0, len(normedSmoothedSplinedProfile1)):
			temp.write(str(i+1)+'\t'+'%.9f' %normedSmoothedSplinedProfile1[i]+'\t'+'%.9f' %normedSmoothedSplinedProfile2[i]+'\n')
		temp.close()

	return normedSmoothedSplinedProfile1, normedSmoothedSplinedProfile2

# given two input arrays, trim zeros from the ends of either so that they end up the same length
def trimReadProfiles(profile1, profile2):

	# profile1: np.array ((), dtype=np.float64); input rate profile 1
	# profile2: np.array ((), dtype=np.float64); input rate profile 2

	# return  : trimmed versions of profiles 1 and 2. Trimming is done such that the profiles maintain the
	#           same total length. If profile1 has a zero at positions 1, 2, 3 then positions 1, 2, 3 are deleted
	#           from profiles 1 and 2.

	#for i in range (0, len(profile1)):
	#	print ('original', '\t', '%.9f' %profile1[i]+'\t'+'%.9f' %profile2[i])

	# trim the front of profile1
	frontTrimProfile1 = np.trim_zeros(profile1, 'f')

	if frontTrimProfile1.size == 0:
		return None, None

	# if profile1 was trimmed then trim profile2 to match
	if len(frontTrimProfile1) != len(profile1):
		profile2 = profile2[-len(frontTrimProfile1):]

	if profile2.size == 0:
		return None, None

	#for i in range (0, len(frontTrimProfile1)):
	#	print ('after profile1 front trimming', frontTrimProfile1[i], profile2[i])

	# now trim the back of profile1
	frontAndBackTrimProfile1 = np.trim_zeros(frontTrimProfile1, 'b')

	if frontAndBackTrimProfile1.size == 0:
		return None, None

	# if profile1 was trimmed again, trim profile2 to match
	if len(frontAndBackTrimProfile1) != len(frontTrimProfile1):
		profile2 = profile2[:len(frontAndBackTrimProfile1)]

	if profile2.size == 0:
		return None, None

	#for i in range (0, len(frontAndBackTrimProfile1)):
	#       print ('after profile1 front and back trimming', frontAndBackTrimProfile1[i], profile2[i])

	# now trim based on the resulting profile 2
	frontTrimProfile2 = np.trim_zeros(profile2, 'f')

	if frontTrimProfile2.size == 0:
		return None, None

	if len(frontTrimProfile2) != len(profile2):
		frontAndBackTrimProfile1 = frontAndBackTrimProfile1[-len(frontTrimProfile2):]

	if frontAndBackTrimProfile1.size == 0:
		return None, None

	#for i in range (0, len(frontAndBackTrimProfile1)):
	#       print ('after profile2 front trimming', frontAndBackTrimProfile1[i], frontTrimProfile2[i])

	frontAndBackTrimProfile2 = np.trim_zeros(frontTrimProfile2, 'b')

	if frontAndBackTrimProfile2.size == 0:
		return None, None

	if len(frontAndBackTrimProfile2) != len(frontTrimProfile2):
		frontAndBackTrimProfile1 = frontAndBackTrimProfile1[:len(frontAndBackTrimProfile2)]

	if frontAndBackTrimProfile1.size == 0:
		return None, None

	# make sure the trimmed profiles are the same length as one another - this is a necessary condition for correct trimming
	if len(frontAndBackTrimProfile1) != len(frontAndBackTrimProfile2):
		print ('Trimmed rate profiles must be the same length:\nprofile1 length='+str(len(profile1))+', profile2 length='+str(len(profile2)))
		sys.exit()

	#for j in range (0, len(frontAndBackTrimProfile1)):
	#	print ('Trimmed profiles:', frontAndBackTrimProfile1[j], frontAndBackTrimProfile2[j])

	return frontAndBackTrimProfile1, frontAndBackTrimProfile2

# spline rate profile positions where there is zero data
def splineProfile(profile):

	# profile: input profile to be splined

	# return : the splined profile, where splining is performed over positions where
	#          the reads or rates are zero

	xvals = np.arange(1.0, np.float64(len(profile)+1), 1.0)

	# want to spline over regions with 0.0 rates, but if you use a univariate spline you will always just
	# return zero for these positions - so, make a new array of all non-zero elements, spline that, and then
	# evaluate it for the original range of xvalues
	indices= np.where(profile != 0.0)

	# check to see if we need to spline - there may be no intra-alignment gaps, in which case no spline is required
	if len(indices) == len(profile):
		return profile

	# otherwise, need to spline
	xvalsNonZero   = xvals[indices]
	profileNonZero = profile[indices]

	# produce the spline while keeping all non-zero positions fixed (seems to work but with a minor floating point error)
	spline = scipy.interpolate.InterpolatedUnivariateSpline(xvalsNonZero, profileNonZero)

	# get spline results at all locations, including those with zero reads
	splinedProfile = spline(xvals)

	return splinedProfile

# function to compute a running average
def running_mean(x, N):

	# x      : np.array(); input array
	# N      : int; desired averaging window size

	# returns: result of averaging over x with window size N

	cumsum = np.cumsum(np.insert(x, 0, 0))

	return (cumsum[N:] - cumsum[:-N]) / float(N)

# function to normalize area under curve, mean, and standard deviation of profile
def normProfile(profile):

	# profile: rate profile to be normalized

	# return : normalized rate profile

	newMean = 5.0
	newStd  = 1.0

	# must make sure the area under the curve is the same before applying the EMD, so set the mean and SD and then normalize area
	tempProfile = newMean+((profile-np.mean(profile))*(newStd/np.std(profile, ddof=1)))
	normedProfile = tempProfile/scipy.integrate.simps(tempProfile)

	return normedProfile

# function that actually carries out comparisons between related and unrelated translation rate profiles
def runComparison(alignedProfiles, domains, outdir, n, k, nComp=19, domains2=None):

	# alignedProfiles: dictionary, alignedProfiles[d1+':'+d2] = [processedProfile1, processedProfile2, ]
	# domains: dictionary; keys are 'orfID;resStart-resEnd', values are the corresponding Domain objects; used to select domains for comparisons
	# outdir: str; directory to which output will be written
	# n; str; label for the output file, e.g. rankings1/1_pairwise.txt
	# k: int; determines whether we use to first or second related domain profile for comparisons
	# domains2 (optional): see function arguments for compareDomains
	# nComp (optional, default=19): int; number of unrelated domains to get for each related domain pair

	# returns: nothing, but runs comparisons between domains and writes them to file in outdir

	if outdir[-1] != '/':
		outdir += '/'
	if k not in [0,1]:
		print ('The value of k must be either 0 or 1!')
		return
	output = open(outdir+n+'_pairwise.txt', 'w')
	c1 = 1
	for pair in alignedProfiles:
		compDomains = chooseUnrelatedDomains(pair.split(':')[k], domains, len(alignedProfiles[pair][0]), verbose=False, domains2=domains2)
		c2 = 0
		if compDomains is None:
			print ('Run '+str(n)+': failed to select', nComp, 'unrelated domains for comparison to', pair.split(':')[k])
			continue
		relatedResult = slowMidFast(alignedProfiles[pair][0], alignedProfiles[pair][1])
		output.write(str(c1)+'\t'+str(c2)+'\trelated\t'+pair+'\t'+'%.9f' %relatedResult+'\n')
		for d in compDomains:
			profile1, profile2 = d[1], d[2]
			slowMidFastResult = slowMidFast(normProfile(profile1[:,1]), normProfile(profile2[:,1]))
			c2 += 1
			output.write(str(c1)+'\t'+str(c2)+'\tnot    \t'+pair.split(':')[k]+':'+d[0]+'\t'+'%.9f' %slowMidFastResult+'\n')
		c1 += 1
		# limit the number of pairs run for testing purposes
		#if c1 == 6:
		#	sys.exit()
	output.close()
	rankResults(outdir, outdir+n+'_pairwise.txt')
	return

# function to provide a list of unrelated domains for comparisons
def chooseUnrelatedDomains(domain, domains1, length, N=19, domains2=None, attempts=50000, maxSizeDiff=25, verbose=False):

	# domain: str; identifier of the domain to which comparisons will be made - of the form 'YEL066W;32-173'
	# domains: dictionary; keys are of the form 'YEL066W;32-173', values are corresponding Domain class objects; used to select domains/profiles for comparisons
	#          the random selections will be drawn from this dictionary
	# length: int; length of the aligned profiles for the pair of related domains to which random unrelated domains will be compared
	# N: (optional, default=19) int; the number of unrelated domains to select for domain
	# domains2: (optional, default=None) dictionary; the dictionary of domain objects from which domain is pulled if it differs from domains1
	#           domains2 must be indexible by the keys used for domain
	# attempts: (optional, default=10000) int; the number of times to make random selections before giving up and skipping to the next domain
	# maxSizeDiff: (optional, default=25) int; the maximum size difference to allow between domains
	# verbose: (optional, default=False) Boolean; whether or not to print additional non-essential information to screen and/or to file

	# returns: a list of N tuples, one for each non-paralog Domain object. Format is [selected non-paralog Domain object, alignedProfile1, alignedProfile2, ] in which
	#          alignedProfile1 corresponds to the paralog profile and alignedProfile2 corresponds to the non-paralog profile

	if domains2 is None: # then draw unrelated profiles from the same dictionary as the related profiles
		domains2 = domains1
	unrelated = []
	attempt  = 0
	alreadyChosen = []
	domainList = []
	for key in domains1:
		domainList.append(key)
	while len(alreadyChosen) < N and attempt < attempts:
		choice = domainList[random.randint(0, len(domainList)-1)]
		if choice in alreadyChosen:
			attempt += 1
			continue
		chosenDomain = domains1[choice]
		if chosenDomain.family != domains2[domain].family and chosenDomain.superfamily != domains2[domain].superfamily: # make sure they are not in the same family or superfamily first
			if chosenDomain.SCOPclass != domains2[domain].SCOPclass: # now make sure they are not in the same SCOP class
				if chosenDomain.singleDomainProfile is not None: # make sure the singleDomainProfile exists
					if np.abs((chosenDomain.resEnd-chosenDomain.resStart+1)-(domains2[domain].resEnd-domains2[domain].resStart+1)) <= maxSizeDiff: # make sure the two domains are roughly the same size
						alignedProfile1, alignedProfile2 = alignSingleDomainProfiles(domains2[domain].singleDomainProfile, chosenDomain.singleDomainProfile, length) # align the singleDomainProfile arrays
						if alignedProfile1 is not None and alignedProfile2 is not None: # and make sure the function alignSingleDomainProfiles gave an acceptable result
							alreadyChosen.append(chosenDomain.orfID+';'+str(chosenDomain.resStart)+'-'+str(chosenDomain.resEnd))
							unrelated.append([choice, alignedProfile1, alignedProfile2, ])
							if verbose:
								temp = open('new_output/alignedSingleDomainProfiles/'+domain.split(';')[0]+'-'+domain.split(';')[1]+'_'+\
                                                                             chosenDomain.orfID+'_'+str(chosenDomain.resStart)+'-'+str(chosenDomain.resEnd)+'.txt', 'w')
								p1 = normProfile(alignedProfile1[:,1])
								p2 = normProfile(alignedProfile2[:,1])
								for i in range (0, len(alignedProfile1)):
									temp.write('%.0f' %alignedProfile1[i,0]+'\t'+'%.9f' %p1[i]+'\t'+'%.9f' %p2[i]+'\n')
								temp.close()
		attempt += 1
	if len(unrelated) < 19:
		return None

	return unrelated

# function to take in two singleDomainProfile arrays and then output the aligned region
def alignSingleDomainProfiles(profile1, profile2, length):

	# profile1: np.array; singleDomainProfile 1
	# profile2: np.array; singleDomainProfile 2
	# length  : int; length of the paralogous domain aligned profiles; must make sure all profiles are exactly the same length

	indices1 = [] # indices of array elements from profile1 to keep
	indices2 = [] # indices of array elements from profile2 to keep

	for i in range (0, len(profile1)):
		for j in range (0, len(profile2)):
			if profile1[i,0] == profile2[j,0]:
				indices1.append(i)
				indices2.append(j)

	if len(indices1) != len(indices2):
		print ('Error: length of indices1 list must equal length of indices2 list.')
		return None, None
	if len(indices1) < 50:
		#print ('Insufficient length of region for comparison.')
		return None, None

	newProfile1 = profile1[indices1, :]
	newProfile2 = profile2[indices2, :]

	# if the profiles are too short there is nothing we can do
	if len(newProfile1) < length:
		return None, None
	# if it is too long, use the first length bins so the compared region is exactly the same size between related and unrelated profiles
	elif len(newProfile1) > length:
		return newProfile1[0:length, :], newProfile2[0:length, :]
	else: # if they are exactly the correct length, nothing needs to be done
		return newProfile1, newProfile2

# function to compare rate profiles based on the basis of the 'three bin' method
def slowMidFast(profile1, profile2):

	# profile1: np.array, rate profile of first paralog in pair
	# profile2: np.array, rate profile of second paralog in pair

	# returns : the fraction of positions that have the same low/mid/high classification in both profile1 and profile2

	binnedProfile1 = np.zeros((len(profile1)), dtype=int)
	binnedProfile2 = np.zeros((len(profile2)), dtype=int)
	same = 0 # number of positions at which binnedProfile1 and binnedProfile2 match
	# produce binned profile1
	p1, p2 = np.percentile(profile1, 33.0), np.percentile(profile1, 66.0)
	for i in range (0, len(profile1)):
		if profile1[i] < p1:
			binnedProfile1[i] = 1
		elif profile1[i] >= p1 and profile1[i] < p2:
			binnedProfile1[i] = 2
		elif profile1[i] >= p2:
			binnedProfile1[i] = 3
		else:
			print ('Failed to correctly select percentile bin for profile1.')
			sys.exit()

	# produce binned profile2
	p1, p2 = np.percentile(profile2, 33.0), np.percentile(profile2, 66.0)
	for i in range (0, len(profile2)):
		if profile2[i] < p1:
			binnedProfile2[i] = 1
		elif profile2[i] >= p1 and profile2[i] < p2:
			binnedProfile2[i] = 2
		elif profile2[i] >= p2:
			binnedProfile2[i] = 3
		else:
			print ('Failed to correctly select percentile bin for profile2.')
			sys.exit()

	if len(binnedProfile1) != len(binnedProfile2):
		print ('min/mid/max profiles for paralogs must have the same length;',len(binnedProfile1),'!=',len(binnedProfile2))

	for i in range (0, len(binnedProfile1)):
		if binnedProfile1[i] == binnedProfile2[i]:
			same += 1

	return np.float64(same)/np.float64(len(binnedProfile1))

# function to read in *_pairwise.txt files and rank related and unrelated proteins against one another
def rankResults(outdir, inFile):

	# outdir: str, path to directory in which data will be written
	# inFile: str, path to the *_pairwise.txt file to be ranked

	# returns: the number of times related domains ranked in the each position - printed to screen
	#          and written to file
	if outdir[-1] != '/':
		outdir += '/'
	data = readFile(inFile)
	prefix = inFile.split('/')[-1].split('_')[0]
	dataList = []
	refList = []
	results = np.zeros((20))
	ind = -1
	for i in range (0, len(data), 20):
		dataList = []
		refList = []
		for j in range (i, i+20):
			line = data[j].strip().split()
			dataList.append(line[ind])
			refList.append(line[2])
		dataList = np.array(dataList)
		refList = np.array(refList)
		indices = dataList.argsort()
		sortedResult = refList[indices]
		for j in range (0, len(sortedResult)):
			if sortedResult[j] == 'related':
				results[j] += 1
	ofile = open(outdir+prefix+'_slowMidFast.txt', 'w')
	print (outdir+prefix+'_slowMidFast.txt')
	print ('----------------------------')
	temp = np.flip(results)
	for i in range (0, len(temp)):
		print (str(i+1)+'\t'+'%.0f' %temp[i])
		ofile.write(str(i+1)+'\t'+'%.0f' %temp[i]+'\n')
	ofile.close()
	print ()
	return

def getSummary(outdir, N):

	# outdir: path to directory in which data to be averaged are placed]
	# N: the number of different runs that are to be averaged

	# returns: nothing, but writes a summary file to outdir

	# make a summary file for all of the different comparison runs
	if outdir[-1] != '/':
		outdir += '/'
	for test in ['slowMidFast']: #, 'earthMover', 'earthMover2']:
		results = np.zeros((N, 20))
		for i in range (0, N):
			data = np.loadtxt(outdir+str(i+1)+'_'+test+'.txt', dtype=np.float64)
			for l in range (0, len(data)):
				results[i, l] = data[l,1]
		output = open(outdir+test+'_summary.txt', 'w')
		for k in range (0, 20):
			output.write('%.0f' %(k+1)+'\t'+'%.9f' %np.mean(results[:,k])+'\t'+'%.9f' %np.std(results[:,k], ddof=1)+'\t'+str(len(results[:,k]))+'\t'+'%.0f' %(np.sum(results[0,:]))+'\n')
		output.close()

	return

# function to read a file into memory as a list of lines
def readFile (fileName):

	# fileName: path to a file

	# returns: the contents of the file at fileName as a list of lines or Boolean False if fileName does not exist

	if os.path.exists(fileName):
		fdata = open(fileName)
		data = fdata.readlines()
		fdata.close()
	else:
		print (fileName, 'could not be found.')
		return False

	return data

# function to generate directories
def mkdir (fPath):

	# fPath: dorectory to be created

	# returns: nothing, but creates a new directory

	# N.B., if fPath already exists it will be deleted and remade as a blank directory

	if os.path.exists(fPath):
		subprocess.run('rm -r '+fPath, shell=True)
	subprocess.run('mkdir '+fPath, shell=True)

	return

# function to ask the user a question and take a yes or not answer
def ask_user(question):

	# question: str; the question which the user will be asked on the command line

	# returns: Boolean

	check = str(input('\n'+question+'? (Y/N): ')).lower().strip()
	try:
		if check[0] == 'y':
			return True
		elif check[0] == 'n':
			return False
		else:
			print('Invalid Input')
			return ask_user()
	except Exception as error:
		print("Please enter valid inputs")
		print(error)
		return ask_user()

def permutationTest(file_list1, file_list2, samples=10000, alpha=0.05):

	# file_list1: list;  of file paths as strings corresponding to the N *_slowMidFast.txt files
	# file_list2: list; a second list of file paths as strings corresponding to the N *_slowMidFast.txt files
	# samples: (optional, default = 10000) int; number of random shufflings to perform
	# alpha: (optional, default = 0.05) float; significance level for permutation test

	# N.B., when file_list1 and file_list2 correspond to files generated for related domains and randomly
	#       selected domains this permutation test has Null Hypothesis "related translation rate profiles
        #       are not significantly more similar to one another than randomly selected profiles than are
	#       random pairs of profiles"
	#
	#       By default a two-sided test is performed

	# this assessment is made in a rank-by-rank fashion
	s1_mean = 0.0
	s2_mean = 0.0
	s1_data = []
	s2_data = []
	all_data = []
	T_obs = 0.0
	greater_abs_diff_count = 0

	# loop over ranks
	for rank in range (0, 20):
		# loop over files in file_list1 and file_list2 and collect data
		all_data = []
		s1_data = []
		s2_data = []
		greater_abs_diff_count = 0
		for file in file_list1:
			temp = np.loadtxt(file)[rank, 1]
			s1_data.append(temp)
		s1_mean = np.mean(s1_data)
		for file in file_list2:
			temp = np.loadtxt(file)[rank, 1]
			s2_data.append(temp)
		s2_mean = np.mean(s2_data)
		T_obs = abs(s1_mean - s2_mean)
		all_data = s1_data + s2_data
		print (rank, s1_mean, s2_mean, T_obs, len(s1_data), len(s2_data), len(all_data))

		# get permutations and output sample means and differences
		for i in range (0, samples):
			np.random.shuffle(all_data)
			randSamp1 = all_data[0:20] # first 20 values
			randSamp2 = all_data[-20:] # last 20 values
			mean1 = np.mean(randSamp1)
			mean2 = np.mean(randSamp2)
			#print (i+1, mean1, mean2, abs(mean1-mean2))
			if abs(mean1-mean2) >= T_obs:
				greater_abs_diff_count += 1
		print ('Rank', rank+1, ': permutation absolute difference was greater than |T_obs| in', greater_abs_diff_count, 'samples out of', '%.5e' %samples, 'random shufflings.')
		#sys.exit()
	return

def str2bool(string):
	if string == 'True':
		return True
	else:
		return False

# this script is not currently set up to be run from the command line - call main to perform analysis
if __name__ == '__main__':
	print ('This script is not currently set up to be run from the command line - load it as a module.')
