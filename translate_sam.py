"""
translate_sam.py

This script translates a SAM (or BAM file) in order to create a pileup of mutations
at each position in the protein sequence.

This script disregards whether reads are paired or unpaired and counts variants at
every position. Users should merge paired reads before using this script (with
something like PEAR), and then run separately on the merged and unmerged alignments.

Matt Rich, 05/2017
"""

import pandas as pd
import pysam
from itertools import product, groupby

def main(sam, wildtype, offset, pileup_type, include_flank):
	#read wildtype FASTA file
	for x in fasta_iter(wildtype):
		wt_name, wt_seq = x

	#make dictionary to store pileup information
	pileup_types = {"aa" : ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
			   				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'],
					"codon" : ["".join(x) for x in product(["A","C","G","T"], repeat=3)]}

	#count by codon rather than position
	pileup = { i/3:{} for i in range(0,len(wt_seq),3) }
	for pos in pileup:
		pileup[pos] = {x:0 for x in pileup_types[pileup_type]}

	#this is so stupid, but I can't figure out the arithmetic to get what I want for
	#for the if not statement below
	last_lookup = {0:0, 1:2, 2:1}

	#open sam/bam file
	if sam.endswith(".bam"):
		samfile = pysam.AlignmentFile(sam, "rb")
	elif sam.endswith(".sam"):
		samfile = pysam.AlignmentFile(sam, "r")
	else:
		raise NameError("file is neither SAM nor BAM")

	#loop through samfile
	for read in samfile.fetch():
		#skip if read has indels or Ns
		if read.get_tag("XO") == 0 and "N" not in read.query_alignment_sequence:
			"""
			1) determine how many wildtype bases we need to add to each end of the
			aligned segment
			2) create new sequence by appending those bases on either end
			3) loop through sequence, adding to pileup
			"""

			new_seq = ""
			first_pos = 0
			if include_flank:
				first_pos = read.reference_start-(read.reference_start+offset)%3
				first = wt_seq[first_pos : read.reference_start]		
				last = wt_seq[read.reference_end : read.reference_end+last_lookup[(read.reference_end+offset)%3]]	
		 		new_seq = first + read.query_alignment_sequence + last
			
			else:
				#if we're not appending wildtype sequences then we just take the
				#in-frame part of the query_alignment_sequence
				
				first_pos = read.reference_start + (3-(read.reference_start+offset)%3)%3
				new_seq = read.query_alignment_sequence[(3-(read.reference_start+offset)%3)%3:\
														-1*((read.reference_end+offset)%3)]			
	
			#weird stuff can happen if the reads start at the very beginning of the
			#sequence, so skip anything where first_pos < 0
			if read.reference_start-offset >= 0:
				#add codons or translate, then add to pileup
				for c in range(0,len(new_seq),3):
					if pileup_type == "codon":
						pileup[(first_pos+c)/3][new_seq[c:c+3]] += 1
					elif pileup_type == "aa":
						pileup[(first_pos+c)/3][translateSequence(new_seq[c:c+3])] += 1
			

	#make a pandas dataframe (for easy printing) from pileup dict
	pd_pileup = pd.DataFrame.from_dict(pileup, 'index')
	if pileup_type == "codon":
		pd_pileup = pd_pileup.reindex_axis(sorted(pd_pileup.columns), axis=1)
	elif pileup_type == "aa":
		pd_pileup = pd_pileup[pileup_types["aa"]]	

	print pd_pileup.to_csv(sep="\t", index_label="codon")

def translateSequence(seq):
	translated_seq = ""
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
				'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
				'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
				'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
				'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
				'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
				'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
				'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
				'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
				'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
				'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
				'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
				'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
				'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
				'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
				'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	return lookup[codon.lower()]

def fasta_iter(fasta_name, noiter=False):
		"""
		given a fasta file. yield tuples of header, sequence
		this by brentp from https://www.biostars.org/p/710/
		"""
		fh = open(fasta_name)
		# ditch the boolean (x[0]) and just keep the header or sequence since
		# we know they alternate.
		faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
		for header in faiter:
				# drop the ">"
				header = header.next()[1:].strip()
				# join all sequence lines to one.
				seq = "".join(s.strip() for s in faiter.next())
				yield header, seq

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--sam', action = 'store', type = str, dest = 'sam', 
		help = "sam or bam file")
	parser.add_argument('--wt', action = 'store', type = str, dest = 'wildtype', 
		help = "fasta file containing wildtype sequence")
	parser.add_argument('--offset', action = 'store', type = int, dest = 'offset',
		help = "offset for reading frame of seq, either 0, 1, or 2 (default=0)", default=0)
	parser.add_argument('--pileup-type', action = 'store', type = str, dest = 'pileup',
		help = "pileup as codon or aa", default="codon")
	parser.add_argument('-q', '--quality', action='store', type = int, dest='qual',
		help = "minimum PHRED quality of mutations, default=20", default = 20)
	parser.add_argument('--include-flank', action = 'store_true', dest = 'include',
		help = "pad reads with flanking wildtype sequence. By default, will shorten\
				reads to first and last full codons.", default=False)
	args = parser.parse_args()
	
	main(args.sam, args.wildtype, args.offset, args.pileup, args.include)	

