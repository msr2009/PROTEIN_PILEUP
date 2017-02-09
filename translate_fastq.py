"""
translate_fastq.py

translates fastq file and outputs either fasta (by default) 
or fastq with sequence qualities either as mean, min, or sum 
of phred scores for codon

Matt Rich, 02/2017
"""

from fastq_tools import read_fastq, print_fastq
from numpy import mean

def main(fastq, offset, qual, nofasta):
		for read in read_fastq(fastq):
			seq = translateSequence(read[1][offset:])
			if qual != None:
				translated_qual = []
				i = 0
				while i <= len(read[2][offset:])-3:
					codon_qual = [ ord(x)-33 for x in read[2][i:i+3] ]
					if qual == "mean":
						translated_qual.append(int(mean(codon_qual)))
					elif qual == "min":
						translated_qual.append(min(codon_qual))
					elif qual == "sum":
						translated_qual.append(sum(codon_qual))
					else:
						raise RuntimeError(qual + \
						" is not a recognized function for quality scores")
						break
					i += 3
				print_fastq((read[0], 
							"".join(seq), 
							"".join([ chr(y+33) for y in translated_qual ])))
			else:
				if not nofasta:
					print ">"+read[0]
				print seq

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
			'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
	try:
		return lookup[codon.lower()]
	except KeyError:
		return "X"

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--fastq', action = 'store', type = str, dest = 'fastq', 
		help = "fastq file to be translated")
	parser.add_argument('--offset', action = 'store', type = int, dest = 'offset',
		help = "sequence offset for translation", default=0)
	parser.add_argument('--quality', action = 'store', type = str, dest = 'quality',
		help = "if quality should be output as FASTQ, \
				how to calculate? (mean, min, sum)", default = None)
	parser.add_argument("--nofasta", action = "store_true", dest = "nofasta",
		help = "do not print FASTA headers (only sequences will be printed). \
				This is overidden by setting '--quality'", default = False)
	args = parser.parse_args()
	
	main(args.fastq, args.offset, args.quality, args.nofasta)	

