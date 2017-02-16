"""
sam_to_fullFQ.py

This script will take a sam file with either paired or unpaired reads,
and a wildtype sequence, outputs a fastq file containing the entire sequence,
with relevant mutations from each read (or read pair)

Matt Rich, 02/2017
"""

import pysam
from itertools import groupby


def main(sam, wildtype, paired, excludeWT):
	#read wildtype FASTA file
	for x in fasta_iter(wildtype):
		wt_name, wt_seq = x

	#open sam/bam file
	if sam.endswith(".bam"):
		samfile = pysam.AlignmentFile(sam, "rb")
	elif sam.endswith(".sam"):
		samfile = pysam.AlignmentFile(sam, "r")
	else:
		raise NameError("file is neither SAM nor BAM")

	#loop through samfile
	for read in samfile.fetch():
		#print first line of FASTQ record
		if not paired:
			#skip if read has indels
			if read.get_tag("XO") == 0:
				new_read = wt_seq[0:read.reference_start] + read.query_alignment_sequence + wt_seq[read.reference_end:]
				new_qual = "".join(["A" for i in wt_seq[0:read.reference_start]]) + \
							"".join([chr(y+33) for y in list(read.query_alignment_qualities)]) + \
							"".join(["A" for j in wt_seq[read.reference_end:]])
				
				if excludeWT == False or new_read != wt_seq:
					printFASTQ(read.query_name, new_read, new_qual)
		
		else:
			read2 = samfile.next()
			#skip if read has indels
			if read.get_tag("XO") == 0 and read2.get_tag("XO") == 0:
				#if read1 is to the left of read2:
				if read.reference_start < read2.reference_start:
					new_read, new_qual = merge_reads(read, read2, wt_seq)
				else:
					new_read, new_qual = merge_reads(read2, read, wt_seq)
				
				if excludeWT == False or new_read != wt_seq:
					printFASTQ(read.query_name, new_read, new_qual)

def merge_reads(r1, r2, wt):
		wt_q = 63
		#do reads overlap, but the overlap was too short for PEAR to deal with?
		overlap = r1.reference_end - r2.reference_start
		if overlap < 0:
			overlap = 0

		new_read = wt[0:r1.reference_start] + \
			r1.query_alignment_sequence + \
			wt[r1.reference_end:r2.reference_start] + \
			r2.query_alignment_sequence[overlap:] + \
			wt[r2.reference_end:]
		
		new_qual = "".join([chr(wt_q) for i in wt[0:r1.reference_start]]) + \
			"".join([chr(y+33) for y in list(r1.query_alignment_qualities)]) + \
			"".join([chr(wt_q) for j in wt[r1.reference_end:r2.reference_start]]) + \
			"".join([chr(y+33) for y in list(r2.query_alignment_qualities)[overlap:]]) + \
			"".join([chr(wt_q) for j in wt[r2.reference_end:]])
		
		return (new_read, new_qual)

def printFASTQ(name, seq, qual):
		print "@"+name
		print seq
		print "+"
		print qual


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
	parser.add_argument('--paired', action = 'store_true', dest = 'paired', 
		default = False, help = "sam file contains paired reads")
	parser.add_argument('--wt', action = 'store', type = str, dest = 'wildtype', 
		help = "fasta file containing wildtype sequence")
	parser.add_argument('--exclude-wt', action = 'store_true', dest = 'excludeWT', 
		help = "only print mutated reads", default = False)
	args = parser.parse_args()
	
	main(args.sam, args.wildtype, args.paired, args.excludeWT)	

