"
output_codon_mutations.py

Taking sam file as input, finds all mutations compared to reference 
and outputs all codons that contain a mutation (along with the read name)

output:
		readname	codon_pos	ref_codon	alt_codon

NB: this code requires that the reference starts in frame. 

Matt Rich, 05/2017
"""

from argparse import ArgumentParser
import pysam

def output_codon_mutations(sam, offset):
	for s in sam:
		if s.endswith(".bam"):
			samfile = pysam.AlignmentFile(s, "rb")
		elif s.endswith(".sam"):
			samfile = pysam.AlignmentFile(s, "r")
		else:
			raise NameError("{} is neither SAM nor BAM, skipping".format(s))
		
		choose_offset = {0:0, 1:2, 2:1}
		#loop through samfile
		for read in samfile.fetch():	
			if read.is_unmapped == False and read.get_tag("XO") == 0:
				read_offset = choose_offset[read.reference_start%3]
				#use get_aligned_pairs to find mutations
				aligned_read = read.get_aligned_pairs(with_seq = True, matches_only = True)
				for i in range(read_offset, len(aligned_read)-3, 3):
   		            #lowercase bases = mutation in read compared to reference
					#so, if a codon has a lowercase base, it's a mutant codon
					codon = zip(*aligned_read[i:i+3])[2]
					if True in [x.islower() for x in codon]:
						#get codon from read based on aligned_read position
						mutant_codon = read.query_sequence[aligned_read[i][0]:aligned_read[i][0]+3]
						print "\t".join([read.query_name, str(int(aligned_read[i][1]/3.0)),\
										"".join(codon).upper(), mutant_codon])

if __name__ == "__main__":

	parser = ArgumentParser()
	parser.add_argument('sam', action = 'store', nargs="+",
		help = "sam files to be counted (bam also supported)")
	parser.add_argument('--offset', action = 'store', type = int, dest = "offset",
		help = "offset to add to codon number", default=1)

	args = parser.parse_args()
	
	output_codon_mutations(args.sam, args.offset)	

