#!/usr/bin/env python


import argparse

# Inputs: SAM file, output file path, UMI file
# Outputs: Deduped SAM file, plus a summary of the following:
# Number of header lines and unique reads in output, count of wrong UMIs, count of removed duplicates
# Number of reads per chromosome


# Following peer review guidance of breaking code up into functions. Taking inspiration on functions from psuedocode of lyxxng

def get_args():
	'''Function creates and returns argparse object to take in 
	an input file, an output file, a umi list file, and help options'''
	parser=argparse.ArgumentParser()
	parser.add_argument("-f","--file", help="required filepath leading to a sorted SAM file", required=True)
	parser.add_argument("-o","--outfile", help="required filepath for desired location of output SAM file", required=True)
	parser.add_argument("-u","--umi", help="required filepath leading to a file with a list of UMI", required=True)
	return parser.parse_args()
args = get_args()


def get_umi_list(umi_file:str):
	'''Function that takes the UMI file path as input 
	and returns a set containing all unique UMIs'''
	UMIs = []

	with open(umi_file, "r") as umi_doc:
		for line in umi_doc:
			line = line.strip("\n")
			UMIs.append(line)
	return UMIs 				# later use UMI_list = get_umi_list(args.umi)


def read_alignment(alignment:str):
	'''Function that, when given a read, pulls out the chrom (RNAME), 
	strand (from bitwise flag), 1-base left-most position (POS), and CIGAR string'''
	read = alignment.strip("\n").split()
	UMI = (read[0].split(":"))[-1]
	FLAG = int(read[1])
	RNAME = read[2] # chromosome
	POS = int(read[3])
	CIGAR = read[5]
	return UMI, FLAG, RNAME, POS, CIGAR


def parse_CIGAR_string(CIGAR_str:str):
	''' Function that takes the CIGAR string and converts
	it to a list that splits up the numbers and characters'''
	# Input: CIGAR string (e.g., "6S64M")
	# Output: Parsed CIGAR string as list (e.g., ["6", "S", "64", "M"])

	temp_C_str = "" 													# temporary str used to combining adjacent digits (so CIGAR list reads 64M as [64,M] and not [6,4,M])
	C_list = [] 														# CIGAR list containing elements of original CIGAR string, split by letter and number
	
	for chr in CIGAR_str:
		if chr.isnumeric():
			temp_C_str = temp_C_str + chr 								# add digit to temporary string
		elif chr.isalpha() or chr == "=":
			C_list.append(temp_C_str)									# append existing temporary string (i.e. number) to output list
			temp_C_str = chr											# reset temporary string to now be equal to current letter
			C_list.append(temp_C_str)									# append temporary string (i.e. letter) to output list
			temp_C_str = ""												# clear temporary string
		else:
			raise Exception("Input is not an acceptable CIGAR string")
	return C_list

def get_five_prime_start(strand:str, pos:int, CIGAR:list):
	'''Function that takes the bitwise flag (or strand), 1-based leftmost mapping position 
	(first matching base position), and parsed CIGAR list, and calculates the 5' start of read'''	
	left_clip = 0
	if CIGAR[1] in ("S", "H"):
		left_clip = int(CIGAR[0]) 					# left_clip != 0 only when CIGAR has soft-clipping, so we can subtract from all conditions whether or not soft-clipping is present
	if strand == "Reverse": 							# if sequence is on reverse strand
		read_size = 0
		for i in range(0, len(CIGAR), 2):
			if CIGAR[i + 1] != "I":
				read_size += int(CIGAR[i])
		start_pos = pos + read_size - left_clip 
	else: 											# if sequence on forward strand
		start_pos = pos - left_clip 				# if soft-clipped, will subtract amount from pos. If not clipped, will subtract 0
	return start_pos

def deduping(input_file, output_file, umi_set):
	''' deduper function, parses each line from input SAM file
    and writes unique sequences to output SAM file'''
    # Inputs: Input and output SAM file paths, UMI set
    # Output: Output file with unique SEQs
	
	input, output = open(input_file, "r"), open(output_file,"w")
	
	current_chrom = 0
	seen_reads = set() # set, each entry is umi, chromosome, strand, and start site

	# Deduping analysis data
	header_lines = 0
	uniq_reads = 0
	uniq_reads_check = 0
	wrong_umis = 0
	removed_duplicates = 0

	for line in input:
		if line.startswith("@"): # indicates header lines
			output.write(line)
			header_lines += 1
			continue
		UMI, FLAG, chrom, POS, CIGAR = read_alignment(line)
		if ((FLAG & 4) == 4): # if sequence is unmapped
			continue
		if UMI not in umi_set:
			wrong_umis += 1
			continue
		if chrom != current_chrom:
			print(f"{current_chrom}\t{len(seen_reads)}")
			uniq_reads_check += len(seen_reads)
			
			current_chrom = chrom
			seen_reads = set()

		strand = "Reverse" if ((FLAG & 16) == 16) else "Forward"
		CIGAR_list = parse_CIGAR_string(CIGAR)
		start_site = get_five_prime_start(strand, POS, CIGAR_list)

		current_sequence = (UMI, chrom, strand, start_site)
		if current_sequence not in seen_reads:
			seen_reads.add(current_sequence)
			output.write(line)
			uniq_reads += 1
		else:
			removed_duplicates +=1
	
	print(f"{current_chrom}\t{len(seen_reads)}")
	uniq_reads_check += len(seen_reads)

	input.close()
	output.close()
	print(header_lines, uniq_reads,wrong_umis,removed_duplicates)
	print("Deduping completed :)")



def main():	
	''' Get UMI list then run deduping function '''
	
	umi_set = get_umi_list(args.umi)
	deduping(args.file, args.outfile, umi_set)
	return

if __name__ == "__main__":
	main()

''''
For each line in input_file:
	# 1. skip unmapped sequences (bitwise flag has 0x4 set)
	# 2. skip umis not in umi_set
	3. if chrom != current_chrom:
		add chrom to chroms_done
		clear seen_reads
		set current_chrom to chrom
	4. evict 
	5. identify soft-clipping-adjusted 5' position for read using pos, strand, and cigar_string
	6. test if current read is a duplicate by comparing adjusted 5' position, strand, and umi
		if not duplicate, add dictionary key:value or expand existing value with other strand/umi. Then write line to output file
'''
