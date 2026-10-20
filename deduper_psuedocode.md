# Alex Bangs - Reference-Based PCR Duplicate Removal Psuedocode

## The Problem

PCR amplification during library prep makes multiple copies of the same DNA strand, which translates to PCR duplication when multiple of these copies are carried over into sequencing. These duplicates inflate the actual coverage of genetic regions, disrupting analysis of things like gene transcription rates when performing mRNA-seq. We want to remove these duplicates so that we only have one read per original DNA strand.

We can identify these duplicates because they will have the same
- Chromosome (RNAME in SAM header)
- Strand (in FLAG)
- 5'start position (read CIGAR string to accounting for soft-clipping)
- and Unique Molecular Identifier (UMI, in the QNAME)

## Psuedocode

I think it makes sense to use a dictionary to record our duplicates, and add any entries to output file if they aren't found in our dictionary

```
use samtools to sort input file by chromosome and position

set up argparse with arguments that feed in input file, output file, umis.txt

umi_set = args.umi_file

seen_reads = {key = soft-clipping-adjusted 5' position: value = {(strand, umi),}}
    All reads for a chromosome are added, used to determine whether to write read to output file or if it's already been seen.
    If a read has the same 5' position as an entry, it will be excluded from output file if it has the same strand and umi, but will be added (to dictionary value as another set and) and written out if the strand or umi are different.
current_chrom = int 
    Tracks current chromosome. When chromosome changes, seen_reads is reset
chroms_done = {set}
    Tracks chromosomes completed



open input_file for reading, open output_file for writing

use samtools to sort input_file by chromosome and position

For each line in input_file:
    if line starts with "@": (indicates header)
        write to ouput_file
        continue
    
    split line by tabs and assign to list called fields
    qname = fields[0], flag = int(fields[1]), chrom = fields[2], pos = int(fields[3]), cigar_string = fields[5]

    1. skip unmapped sequences (bitwise flag has 0x4 set)
    2. skip umis not in umi_set
    3. if chrom != current_chrom:
        add chrom to chroms_done
        clear seen_reads
        set current_chrom to chrom
    4. evict 
    5. identify soft-clipping-adjusted 5' position for read using pos, strand, and cigar_string
    6. test if current read is a duplicate by comparing adjusted 5' position, strand, and umi
        if not duplicate, add dictionary key:value or expand existing value with other strand/umi. Then write line to output file
```

## Edge Cases

- Unsorted input: detect and exit with error message.
- Chromosome names that are not numeric (`chrX`, scaffolds): treat as strings.
- Reads with `*` CIGAR (unmapped): skipped.
- UMI at wrong position / QNAME with no colon: treat as invalid.