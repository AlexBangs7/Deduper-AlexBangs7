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

Shorthand used below: **A** = `AACGCCAT`, **B** = `AAGGTACC`, **C** = `GAACAGGT`.
QNAME prefix for all reads: `NS500451:154:HWKTMBGXX:1:11101:`

### 3.1 `unit_tests/input.sam` (16 reads, coordinate sorted)

Columns are tab-separated.

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:2	LN:1000
@SQ	SN:3	LN:1000
NS500451:154:HWKTMBGXX:1:11101:1000:1001:AACGCCAT	0	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1002:AACGCCAT	0	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1003:AAGGTACC	0	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1004:AACGCCAT	16	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1005:AACGCCAT	16	2	100	255	8M2S	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1006:AACGCCAT	16	2	100	255	8M	*	0	0	ACGTACGT	IIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1007:AACGCCAT	0	2	103	255	3S7M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1008:AACGCCAT	0	2	105	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1009:AACGCCAA	0	2	120	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1010:GTNACCTA	0	2	130	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1011:AACGCCAT	0	3	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1012:AACGCCAT	0	3	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1013:AAGGTACC	0	3	200	255	5M100N5M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1014:AAGGTACC	0	3	200	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1015:GAACAGGT	16	3	300	255	5M100N5M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1016:GAACAGGT	16	3	400	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
```

### 3.2 `unit_tests/expected_output.sam` (8 reads kept)

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:2	LN:1000
@SQ	SN:3	LN:1000
NS500451:154:HWKTMBGXX:1:11101:1000:1001:AACGCCAT	0	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1003:AAGGTACC	0	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1004:AACGCCAT	16	2	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1006:AACGCCAT	16	2	100	255	8M	*	0	0	ACGTACGT	IIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1008:AACGCCAT	0	2	105	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1011:AACGCCAT	0	3	100	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1013:AAGGTACC	0	3	200	255	5M100N5M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
NS500451:154:HWKTMBGXX:1:11101:1000:1015:GAACAGGT	16	3	300	255	5M100N5M	*	0	0	ACGTACGTAC	IIIIIIIIII	NH:i:1
```

### 3.3 Test case descriptions (content for `unit_tests/README.md`)

| QNAME suffix | Chr | Strand | POS | CIGAR | UMI | Adjusted 5' | Result | Reason |
|---|---|---|---|---|---|---|---|---|
| `1001` | 2 | + | 100 | 10M | A | 100 | **Keep** | First read at this locus |
| `1002` | 2 | + | 100 | 10M | A | 100 | Remove | Exact duplicate of 1001 |
| `1003` | 2 | + | 100 | 10M | B | 100 | **Keep** | Same position/strand, **different UMI** = different molecule |
| `1004` | 2 | - | 100 | 10M | A | 109 | **Keep** | Same POS and UMI as 1001, but **opposite strand** |
| `1005` | 2 | - | 100 | 8M2S | A | 109 | Remove | Reverse strand; right soft-clip makes 5' = 100+8-1+2 = 109, same as 1004 |
| `1006` | 2 | - | 100 | 8M | A | 107 | **Keep** | Same POS/strand/UMI as 1004, but **different 5' end** (shorter read) |
| `1007` | 2 | + | 103 | 3S7M | A | 100 | Remove | Forward; left soft-clip makes 5' = 103-3 = 100, same as 1001 |
| `1008` | 2 | + | 105 | 10M | A | 105 | **Keep** | Same UMI/strand, **different position** |
| `1009` | 2 | + | 120 | 10M | `AACGCCAA` | n/a | Discard | UMI not in list (1 base error) |
| `1010` | 2 | + | 130 | 10M | `GTNACCTA` | n/a | Discard | UMI contains `N`, not in list |
| `1011` | 3 | + | 100 | 10M | A | 100 | **Keep** | Same POS/UMI/strand as 1001 but **different chromosome** |
| `1012` | 3 | + | 100 | 10M | A | 100 | Remove | Duplicate of 1011 |
| `1013` | 3 | + | 200 | 5M100N5M | B | 200 | **Keep** | Spliced read; forward 5' unaffected by `N` |
| `1014` | 3 | + | 200 | 10M | B | 200 | Remove | Duplicate of 1013 (different CIGAR, same 5' start) |
| `1015` | 3 | - | 300 | 5M100N5M | C | 409 | **Keep** | Reverse spliced: 300+5+100+5-1 = 409 (`N` consumes reference) |
| `1016` | 3 | - | 400 | 10M | C | 409 | Remove | Different POS (400), but same 5' end (409) as 1015: duplicate. Also tests the eviction window |



## Edge Cases

- Unsorted input: detect and exit with a clear error message.
- Chromosome names that are not numeric (`chrX`, scaffolds): treat as strings.
- Reads with `*` CIGAR (unmapped): skipped.
- UMI at wrong position / QNAME with no colon: treat as invalid.
- Blank last line / Windows line endings: strip whitespace before parsing.
- Empty input (header only): output is just the header.
- 