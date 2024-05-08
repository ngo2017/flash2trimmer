## Introduction

Flash2trimmer trims adapter sequences and library type-specific unwanted sequences from paired-end sequence data. 
It is designed for **haplotype-resolved NGS analysis**, which prefers to trim unwanted sequences completely while keeping the read length as long as possible.

Flash2trimmer finds the overlap of paired-end reads by FLASH2 and adaptor sequences from the non-overlapping outie segments by cutadapt. 
Thus, Flash2trimmer can trim not only very short adaptor sequences, but also unwanted insert sequences and their complementary strands, such as an error-prone stretch caused by random hexamer priming, a template switch oligo, and barcodes.

## Requirements

* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [FLASH2](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [mawk](https://invisible-island.net/mawk/) Recommended, also possible with built-in awk.

## Getting started

* ChIP-seq (NEBNext Ultra II DNA Library Prep Kit)  
No trimming other than adaptor sequences
``` shell
mawk -f flash2trimmer.awk -- read1.fq.gz read2.fq.gz 2> log.txt | pigz > interleaved.fq.gz
```

* Stranded RNA-seq (NEBNext Ultra II Directional RNA Library Prep / Truseq stranded mRNA)  
Trim the first 2 bp of read1 and its complementary strand of read2, an error-prone stretch caused by random hexamer priming during library preparation
``` shell
mawk -f flash2trimmer.awk -- -1 2 read1.fq.gz read2.fq.gz 2> log.txt | pigz > interleaved.fq.gz
```

* SMARTer Stranded Total RNA-Seq Kit v2  
Trim the first 2 bp of read1 (derived from random hexamers), the first 3 bp of read2 (derived from a template switch oligo), and their complementary strands of the reads
``` shell
mawk -f flash2trimmer.awk -- -1 2 -2 3 read1.fq.gz read2.fq.gz 2> log.txt | pigz > interleaved.fq.gz
```

* help
``` shell
mawk -f flash2trimmer.awk -- -h
```

## Limitations
* The lengths of read1 and read2 in each pair must be equal
* Not allow piped input
