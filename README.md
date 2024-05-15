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
* [pigz](https://zlib.net/pigz/) for flash2trimmertk.rb
* [bwa](https://github.com/lh3/bwa) for flash2trimmertk.rb. I use bwa-0.7.10 as bwa710 for historical reasons.

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
```
mawk -f flash2trimmer.awk -- [options] <fastq1> <fastq2>
mawk -f flash2trimmer.awk -- [options] <interleaved fastq.gz>
  version:0.8.0

Options
  -m/--min LENGTH               Discard reads shorter than LENGTH [30]
  -1/--read1 NUM                Trim NUM bases from 5' end of read1 after adaptor trimming [0]
  -2/--read2 NUM                Trim NUM bases from 5' end of read2 after adaptor trimming [0]
  -a/--adaptors Opt             Cutadapt adaptor sequence option [-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
...
```

## Limitations
* The lengths of read1 and read2 in each pair must be equal
* Not allow piped input

---
## Helper scripts, flash2trimmertk.rb
Compared with flash2trimmer.awk, flash2trimmertk.rb works faster when many reads contain adaptor sequences and is more robust even when one of the reads is of poor quality.
This is achieved by detecting longer adaptor sequences without finding the overlap of paired-end reads before processing them in flash2trimmer.
Optionally, flash2trimmertk.rb discards reads mapped to a decoy database consisting of rDNA etc after trimming.


>With PREFIX as {DIR}/{basename(DIR)}, the source fastq files should be {PREFIX}\_?.fastq.gz or {PREFIX}_?.fq.gz. 
* ChIP-seq  
``` ruby
ruby flash2trimmertk.rb DIR
```

* Stranded RNA-seq  
``` ruby
ruby flash2trimmertk.rb -1 2 --chrMR decoyDB DIR
```

* SMARTer Stranded Total RNA-Seq Kit v2  
``` ruby
ruby flash2trimmertk.rb -1 2 -2 3 --chrMR decoyDB DIR
```

* help
``` ruby
ruby flash2trimmertk.rb -h
```

>The bwa-indexed decoyDB consists of the following sequences
``` 
For mouse
>mtDNA UCSC chrM
>rDNA BK000964.1 TPA_exp: Mus musculus ribosomal DNA, complete repeating unit
>5SrRNA NR_030686.1 Mus musculus 5S RNA (Rn5s), ribosomal RNA
``` 

>Output consists of the following files:
``` 
- {PREFIX}.filtered.fq.gz             Adaptor trimmed reads. Interleaved fastq. Paired-end read order is not preserved
- {PREFIX}.cutadapt.log.txt           The 1st stage cutadapt log
- {PREFIX}.readlen_aligner.log.txt    Table of read length after the 1st stage
- {PREFIX}.flash2trimmer.log.txt      The 2nd stage flash2trimmer log
- {PREFIX}.chrMRfilter2.log.txt       Statics for mapping to the decoy database
``` 

