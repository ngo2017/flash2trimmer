# Change log

## flash2trimmer.awk

### ver0.3.0, 2019.05.07
* allow to set old adaptor sequences
* compatible with too many N reads

### ver0.4.0, 2020.03.20
* retain fastq comments

### ver0.5.0, 2021.06.23
* add --any option to fix a bug that prevented trimming adaptors if the adaptor sequence was only found in one of the mates.

### ver0.6.0, 2022.04.12
* add --rrbs option

### ver0.7.0, 2023.10.15
* added --origlen and --mergeseq options

### ver0.8.0, 2024.05.08
* public release
* --any option as default and remove it from the command line option
* change default NUM of -O/--cutadaptO to 1


## chrMRfilter2.awk

### v0.1.0
* initial release combined from bwasplice2_chrMRbwamem.awk and bwasplice2_chrMRfilter.awk

### v0.1.1, 2017.06.29
* supporting an interleaved paired-end fastq from stdin
* add -b/--trim option

### v0.1.2, 2018.01.11
* fixed the usage

### v0.2.0, 2022.04.26
* add -T option

### v0.3.0, 2024.05.09
* public release
* renamed to chrMRfilter2.awk from bwasplice2_chrMRfilter2.awk
* add -w/--bwa option to select BWA command. The default is 'bwa710' if exist, otherwise 'bwa'.


## flash2trimmertk.rb

### v0.3.0, 2023.11.26
* remove -s/--srclen option
* add -z/--stdout option

### v0.4.0, 2023.12.12
* allow to use stdin fastq as '-'
* add -d/--dir option

### v0.4.1, 2024.01.31
* fixed the bug introduced in ver0.3.0 to remove the -s/--srclen option.  
  If an adaptor sequence is found in read1, a longer adaptor sequence found in read2 is not removed.  
  But such reads are very rare, ~7 reads / million trimmed reads.

### v0.4.2, 2024.03.22
* fixed a bug when no trimmed reads are found at the 1st stage (even if trimmed reads is found, their length are shorter than the minimum output length)

### v0.4.3, 2024.04.02
* align read length in pairs when not trimming

### v0.5.0, 2024.04.11
* add -f/--first and --times options
* fixed a bug where Sys.glob in old Ruby did not guarantee ascending order of filenames, so the order of read1 and read2 files could not be predicted.  
  https://bugs.ruby-lang.org/issues/8709  
  https://qiita.com/devzooiiooz/items/43da78f1c3c5c0552710  
  https://stackoverflow.com/questions/6220753/does-dir-glob-guarantee-the-order  
  https://teratail.com/questions/105872  
  Sys.glob in Ruby 3.0 returns sorted filenames by default.  
  'ls' in Shell returns sorted filenames. Try ls -f, which returns unsorted filenames
* add the column 'trimmed' to readlen_aligner.log.txt, which means if the read was trimmed

### v0.6.0, 2024.05.13
* public release
* rename options -5/--trim5 and -3/--trim3 to -1/--read1 and -2/--read2, respectively
* add -t/--threads option
* fixed a bug that prevents processing of directories not under the current directory
* fixed a bug for the format of readlen_aligner.log.txt
* change to use helper scripts in the same directory




