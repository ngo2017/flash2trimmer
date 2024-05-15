BEGIN {
  version = "0.3.0";
  program_name = "chrMRfilter2";
  FS="\t";
  coverage = 0.75; db = ""; paired = 0; n_thread = 4;
  stderr = 0; trim1 = 0; min_score = 0;
  bwa_cmd = "";

  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-d" || ARGV[i] == "--db"){
      if (system("test -r " ARGV[i+1])) raise(ARGV[i+1] " is not found");
      delete ARGV[i];      
      db = ARGV[++i];
      delete ARGV[i];
    } else if (ARGV[i]=="-w" || ARGV[i] == "--bwa"){
      delete ARGV[i];
      bwa_cmd = ARGV[++i];
      delete ARGV[i];
    } else if (ARGV[i]=="-c" || ARGV[i] == "--coverage"){
      delete ARGV[i];
      coverage = ARGV[++i];
      delete ARGV[i];
    } else if (ARGV[i]=="-p" || ARGV[i] == "--paired"){
      paired = 1;
      delete ARGV[i];
    } else if (ARGV[i]=="-b" || ARGV[i] == "--trim"){
      delete ARGV[i];
      trim1 = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-T"){
      delete ARGV[i];
      min_score = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-t" || ARGV[i] == "--thread"){
      delete ARGV[i];      
      n_thread = ARGV[++i];
      delete ARGV[i];
    } else if (ARGV[i]=="-e" || ARGV[i] == "--error"){
      stderr = 1;
      delete ARGV[i];
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }
  if (db == "") raise("No database");
  if (bwa_cmd == ""){
    if (("which bwa710 2> /dev/null" | getline) > 0){
      bwa_cmd = "bwa710";
    } else if (("which bwa 2> /dev/null" | getline) > 0){
      bwa_cmd = "bwa";
    } else {
      raise("BWA is not found");
    }
  }
  split("", report); qname = "";  
  mempipe = sprintf("%s mem -t %s %s %s %s %s %s %s", 
                    bwa_cmd, 
                    n_thread,
                    paired ? "-P -S" (ARGV[i]=="" ? " -p" : "") : "", 
                    min_score > 0 ? ("-T " min_score) : "",
                    db, 
                    ARGV[i]=="" ? "-" : ARGV[i],
                    paired ? ARGV[i+1] : "",
                    stderr ? "" : "2> /dev/null");
		    
  while((mempipe | getline) > 0){
    if ($0 ~ /^@/ || int($2/256) > 0){ continue; } # skip header or supplementaly/non-primary alignment
    if (int($2/4)%2 == 0 && $14 !~ /^AS/) raise("AS is not found " $0); # mapped but no AS
    ## keep rname as a filtering flag if unmapped or AS/seqLen > coverage
    rname = ( int($2/4)%2 || substr($14, 6)/length($10) <= coverage ? "" : $3);  

    if ($2%2 == 0){ # single end
      if (rname=="") outseq($1, $10, $11, int($2/16)%2, trim1);
    } else {
      if (int($2/64)%2){ # paired-end read1
        if (qname != "") raise("Not paired " $0);
        qname = $1; seq1 = $10; qual1 = $11; rev1 = int($2/16)%2; rname1 = rname;
        continue;
      }
      if (int($2/128)%2){ # paired-end read2
        if (qname != $1) raise("Not paired " $0);
        qname = "";
        if (rname1 != "") rname = rname1;  # RNAME mapped by read1 is prefered
        if (rname == ""){
          outseq($1 "/1", seq1, qual1, rev1, trim1);
          outseq($1 "/2", $10, $11, int($2/16)%2);
        } 
      }
    }
    report[rname]++;
  }
}

function outseq(qname, seq, qual, rev, trim){
  if (rev){
    seq = complement_seq(seq);
    qual = reverse(qual);
  }
  if (trim > 0){
    seq = substr(seq, trim+1);
    qual = substr(qual, trim+1);
  }
  print "@" qname;
  print seq;
  print "+";
  print qual;
}
function reverse(str,  new_str, i){
  new_str = "";
  for (i=length(str); i>0; i--)
    new_str = new_str substr(str, i, 1);
  return(new_str);
}
function complement_seq(seq){
  seq = reverse(seq);
  gsub(/A/, "_", seq); gsub(/T/, "A", seq); gsub(/_/, "T", seq);
  gsub(/C/, "_", seq); gsub(/G/, "C", seq); gsub(/_/, "G", seq);
  return(seq);
}

function usage(){
  print sprintf("mawk -f %s.awk -- [options] [fastq1.gz] [fastq2.gz]", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print "Filtering reads mapped to the decoy DB" > "/dev/stderr";
  print > "/dev/stderr";
  print "options" > "/dev/stderr";
  print         "  -d/--db Database :     BWA indexed decoy database. Required" > "/dev/stderr";  
  print         "  -p/--paired :          Paired end mode." > "/dev/stderr";  
  print         "  -w/--bwa Command :     BWA command ['bwa710' if exist, otherwise 'bwa']" > "/dev/stderr";  
  print         "  -T Score :             BWA MEM minimum score option [default]" > "/dev/stderr";  
  print sprintf("  -t/--thread NUM :      Number of BWA MEM threads [%d]", n_thread) > "/dev/stderr";  
  print > "/dev/stderr";
  print sprintf("  -c/--coverage Float :  Coverage threshold. Discard reads with [the alignment score, AS]/[read length] > coverage [%s]", coverage) > "/dev/stderr";  
  print sprintf("  -b/--trim NUM :        Trim NUM bp from 5' end of read1 after filtering [%d]", trim1) > "/dev/stderr";  
  print         "  -e/--error :           Output BWA MEM stderr" > "/dev/stderr";  
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}

END {
  close(mempipe);
  n_filtered = 0;
  for (key in report){
    if (key != "") n_filtered += report[key];
    else n_left = report[key];
  }
  
  if (n_filtered + n_left > 0){
    print sprintf("# %s,  version %s", program_name, version) > "/dev/stderr";
    for (key in report)
      if (key != "") print sprintf("# %s: %d", key, report[key]) > "/dev/stderr";
    print sprintf("# fastq records: %d", n_left+n_filtered) > "/dev/stderr";
    print sprintf("# Num. reads left: %d", n_left) > "/dev/stderr"
  }

  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
}
