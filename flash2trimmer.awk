BEGIN{
  version = "0.8.0";
  program_name = "flash2trimmer";
  OFS="\t"; 
  # FLASH2 options
  # -c, --to-stdout
  # -m, --min-overlap=NUM [10]
  # -M, --max-overlap=NUM [65]
  # -e, --min-overlap-outie=NUM [35]
  # -x, --max-mismatch-density [0.25]

  minlen = 30;
  flash2_x = 0.25; # flash2 default
  flash2_minov = 30;
  cutadapt_minov = 1;
  cmdlog = "/dev/null"; debug_mode = 0; show_origlen = 0; show_merge_seq = 0;
  adaptor_seqopt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
  split("", trim5); trim5[1]=0; trim5[2]=0;
  rrbs = 0;
  flash2_seq = ""; flash2_qual = "";

  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-m" || ARGV[i] == "--min"){
      delete ARGV[i];
      minlen = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-1" || ARGV[i] == "--read1"){
      delete ARGV[i];
      trim5[1] = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-2" || ARGV[i] == "--read2"){
      delete ARGV[i];
      trim5[2] = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-O" || ARGV[i] == "--cutadaptO"){
      delete ARGV[i];
      cutadapt_minov = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="--rrbs"){
      delete ARGV[i];
      rrbs = 1;
    } else if (ARGV[i]=="-o" || ARGV[i] == "--flash2_minoverlap"){
      delete ARGV[i];
      flash2_minov = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-x" || ARGV[i] == "--flash2_x"){
      delete ARGV[i];
      flash2_x = ARGV[++i]+0;
      delete ARGV[i];
    } else if (ARGV[i]=="-a" || ARGV[i] == "--adaptors"){
      delete ARGV[i];
      adaptor_seqopt = ARGV[++i];
      delete ARGV[i];
    } else if (ARGV[i] == "--origlen"){
      show_origlen = 1;
      delete ARGV[i];
    } else if (ARGV[i] == "--mergeseq"){
      show_merge_seq = 1;
      delete ARGV[i];
    } else if (ARGV[i]=="-e" || ARGV[i] == "--stderr"){
      delete ARGV[i];
      cmdlog = "/dev/stderr";
    } else if (ARGV[i]=="-d" || ARGV[i] == "--debug"){
      delete ARGV[i];
      debug_mode = 1;
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }

  blank_seq = ""; blank_qual = "";
  for (j=1; j<=trim5[2]*(1-rrbs)+1; j++){
    blank_seq = blank_seq "N"; 
    blank_qual = blank_qual "!";
  }

  ## prepare source connection
  fastq1 = ARGV[i]; delete ARGV[i];
  fastq2 = ARGV[++i]; delete ARGV[i];
  cutadapt_opt = sprintf("-O %d %s --interleaved --discard-untrimmed --mask-adapter --pair-filter both", cutadapt_minov, adaptor_seqopt);
  flash2_opt = sprintf("-e %d -m %d -x %s -M 10000 -c -t 1", flash2_minov, flash2_minov, flash2_x);

  if (fastq2==""){
    if (fastq1=="")  usage();
    cutadapt_con = sprintf("seqtk seq -q 100 %s | cutadapt %s - 2> %s", fastq1, cutadapt_opt, cmdlog); 
    flash2con = sprintf("flash2 %s --interleaved-input %s 2> %s", flash2_opt, fastq1, cmdlog); 
    srcfastq = sprintf("zcat %s", fastq1);
  } else {
    cutadapt_con = sprintf("seqtk mergepe %s %s | seqtk seq -q 100 - | cutadapt %s - 2> %s", fastq1, fastq2, cutadapt_opt, cmdlog); 
    flash2con = sprintf("flash2 %s %s %s 2> %s", flash2_opt, fastq1, fastq2, cmdlog); 
    srcfastq = sprintf("seqtk mergepe %s %s", fastq1, fastq2);
  }
  # initialize
  space = " "; for (i=1; i<= 8; i++) space = space space;
  flash2_seqid = ""; flash2_seqlen = 0; 
  cutadapt_seqid = ""; split("", cutadapt_len); # 1: read1, 2: read2
  split("", flash2_dist); split("", remain_dist);
  split("", buf); n_pairs = 0; split("", trim3); split("", fq_comment);
  
  while((srcfastq | getline) > 0){
    n_pairs++;
    fq_comment[1] = (NF>1 ? substr($0, length($1)+1) : "");
    sid = parse_seqid($1);
    for (i=2; i<=8; i++){
      srcfastq | getline;
      if (i%2==0) buf[i/2] = $0; # store 1: seq1, 2: qual1, 3: seq2, 4: qual2
      if (i==5) fq_comment[2] = (NF>1 ? substr($0, length($1)+1) : "");
    }
    slen = length(buf[1]);
    if (slen != length(buf[3])) raise(sid ": Sequences are not equal in length");
    flen = flash2_result(sid);
    cflag = cutadapt_result(sid);

    ## ex. trim5[1]=3, trim5[2]=2
    ## 1  @@@------>     l = slen-flen = -2
    ## 2    <-------**   
    ##
    ## 1    @@@------>   l = slen-flen = 2
    ## 2  <-------**   
    l = -1; # adaptor trimmed length
    if (flen > 0 && flen-(trim5[1]+trim5[2]) < slen &&
        (slen-flen < cutadapt_minov || slen-flen==cflag*cutadapt_len[1] || slen-flen==cflag*cutadapt_len[2])){ # merged by FLASH2 & (adaptor too short or adaptor found in either read)
      l = flen;
    } else if (flen == 0 && cflag && cutadapt_len[1]==cutadapt_len[2] && slen-cutadapt_len[2] < flash2_minov){ # not merged by FLASH2 due to too short overlap but have long adaptor sequences
      l = slen - cutadapt_len[1];
    } 
    if (l >= 0){
      remain_dist[l]++;
      l -= trim5[1]+trim5[2]; 
    } else {
      l = slen;
    }

    if (debug_mode){
      if (cflag || flen > 0) debug(flen > 0 ? flen-(trim5[1]+trim5[2]) : l);
    } else if (l >= minlen){
      if (l != slen && show_merge_seq){
        buf[1] = flash2_seq;
        buf[2] = flash2_qual;
        buf[3] = blank_seq; buf[4] = blank_qual;
      }
      print "@" sid "/1" (show_origlen ? " " slen : fq_comment[1]);
      print substr(buf[1], trim5[1]+1, l);
      print "+";
      print substr(buf[2], trim5[1]+1, l);

      print "@" sid "/2" (show_origlen ? " " slen : fq_comment[2]);
      print substr(buf[3], trim5[2]*(1-rrbs)+1, l+trim5[2]*rrbs);
      print "+";
      print substr(buf[4], trim5[2]*(1-rrbs)+1, l+trim5[2]*rrbs);
    }
  }

  ## Log
  print sprintf("# %s, version %s", program_name, version) > "/dev/stderr";
  "flash2 -version" | getline;
  print sprintf("# %s %s", $0, flash2_opt) > "/dev/stderr";
  "cutadapt --version" | getline;
  print sprintf("# cutadapt v%s %s", $0, cutadapt_opt) > "/dev/stderr";
  print sprintf("# pairs: %d", n_pairs) > "/dev/stderr";
  min_l = 10000; max_l = 0;
  for (k in flash2_dist){
    if (k+0 < min_l) min_l = k+0;
    if (k+0 > max_l) max_l = k+0;
  }
  for (k in remain_dist){
    if (k+0 < min_l) min_l = k+0;
    if (k+0 > max_l) max_l = k+0;
  }
  print "length", "with_adaptor", "flash2_merged" > "/dev/stderr";
  for (i = min_l; i <= max_l; i++) print i, remain_dist[i]+0, flash2_dist[i]+0 > "/dev/stderr";
}

## grobal_variable: flash2dist, flash2_seqlen, flash2_seqid, flash2_seq, flash2_qual
function flash2_result(seqid){
  if (flash2con == "") return(0);
  if (flash2_seqid==""){
    if ((flash2con | getline) > 0){
      flash2_seqid = parse_seqid($1);
      flash2con | getline flash2_seq; flash2_seqlen = length(flash2_seq);
      flash2con | getline; flash2con | getline flash2_qual;  
      flash2_dist[flash2_seqlen]++;
    } else {
      close(flash2con);
      flash2con = ""; flash2_seqid = "";
    }
  }
  if (flash2_seqid == seqid){
    flash2_seqid = "";
    return(flash2_seqlen);
  } else return(0);
}

function cutadapt_result(seqid,  i){
  if (cutadapt_con == "") return(0);
  if (cutadapt_seqid == ""){
    if ((cutadapt_con | getline) > 0){
      cutadapt_seqid = parse_seqid($1);
      for (i=2; i<=8; i++){
        cutadapt_con | getline; 
        if (i%4==2) cutadapt_len[(i+2)/4] = length($1)-match($1, "N+$")+1;
      }
    } else {
      close(cutadapt_con);
      cutadapt_con = ""; cutadapt_seqid = "";
    }
  }
  if (cutadapt_seqid == seqid){
    cutadapt_seqid = "";
    return(1);
  } else return(0);
}
    

## remove /1 or /2 suffix
function parse_seqid(str,   i){
  i = index(str, "/");
  return(substr(str, 2, (i > 0 ? i-1 : length(str))-1));
}

function debug(ol,    seq1, seq2, l){
  seq1 = substr(buf[1], 1, trim5[1]) " " substr(buf[1], trim5[1]+1, ol) " " substr(buf[1], trim5[1]+ol+1);
  seq2 = substr(buf[3], 1, trim5[2]) " " substr(buf[3], trim5[2]+1, ol) " " substr(buf[3], trim5[2]+ol+1);
  seq2 = complement_seq(seq2);
  l = trim5[1]-(slen -ol -trim5[2]);
  print ">" sid;
  print sprintf("# overlapLen: %d, FLAHS2_trimlen: %d, cutadapt_trimlen1: %d, cutadapt_trimlen2: %d", ol, (flen > 0)*(slen-flen), cflag*cutadapt_len[1], cflag*cutadapt_len[2]);
  print substr(space, 1, l < 0 ? -l : 0) seq1;
  print substr(space, 1, l > 0 ? l : 0)  seq2;
  print "";
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
  print sprintf("mawk -f %s.awk -- [options] <fastq1> <fastq2>", program_name) > "/dev/stderr";
  print sprintf("mawk -f %s.awk -- [options] <interleaved fastq.gz>", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print "Options"  > "/dev/stderr";
  print sprintf("  -m/--min LENGTH               Discard reads shorter than LENGTH [%s]", minlen) > "/dev/stderr";
  print sprintf("  -1/--read1 NUM                Trim NUM bases from 5' end of read1 after adaptor trimming [%s]", trim5[1]) > "/dev/stderr";
  print sprintf("  -2/--read2 NUM                Trim NUM bases from 5' end of read2 after adaptor trimming [%s]", trim5[2]) > "/dev/stderr";
  print sprintf("  -a/--adaptors Opt             Cutadapt adaptor sequence option [%s]", adaptor_seqopt) > "/dev/stderr";
  print         > "/dev/stderr";
  print         "  --rrbs                        The -2 option applies to only read1" > "/dev/stderr";
  print sprintf("  -O/--cutadaptO NUM            Cutadapt -O option [%s]", cutadapt_minov) > "/dev/stderr";
  print sprintf("  -o/--flash2_minoverlap NUM    FLASH2 -e & -m option [%s]", flash2_minov ) > "/dev/stderr";
  print sprintf("  -x/--flash2_x FLOAT           FLASH2 -x option [%s]", flash2_x) > "/dev/stderr";
  print         "  --origlen                     Show the length of the original sequence as a fastq comment" > "/dev/stderr";
  print         "  --mergeseq                    Show the FLASH2 merged sequence as read1 if trimmed. Read2 is set to 'N'" > "/dev/stderr";
  print         "  -d/--debug                    Debug mode" > "/dev/stderr";
  print         "  -e/--stderr                   Output FLAHS2/cutadapt log to stderr instead of discarding" > "/dev/stderr";
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}
