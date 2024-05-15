require 'optparse'
require 'tmpdir'
require 'fileutils'

class Tmpfile
  @index = 0
  @tmpdir = nil
  class << self
    def create(suffix=''); "#{tmpdir}/#{@index += 1}_#{suffix}"; end
    def tmpdir; @tmpdir ||= Dir.mktmpdir; end
    def clear;
      FileUtils.remove_entry_secure(@tmpdir) if @tmpdir
      @tmpdir = nil
    end
  end
end

at_exit do
  Tmpfile.clear
end

Version = '0.6.0'
params = {
  :minlen => 30,
  :olen => 8,
  :trim5 => 0,
  :trim3 => 0,
  :trimA => 0,
  :adapt_opt => nil,
  :chrMR => nil,
  :stdout => false,
  :dry => false,
  :dir => nil,
  :adapt_1st => '-a AGATCGGAAGAGC -A AGATCGGAAGAGC',
  :cutadapt_times => '',
  :chrMR_minscore => nil,
  :threads => 8,
}
opts = OptionParser.new(<<HEADER){ |opt|
ruby #{File.basename($0)} [options] <DIR>
cat interleaved.fastq | ruby #{File.basename($0)} [options] -d <DIR> -

  For paired-end reads only.
  At the 1st stage, run cutadapt {-f/--first ADAPTOROPT} -O {-o/--overlap OLEN}
  If adapter sequences are found, additional trimming is performed according to the options -1/--read1, -2/--read2, and -A/--trimA.
  Otherwise, the length of the reads will be aligned to the shorter one if necessary, and then perform the 2nd stage trimming using flash2trimer with the options -1/--read1 and -2/--read2.

  With PREFIX as {DIR}/{basename(DIR)}, the source fastq files should be {PREFIX}_? .fastq.gz or {PREFIX}_? .fq.gz

  Output consists of the following files:
   - {PREFIX}.filtered.fq.gz             Adaptor trimmed reads. Interleaved fastq. Paired-end read order is not preserved
   - {PREFIX}.cutadapt.log.txt           The 1st stage cutadapt log
   - {PREFIX}.readlen_aligner.log.txt    Table of read length after the 1st stage
   - {PREFIX}.flash2trimmer.log.txt      The 2nd stage flash2trimmer log
   - {PREFIX}.chrMRfilter2.log.txt       Statics for mapping to the decoy database

HEADER
  # opt.on('-s', '--src LEN', Integer, "Source read length. Required"){ |v| params[:srclen] = v }
  opt.on('-m', '--min LEN', Integer, "Discard reads shorter than LEN [#{params[:minlen]}]"){ |v| params[:minlen] = v }
  opt.on('-o', '--overlap OLEN', Integer, "Required length of overlap between read and adapter at the 1st stage [#{params[:olen]}]"){ |v| params[:olen] = v }
  opt.on('-1', '--read1 LEN', Integer, "Trim LEN bp from 5' end of read1 after adaptor trimming [#{params[:trim5]}]"){ |v| params[:trim5] = v }
  opt.on('-2', '--read2 LEN', Integer, "Trim LEN bp from 3' end of read2 after adaptor trimming [#{params[:trim3]}]"){ |v| params[:trim3] = v }
  opt.on('-A', '--trimA LEN', Integer, "Trim an additional LEN bp from the end of the adaptor sequence at the 1st stage (for RELACS library) [#{params[:trimA]}]"){ |v| params[:trimA] = v }
  opt.on('-f', '--first ADAPTOROPT', String, "Cutadapt adaptor sequence option at the 1st stage [#{params[:adapt_1st]}]"){ |v| params[:adapt_1st] = v }
  opt.on('-a', '--adaptor OPT', String, "flash2trimmer -a/--adaptors option at the 2nd stage [default]"){ |v| params[:adapt_opt] = v }
  opt.on('--times NUM', Integer, "Cutadapt --times option at the 1st stage"){ |v| params[:cutadapt_times] = "--times #{v}" }
  opt.on('-d', '--dir DIR', String, "Required if the source is '-'"){ |v| params[:dir] = v }
  opt.on('--chrMR DB', String, "Perform chrMRfilter2 filtering with DB after adaptor trimming"){ |v| params[:chrMR] = v }
  opt.on('-T', '--chrMR_minscore THRES', Integer, "BWA MEM minimum score option for chrMRfilter2"){ |v| params[:chrMR_minscore] = v }
  opt.on('-t', '--threads NUM', Integer, "Number of threads for Cutadapt and BWA MEM [#{params[:threads]}]"){ |v| params[:threads] = v }
  opt.on('-z', '--stdout', "Output to stdout"){ params[:stdout] = true }
  opt.on('--dry', "Dry run"){ params[:dry] = true }
  # opt.on('--check_file_order', "for detection of bug in version < 0.5"){ params[:check_file_order] = true }
  opt.parse!(ARGV)
}

case
when ARGV.size != 1
  $stderr.puts opts.to_s
  exit!
#when params[:check_file_order]
#  ARGV.collect{ |f| File.dirname(f) }.uniq.each do |dir|
#    f = [Dir.glob("#{dir}/#{File.basename(dir)}_?.fq.gz"),
#         Dir.glob("#{dir}/#{File.basename(dir)}_?.fastq.gz")].flatten
#    puts dir unless  f == f.sort
#  end
#  exit!
when ARGV.first == '-'
  unless params[:dir] 
    $stderr.puts "-d/--dir option is required" 
    exit!
  end
  oprefix = "#{params[:dir]}/#{File.basename(params[:dir])}"
  cmd = ''
else
  dir = ARGV.shift
  oprefix = "#{dir}/#{File.basename(dir)}"
  fq = [Dir.glob("#{oprefix}_?.fq.gz"), Dir.glob("#{oprefix}_?.fastq.gz") ].flatten.sort
  raise "Not paired-end fastq: #{fq.join(',')}" if fq.size != 2
  cmd = "seqtk mergepe #{fq.join(" ")} | "
end
  
tmpfile1 = Tmpfile.create
tmpfile2 = Tmpfile.create
system("touch #{tmpfile1} #{tmpfile2}")
tl = params[:trim5] + params[:trim3] + params[:trimA]
cmd += <<CMD.strip
cutadapt #{params[:adapt_1st]} -O #{params[:olen]} -j #{params[:threads]} --interleaved #{params[:cutadapt_times]} --action lowercase - 2> #{oprefix}.cutadapt.log.txt | 
mawk 'BEGIN{ OFS="\\n"} 
{ desc1 = $1;
  getline; seq1 = $0;
  getline; getline; qual1 = $0;
  getline; desc2 = $1;
  getline; seq2 = $0;
  getline; getline; qual2 = $0;
  trimmed = (sub("[atgcn]+$", "", seq1) + sub("[atgcn]+$", "", seq2));

  l = (length(seq1) < length(seq2) ? length(seq1) : length(seq2));
  s[length(seq1) "\\t" length(seq2) "\\t" trimmed]++;
  if (l < #{params[:minlen]+tl}) next;
  if (trimmed){
    print desc1, substr(seq1, 1+#{params[:trim5]}, l-#{tl}), "+", substr(qual1, 1+#{params[:trim5]}, l-#{tl}), 
      desc2, substr(seq2, 1+#{params[:trim3]}, l-#{tl}), "+", substr(qual2, 1+#{params[:trim3]}, l-#{tl}) > "#{tmpfile2}";
  } else if (length(seq1) == length(seq2)){
    print desc1, seq1, "+", qual1, desc2, seq2, "+", qual2;
  } else {
    print desc1, substr(seq1, 1, l), "+", substr(qual1, 1, l), desc2, substr(seq2, 1, l), "+", substr(qual2, 1, l);
  }
}
END {
  OFS="\\t";
  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
  print "read1len", "read2len", "trimmed", "n" > "/dev/stderr";
  for (k in s) print k, s[k] > "/dev/stderr";
}' 2> #{oprefix}.readlen_aligner.log.txt |
pigz > #{tmpfile1} &&
mawk -f #{File.dirname(__FILE__)}/flash2trimmer.awk -- -m #{params[:minlen]} #{params[:adapt_opt] ? "-a '#{params[:adapt_opt]}'" : ''} -1 #{params[:trim5]} -2 #{params[:trim3]} #{tmpfile1} 2> #{oprefix}.flash2trimmer.log.txt | 
cat - #{tmpfile2}
CMD
cmd += " | mawk -f #{File.dirname(__FILE__)}/chrMRfilter2.awk -- -t #{params[:threads]} -p -d #{params[:chrMR]} #{params[:chrMR_minscore] ? "-T #{params[:chrMR_minscore]}" : ''} 2> #{oprefix}.chrMRfilter2.log.txt" if params[:chrMR]
cmd += " | pigz -c > #{oprefix}.filtered.fq.gz" unless params[:stdout]

if params[:dry] 
  puts cmd
else
  system(cmd)
end




