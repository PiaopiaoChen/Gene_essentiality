#perl

$bwa="bwa-0.7.17";
$bwa_index = "";
$picard = "picard.jar";
$GATK = "GenomeAnalysisTK.jar";
$GATK_index = "";
$samtools = "samtools-1.8";
$bedtools = "";

my $dir = "";

opendir (DIR, $dir) or die "can't open the directory!";
@dir = readdir DIR;
   foreach $file (@dir) {
	    if  ($file =~m /fastq$/){
        @name = split/\.fastq/,$file;
  #trim
system("cutadapt -a GATC $dir/$name[0]_dpn$name[1].fastq -m 20 --overlap 4 -o $dir/$name[0]_dpn$name[1].cut.fastq");
system("cutadapt -a CATG $dir/$name[0]_nla$name[1].fastq -m 20 --overlap 4 -o $dir/$name[0]_nla$name[1].cut.fastq");
system("cat $dir/$name[0]_dpn$name[1].cut.fastq $dir/$name[0]_nla$name[1].cut.fastq > $dir/$name[0]$name[1].fastq");

  #bwa2-align;
system("$bwa/bwa mem -t 2 $bwa_index $dir/$name[0]$name[1].fastq > $dir/$name[0]$name[1].sam");
system("$samtools/samtools view -bS -h $dir/$name[0]$name[1].sam > $dir/$name[0]$name[1].bam");
  
  #picard sort bam
system("java -Xmx3g -jar $picard SortSam I=$dir/$name[0]$name[1].bam O=$dir/$name[0]$name[1].sort.bam SORT_ORDER=coordinate");
system("$samtools/samtools index $dir/$name[0]$name[1].sort.bam");


  #add filter
system("$samtools/samtools view -b -h -F 4 -F 256 -F 2048 -q 30 $dir/$name[0]$name[1].sort.bam > $dir/$name[0]$name[1].filter30.bam");
system("$samtools/samtools index $dir/$name[0]$name[1].filter30.bam");
system("$samtools/samtools view -h $dir/$name[0]$name[1].filter30.bam -o $dir/$name[0]$name[1].filter30.sam");
}
}

##reads filter;

opendir (DIR, $dir) or die "can't open the directory!";
@dir = readdir DIR;
foreach $file (@dir) {
 if  ($file =~m /filter30\.sam$/){
       @name = split/\.sam/,$file;

open(FI,"$dir/$name[0].sam");
open(SE1,">$dir/$name[0]_filterM.sam");
open(SE2,">$dir/$name[0]_filterM_removed.sam");
$i = 1; 
while($line=<FI>){	 
    chomp $line;
    if ($i < 20){
	   print SE1 $line,"\n";
	 }
	else{
      @all=split/\t/,$line;
	   if (($all[1]==0) && ($all[5] =~ /^\d+M/)){
	     print SE1 $line,"\n";
	   }
	   elsif (($all[1]==16) && ($all[5] =~ /\d+M$/)){
	     print SE1 $line,"\n";
	   }
	   else{
	     print SE2 $line,"\n";
		}

	}
	$i++;
}
	   
system("samtools view -bS -h $dir/$name[0]_filterM.sam > $dir/$name[0]_filterM.bam");
system("samtools index $dir/$name[0]_filterM.bam");
system("$bedtools bamtobed -i $dir/$name[0]_filterM.bam > $dir/$name[0]_filterM.bed");
}
}

