my $dir = "";
my %chr=
		(	chrI => "chr1",chrII => "chr2",chrIII => "chr3",chrIV => "chr4",chrV => "chr5",chrVI => "chr6",chrVII => "chr7",chrVIII => "chr8",chrIX => "chr9",chrX => "chr10",chrXI => "chr11",chrXII => "chr12",chrXIII => "chr13",chrXIV => "chr14",chrXV => "chr15",chrXVI => "chr16",chrM => "chr17"
		);
my %chr2=
		(	chr1 => "chrI",chr2 => "chrII",chr3 => "chrIII",chr4 => "chrIV",chr5 => "chrV",chr6 => "chrVI",chr7 => "chrVII",chr8 => "chrVIII",chr9 => "chrIX",chr10 => "chrX",chr11 => "chrXI",chr12 => "chrXII",chr13 => "chrXIII",chr14 => "chrXIV",chr15 => "chrXV",chr16 => "chrXVI",chr17 => "chrM"
		);
		
opendir (DIR, $dir) or die "can't open the directory!";
@dir = readdir DIR;

foreach $file (@dir) {
    if  ($file =~m /filterM\.bam$/){
	   print $file,"\n";
	   @name = split/_filterM/,$file;
	   @name2 = split/\./,$file;
	    open(FI,"$dir/$name[0]_filterM.bed");
	    open(SE,">$dir/$name[0]_filterM.bed2");
		print SE 'track name="',$name2[0],'"',' useScore=1',"\n";

		while($line=<FI>){
		   chomp $line;
		   @all=split/\t/,$line;
		   if ($all[5] eq "+"){
		       print SE $chr{$all[0]},"\t",$all[1]+1,"\t",$all[1]+2,"\t","1","\n";
		    }
		   elsif($all[5] eq "-"){
		       print SE $chr{$all[0]},"\t",$all[2]+1,"\t",$all[2]+2,"\t","0","\n";
			}
		}
	system("sort -n -k 4 -k 1.4 -k 2 $dir/$name[0]_filterM.bed2 > $dir/$name[0]_filterM.bed3");
		
		open(FI,"$dir/$name[0]_filterM.bed3");
	    open(SE,">$dir/$name[0]_filterM.bed4");
        $i=1;
		$j=0;
		$loc=0;
		
		while($line=<FI>){
		   chomp $line;
		   @all=split/\t/,$line;
		   if ($i == 1){
		       print SE $line,"\n";
		       
		    }
		   elsif ($i == 2){
		       $mark=$all[1];
			   $loc = $all[1]+$loc;
			   $chr=$all[0];
			   $j = $j + 1;
		    }
		   elsif(abs($all[1] - $mark) <= 2){
		       $j = $j + 1;
			   $loc = $all[1]+$loc;
			   $chr=$all[0];  
		    }
		   elsif(abs($all[1] - $mark) > 2){
		        $mean = $loc/$j;
		        print SE $chr,"\t",sprintf("%.0f", $mean),"\t",sprintf("%.0f", $mean)+1,"\t",$j,"\n";
				$mark = $all[1];
		        $j=0;
		        $loc=0;
				$loc = $all[1]+$loc;
				$chr=$all[0];
				$j = $j + 1;	
		    }
		   $i++;
        }
	system("sort -n -k 1.4 -k 2 -k 3 $dir/$name[0]_filterM.bed4 > $dir/$name[0]_filterM.bed5");		
	
    
	 open(FI,"$dir/$name[0]_filterM.bed5");
	 open(SE,">$dir/$name[0]_filterM.bed6");
		%hash={};
        $i=1;
		while($line=<FI>){
		   chomp $line;
		   @all=split/\t/,$line;
		   if ($i == 1){
		       print SE $line,"\n";
		    }
		   elsif($i > 1){
		       $loc = $all[0]."_".$all[1]."_".$all[2];
		       $hash{$loc} = $hash{$loc} + $all[3];
			}
			$i++;
		}
	
	for $a(sort keys %hash){
	   @all=split/_/,$a;
	   $num=$hash{$a} * 20+100;
	   if($num > 100){
	     print SE $all[0],"\t",$all[1],"\t",$all[2],"\t",".","\t",$num,"\n";
	   }
	}
	
	system("sort -n -k 1.4 -k 2 -k 3 $dir/$name[0]_filterM.bed6 > $dir/$name[0]_filterM.bed7");
	open(FI1,"$dir/$name[0]_filterM.bed7");
	open(SE1,">$dir/$name[0]_filterM_show.bed");
	open(SE2,">$dir/$name[0]_filterM_num.bed");
	$i=1;
		while($line=<FI1>){
		   chomp $line;
		   @all=split/\t/,$line;
		   if ($i == 1){
		       print SE1 $line,"\n";
			   print SE2 $line,"\n";
		    }
		   elsif($i > 1){
		       $num=($all[4]-100)/20;
		       print SE1 $chr2{$all[0]},"\t",$all[1],"\t",$all[2],"\t",$all[3],"\t",$all[4],"\n";
			   print SE2 $chr2{$all[0]},"\t",$all[1],"\t",$all[2],"\t",$all[3],"\t",$num,"\n";
			}
			$i++;
		}	
	    system("rm $dir/$name[0]_filterM.bed2");
	    system("rm $dir/$name[0]_filterM.bed4");	
        system("rm $dir/$name[0]_filterM.bed5");
        system("rm $dir/$name[0]_filterM.bed6");
        system("rm $dir/$name[0]_filterM.bed7");
				
		
	}
	
}
