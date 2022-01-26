my $dir = "";
my $dir2 = "";

my @list = ("MBE_YJM975","MBE_RM11","MBE_DBVPG1107","MBE_WE372","MBE_DBVPG4651","MBE_322134S","MBE_S288C","MBE_YJM326","MBE_YJM454","MBE_CLIB219","MBE_YPS128","MBE_UWOPS83-787.3","NC_HN8","NC_BJ14","NC_BJ20","NC_HN16");
my @list2 = ("s2","s8","s10","s17","s19","s28","s31","s39","s45","s46","s50","s58","s67","s71","s72","s75");


############length of the longest transposon-tolerant coding region within the gene
$z=0;

foreach $b (@list){
    %hash=();
    open(FI,"$dir/5_2_${b}_repeat_specific_gene");
    while($line=<FI>){
	   chomp $line;
	   @all=split/\t/,$line;
	   $hash{$all[0]}=$all[2];
    }
    close(FI);

  open(FI,"$dir/all_gene_loc_length_exon");
  open(SE,">$dir2/4_$list2[$z]_max_length_wo_tn");
  while($line=<FI>){
	   chomp $line;
	   @all=split/\t/,$line;
	   $i=1;
	   $max=0;
	   @value= split/_/,$hash{$all[4]};
	   $len =  @value;
       while($i <= $all[1]){
         @num=split/-/,$all[$i+7];
		 $start = $num[0]+1;
		 $end = $num[1]-1;
	     $j=$num[0]+1;
		 
	        open(IN,"$dir2/3_$list2[$z]_tn_hit_location_in_gene");
		    while($line2=<IN>){
	           chomp $line2;
	           @all2=split/\t/,$line2;
			    if(($all[6] eq $all2[0]) && ($all2[1] > $start) && ($all2[1] < $end)){ 
			      $diff1=$all2[1]-$j-1;
	           
				  $v=0;
				  $lost=0;
				  while($v <= $len){
   				    if(($value[$v] >= $j) && ($value[$v] <= $all2[1])){
					  $lost++;
					}
					$v++;
				  }
				  $diff=$diff1-$lost;
				 
				  if($diff > $max){
				    $max=$diff;
					$max_left=$j;
					$max_right=$all2[1];
				  } 
				  $j = $all2[1];
				}
			}
				    $last1=$end-$j-1;
				      $v=0;
				      $lost=0;
				      while($v <= $len){
   				       if(($value[$v] >= $j) && ($value[$v] <= $end)){
					    $lost++;
					   }
					   $v++;
				      }
				   $last=$last1-$lost;       
				   
				   if($last > $max){
				     $max=$last;
					 $max_left=$j;
					 $max_right=$end;
				   } 
				
            
			if ($max==0){
			    $max=$end-$start;
			}
			
            close(IN);
			$i++;
		}
		
		print SE $all[4],"\t",$max,"\t",$max_left,"\t",$max_right,"\n";
    }
	 $z++;
}


############transposon density in the gene divided by the transposon density in the surrounding 10,000 noncoding nucleotides (aka neighborhood index)


my @list1 = ("MBE_YJM975","MBE_YJM975","MBE_RM11","MBE_DBVPG1107","MBE_WE372","MBE_DBVPG4651","MBE_322134S","MBE_322134S","MBE_322134S","MBE_S288C","MBE_YJM326","MBE_YJM454","MBE_CLIB219","MBE_CLIB219","MBE_CLIB219","MBE_YPS128","MBE_YPS128","MBE_UWOPS83-787.3","MBE_UWOPS83-787.3","NC_HN8","NC_HN8","NC_BJ14","NC_BJ20","NC_HN16");
my @list2 = ("s2_1","s2_2","s8","s10","s17","s19","s28_1","s28_2","s28_3","s31","s39","s45","s46_1","s46_2","s46_3","s50_1","s50_2","s58_1","s58_2","s67_1","s67_2","s71","s72","s75");
$i=0;
foreach $b (@list1){
  %hash1=();

open(FI,"$dir/9_${b}_neighbour10000_noncoding");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
   if (($all[3] eq "yes")&&($all[4] eq "yes")){
       $loc = $all[0]."_".$all[1];
       $hash1{$loc}=$all[6];
	}
}
close(FI);


open(IN,"/$list2[$i].filter30_filterM_num.bed");
open(SE,">/5_$list2[$i]_noncoding_neighbour10000");

while($line=<IN>){
   chomp $line;
   @all=split/\t/,$line;
    $loc = $all[0]."_".$all[1];
	if (exists $hash1{$loc}){
	   print SE $line,"\t",$hash1{$loc},"\n";
	}
	else{
	   print SE $line,"\n";
	}
}

print $b,"\t",$list2[$i],"\n";
$i++;
}

system("cat $dir/5_s2_1_noncoding_neighbour10000 $dir/5_s2_2_noncoding_neighbour10000 > $dir/5_s2_noncoding_neighbour10000");

system("cat $dir/5_s28_1_noncoding_neighbour10000 $dir/5_s28_2_noncoding_neighbour10000 $dir/5_s28_3_noncoding_neighbour10000 > $dir/5_s28_noncoding_neighbour10000");

system("cat $dir/5_s46_1_noncoding_neighbour10000 $dir/5_s46_2_noncoding_neighbour10000 $dir/5_s46_3_noncoding_neighbour10000 > $dir/5_s46_noncoding_neighbour10000");

system("cat $dir/5_s50_1_noncoding_neighbour10000 $dir/5_s50_2_noncoding_neighbour10000 > $dir/5_s50_noncoding_neighbour10000");

system("cat $dir/5_s58_1_noncoding_neighbour10000 $dir/5_s58_2_noncoding_neighbour10000 > $dir/5_s58_noncoding_neighbour10000");

system("cat $dir/5_s67_1_noncoding_neighbour10000 $dir/5_s67_2_noncoding_neighbour10000 > $dir/5_s67_noncoding_neighbour10000");

my @list = ("s2","s8","s10","s17","s19","s28","s31","s39","s45","s46","s50","s58","s67","s71","s72","s75");

foreach $b (@list){
 %hash=();

open(FI,"$dir/5_${b}_noncoding_neighbour10000");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;	
   $len = @all;   
   if ($len > 5){
       @gene=split/_/,$all[5];
	   $num_gene = @gene;
	   $i=0;
	   while($i < $num_gene){
			   $hash{$gene[$i]}++;
			   $i++;
	   }
	}
}
close(FI);


my $dir2 = "";
open(FI,"$dir2/all_gene_loc_length_exon");
open(SE,">$dir/6_${b}_number_of_tn_in_noncoding_neighbour10000");
while($line=<FI>){
  chomp $line;
		   @all=split/\t/,$line;
		   if(exists $hash{$all[4]}){
		      print SE $all[4],"\t",$hash{$all[4]},"\n";
			}
			else{
			  print SE $all[4],"\t","0","\n";
			}
}
close(FI);
close(SE);

}
	
####number of transposons within the 100 nucleotides upstream of the gene


my @list1 = ("MBE_YJM975","MBE_YJM975","MBE_RM11","MBE_DBVPG1107","MBE_WE372","MBE_DBVPG4651","MBE_322134S","MBE_322134S","MBE_322134S","MBE_S288C","MBE_YJM326","MBE_YJM454","MBE_CLIB219","MBE_CLIB219","MBE_CLIB219","MBE_YPS128","MBE_YPS128","MBE_UWOPS83-787.3","MBE_UWOPS83-787.3","NC_HN8","NC_HN8","NC_BJ14","NC_BJ20","NC_HN16");
my @list2 = ("s2_1","s2_2","s8","s10","s17","s19","s28_1","s28_2","s28_3","s31","s39","s45","s46_1","s46_2","s46_3","s50_1","s50_2","s58_1","s58_2","s67_1","s67_2","s71","s72","s75");
$i=0;
foreach $b (@list1){
  %hash1=();

open(FI,"$dir/11_${b}_upstream100");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
   if (($all[3] eq "yes")&&($all[4] eq "yes")){
       $loc = $all[0]."_".$all[1];
       $hash1{$loc}=$all[5];
	}
}
close(FI);


open(IN,"/$list2[$i].filter30_filterM_num.bed");
open(SE,">/7_$list2[$i]_upstream100");

while($line=<IN>){
   chomp $line;
   @all=split/\t/,$line;
    $loc = $all[0]."_".$all[1];
	if (exists $hash1{$loc}){
	   print SE $line,"\t",$hash1{$loc},"\n";
	}
	else{
	   print SE $line,"\n";
	}
}

print $b,"\t",$list2[$i],"\n";
$i++;
}

system("cat $dir/7_s2_1_upstream100 $dir/7_s2_2_upstream100 > $dir/7_s2_upstream100");

system("cat $dir/7_s28_1_upstream100 $dir/7_s28_2_upstream100 $dir/7_s28_3_upstream100 > $dir/7_s28_upstream100");

system("cat $dir/7_s46_1_upstream100 $dir/7_s46_2_upstream100 $dir/7_s46_3_upstream100 > $dir/7_s46_upstream100");

system("cat $dir/7_s50_1_upstream100 $dir/7_s50_2_upstream100 > $dir/7_s50_upstream100");

system("cat $dir/7_s58_1_upstream100 $dir/7_s58_2_upstream100 > $dir/7_s58_upstream100");

system("cat $dir/7_s67_1_upstream100 $dir/7_s67_2_upstream100 > $dir/7_s67_upstream100");


my @list = ("s2","s8","s10","s17","s19","s28","s31","s39","s45","s46","s50","s58","s67","s71","s72","s75");

foreach $b (@list){
 %hash=();

open(FI,"$dir/7_${b}_upstream100");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;	
   $len = @all;   
   if ($len > 5){
       @gene=split/_/,$all[5];
	   $num_gene = @gene;
	   $i=0;
	   while($i < $num_gene){
			   $hash{$gene[$i]}++;
			   $i++;
	   }
	}
}
close(FI);


my $dir2 = "";
open(FI,"$dir2/all_gene_loc_length_exon");
open(SE,">$dir/8_${b}_number_of_tn_upstream100");
while($line=<FI>){
  chomp $line;
		   @all=split/\t/,$line;
		   if(exists $hash{$all[4]}){
		      print SE $all[4],"\t",$hash{$all[4]},"\n";
			}
			else{
			  print SE $all[4],"\t","0","\n";
			}
}
close(FI);
close(SE);

}
	

####number of gene segments without transposon insertion, each segment being one tenth of the coding sequence of the gene 
my @list1 = ("MBE_YJM975","MBE_YJM975","MBE_RM11","MBE_DBVPG1107","MBE_WE372","MBE_DBVPG4651","MBE_322134S","MBE_322134S","MBE_322134S","MBE_S288C","MBE_YJM326","MBE_YJM454","MBE_CLIB219","MBE_CLIB219","MBE_CLIB219","MBE_YPS128","MBE_YPS128","MBE_UWOPS83-787.3","MBE_UWOPS83-787.3","NC_HN8","NC_HN8","NC_BJ14","NC_BJ20","NC_HN16");
my @list2 = ("s2_1","s2_2","s8","s10","s17","s19","s28_1","s28_2","s28_3","s31","s39","s45","s46_1","s46_2","s46_3","s50_1","s50_2","s58_1","s58_2","s67_1","s67_2","s71","s72","s75");

$i=0;
foreach $b (@list1){
  %hash1=();

open(FI,"$dir/13_${b}_with_gene_bin10");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
   if (($all[3] eq "yes")&&($all[4] eq "yes")){
       $loc = $all[0]."_".$all[1];
       $hash1{$loc}=$all[5];
	}
}
close(FI);


open(IN,"/$list2[$i].filter30_filterM_num.bed");
open(SE,">/12_$list2[$i]_with_gene_bin10");

while($line=<IN>){
   chomp $line;
   @all=split/\t/,$line;
    $loc = $all[0]."_".$all[1];
	if (exists $hash1{$loc}){
	   print SE $line,"\t",$hash1{$loc},"\n";
	}
	else{
	   print SE $line,"\n";
	}
}

print $b,"\t",$list2[$i],"\n";
$i++;
}

system("cat $dir/12_s2_1_with_gene_bin10 $dir/12_s2_2_with_gene_bin10 > $dir/12_s2_with_gene_bin10");

system("cat $dir/12_s28_1_with_gene_bin10 $dir/12_s28_2_with_gene_bin10 $dir/12_s28_3_with_gene_bin10 > $dir/12_s28_with_gene_bin10");

system("cat $dir/12_s46_1_with_gene_bin10 $dir/12_s46_2_with_gene_bin10 $dir/12_s46_3_with_gene_bin10 > $dir/12_s46_with_gene_bin10");

system("cat $dir/12_s50_1_with_gene_bin10 $dir/12_s50_2_with_gene_bin10 > $dir/12_s50_with_gene_bin10");

system("cat $dir/12_s58_1_with_gene_bin10 $dir/12_s58_2_with_gene_bin10 > $dir/12_s58_with_gene_bin10");

system("cat $dir/12_s67_1_with_gene_bin10 $dir/12_s67_2_with_gene_bin10 > $dir/12_s67_with_gene_bin10");



my @list = ("s2","s8","s10","s17","s19","s28","s31","s39","s45","s46","s50","s58","s67","s71","s72","s75");


foreach $b (@list){
 %hash=();

open(FI,"$dir/12_${b}_with_gene_bin10");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;	
   $len = @all;   
   if ($len > 5){
       @gene=split/_/,$all[5];
	   $num_gene = @gene;
	   $i=0;
	   while($i < $num_gene){
			   $hash{$gene[$i]}++;
			   $i++;
	   }
	}
}
close(FI);


my $dir2 = "";
open(FI,"$dir2/all_gene_loc_length_exon");
open(SE,">$dir/13_${b}_number_of_tn_bin10");
while($line=<FI>){
  chomp $line;
		@all=split/\t/,$line;
		print SE $all[4];
		for($f = 1 ;$f < 11; $f = $f + 1 ){
		  $ke = $all[4].":".$f;
		    if(exists $hash{$ke}){
		      print SE "\t",$hash{$ke};
			}
			else{
			  print SE "\t","0";
			}
		}
		print SE "\n";
}
close(FI);
close(SE);

}
	















       