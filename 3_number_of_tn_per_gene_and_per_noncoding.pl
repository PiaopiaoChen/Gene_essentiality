my $dir = "/gpfs/accounts/lsa_root/lsa/piaopiao/transposon/data/combined/genome";

my @list1 = ("MBE_YJM975","MBE_YJM975","MBE_RM11","MBE_DBVPG1107","MBE_WE372","MBE_DBVPG4651","MBE_322134S","MBE_322134S","MBE_322134S","MBE_S288C","MBE_YJM326","MBE_YJM454","MBE_CLIB219","MBE_CLIB219","MBE_CLIB219","MBE_YPS128","MBE_YPS128","MBE_UWOPS83-787.3","MBE_UWOPS83-787.3","NC_HN8","NC_HN8","NC_BJ14","NC_BJ20","NC_HN16");
my @list2 = ("s2_1","s2_2","s8","s10","s17","s19","s28_1","s28_2","s28_3","s31","s39","s45","s46_1","s46_2","s46_3","s50_1","s50_2","s58_1","s58_2","s67_1","s67_2","s71","s72","s75");

$i=0;
foreach $b (@list1){
  %hash1=();
  %hash2=();
open(FI,"$dir/5_${b}_with_gene");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
   if (($all[3] eq "yes")&&($all[4] eq "yes")){
       $loc = $all[0]."_".$all[1];
       $hash1{$loc}=$all[5];
	}
}
close(FI);

open(FI,"$dir/7_${b}_with_noncoding");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
    if (($all[3] eq "yes")&&($all[4] eq "yes")){
       $loc = $all[0]."_".$all[1];
       $hash2{$loc}=$all[5];
	}
}
close(FI);

open(IN,"$list2[$i].filter30_filterM_num.bed");
open(SE1,">1_$list2[$i]_with_gene");
open(SE2,">1_$list2[$i]_with_noncoding");
while($line=<IN>){
   chomp $line;
   @all=split/\t/,$line;
    $loc = $all[0]."_".$all[1];
	if (exists $hash1{$loc}){
	   print SE1 $line,"\t",$hash1{$loc},"\n";
	}
	else{
	   print SE1 $line,"\n";
	}
	
	if (exists $hash2{$loc}){
	   print SE2 $line,"\t",$hash2{$loc},"\n";
	}
	else{
	   print SE2 $line,"\n";
	}
}

print $b,"\t",$list2[$i],"\n";
$i++;
}

system("cat $dir/1_s2_1_with_gene $dir/1_s2_2_with_gene > $dir/1_s2_with_gene");
system("cat $dir/1_s2_1_with_noncoding $dir/1_s2_2_with_noncoding > $dir/1_s2_with_noncoding");

system("cat $dir/1_s28_1_with_gene $dir/1_s28_2_with_gene $dir/1_s28_3_with_gene > $dir/1_s28_with_gene");
system("cat $dir/1_s28_1_with_noncoding $dir/1_s28_2_with_noncoding $dir/1_s28_3_with_noncoding > $dir/1_s28_with_noncoding");

system("cat $dir/1_s46_1_with_gene $dir/1_s46_2_with_gene $dir/1_s46_3_with_gene > $dir/1_s46_with_gene");
system("cat $dir/1_s46_1_with_noncoding $dir/1_s46_2_with_noncoding $dir/1_s46_3_with_noncoding > $dir/1_s46_with_noncoding");

system("cat $dir/1_s50_1_with_gene $dir/1_s50_2_with_gene > $dir/1_s50_with_gene");
system("cat $dir/1_s50_1_with_noncoding $dir/1_s50_2_with_noncoding > $dir/1_s50_with_noncoding");

system("cat $dir/1_s58_1_with_gene $dir/1_s58_2_with_gene > $dir/1_s58_with_gene");
system("cat $dir/1_s58_1_with_noncoding $dir/1_s58_2_with_noncoding > $dir/1_s58_with_noncoding");

system("cat $dir/1_s67_1_with_gene $dir/1_s67_2_with_gene > $dir/1_s67_with_gene");
system("cat $dir/1_s67_1_with_noncoding $dir/1_s67_2_with_noncoding > $dir/1_s67_with_noncoding");



my @list = ("s2","s8","s10","s17","s19","s28","s31","s39","s45","s46","s50","s58","s67","s71","s72","s75");

foreach $b (@list){
 %hash1=();
 %hash1a=();
open(FI,"$dir/1_${b}_with_gene");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;	
   $len = @all;   
   if ($len > 5){
       @gene=split/_/,$all[5];
	   $num_gene = @gene;
	   $i=0;
	   while($i < $num_gene){
			   $hash1{$gene[$i]}++;
			   $hash1a{$gene[$i]} += $all[4];
			   $i++;
	   }
	}
}
close(FI);


 %hash2=();
 %hash2a=();
open(FI,"$dir/1_${b}_with_noncoding");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;	
   $len = @all;   
   if ($len > 5){
	      $hash2{$all[5]}++;
		  $hash2a{$all[5]} += $all[4];
	}
}
close(FI);


my $dir2 = "";

open(FI,"$dir2/all_gene_loc_length_exon");
open(SE,">$dir/2_${b}_number_of_tn_in_gene");
while($line=<FI>){
  chomp $line;
		   @all=split/\t/,$line;
		   if(exists $hash1{$all[4]}){
		      print SE $all[4],"\t",$hash1{$all[4]},"\t", $hash1a{$all[4]},"\n";
			}
			else{
			  print SE $all[4],"\t","0","\t","0","\n";
			}
}
close(FI);
close(SE);

open(FI,"$dir2/noncoding_region_1");
open(SE,">$dir/2_${b}_number_of_tn_in_noncoding");
while($line=<FI>){
           chomp $line;
           @all=split/\t/,$line;
		    if(exists $hash2{$all[0]}){
		      print SE $all[0],"\t", $hash2{$all[0]},"\t",$hash2a{$all[0]},"\n";
			}
			else{
			  print SE $all[0],"\t","0","\t","0","\n";
			}
}


}
















