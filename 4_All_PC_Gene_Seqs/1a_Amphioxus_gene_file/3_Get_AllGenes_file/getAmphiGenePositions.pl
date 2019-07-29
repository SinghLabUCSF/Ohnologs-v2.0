# To read the gtf file and get start end coordinates and gene ids for all the genes
# These I sort manually and use  for OHNOLOG run
use strict;
use warnings;

open GTF, '..\0_JGI_downloads\Bfloridae_v1.0_FilteredModelsMappedToAssemblyv2.0.gff' or die $!;
open OUT, '>AllGenes_Amphi_JGI_Unsorted.txt';

print OUT "Amphioxus Gene Id	Orientation	Scaffold	Start	End\n";

foreach (<GTF>){
	
	my @line = split "\t", $_;
	
	if ($line[2] eq 'mRNA'){
		$line[8] =~ /mRNA (\d+)\;.+\n/g;
		
		if ($line[6] eq '+'){$line[6] = '+1'}
		if ($line[6] eq '-'){$line[6] = '-1'}
		
		
		print OUT "$1\t$line[6]\t$line[0]\t$line[3]\t$line[4]\n";
	}
}



