# 09-June-2010
# This is to get all the Protein ids from filter model file for V02 amphioxus JGI
# Modified on 2016-09-29

use strict;
use warnings;

####### INPUT ##########

# Path and name of filter model file for V2
my $filterModelFile = '..\0_JGI_downloads\Bfloridae_v1.0_FilteredModelsMappedToAssemblyv2.0.gff';
# Path and Name of Output File
my $outfile = 'v2_Proteins_Ids_Scaffold.txt';

open FH1, $filterModelFile or die $!; 
open FH2, ">$outfile" or die $!; 

foreach (<FH1>){
	
	my @line = split "\t", $_;
	if ($line[2] eq 'mRNA'){
		$line[8] =~ /mRNA (\d+)\;.+\n/g;
		print FH2 "$1\t$line[0]\n";
	}
}
close (FH1);
close (FH2);