# To assign new ids to peptides from version 1 of amphioxus assembly
# This new Id includes the location of the peptide on version 2 assembly 
# taken from gtf mapping file.

use strict;
use warnings;

open FH, "v2_Proteins_Ids_Scaffold.txt" or die $!;
my @scafIds = <FH>;
close(FH);
my %scaffoldHash;

foreach (@scafIds){
	
	my ($id, $scaf) = split " ", $_, 2;
	$scaf =~s/\n//g; 
	$scaffoldHash{$id} = $scaf;
}


local $/ = '>';
open FH2, '..\JGI_downloads\proteins.Brafl1.fasta' or die $!;
my @proteins = <FH2>;
close(FH2);
shift(@proteins);

open FH3, ">Aphioxus_Proteins_with_ScafIds.fasta" or die $!;
open FH4, ">Aphioxus_Proteins_Header_Mapping.txt" or die $!;

foreach (@proteins){
	
	my @lines = split "\n", $_;
	my $header = shift @lines;
	$header =~s/\n//g; 
	if ($lines[-1] =~/\>/g){pop @lines}
	
	$header =~ /jgi\|Brafl1\|(\d+)|.+/g;
	my $peptideId = $1;
	if (exists $scaffoldHash{$peptideId}){
		
		print FH3 ">$peptideId|$scaffoldHash{$peptideId}\n";
		print FH3 join "\n", @lines;
		print FH3 "\n";
	}
	else {print "$peptideId Not found\n";}
	print FH4 "$peptideId|$scaffoldHash{$peptideId}\t\t$header\n";
}

