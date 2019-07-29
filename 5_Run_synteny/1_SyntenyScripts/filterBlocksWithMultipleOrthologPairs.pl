# Filter the block file to remove the blocks which have less genes in the window.
# I filter blocks whcih have at at least defined number of genes within the window
# To count the synteny support I also use the approach to merge SSD's and define a 
# "color" parameter.
# 
# Input : Block file
# Output: Anchors belonging to multiple blocks
#

use strict;
use warnings;
use diagnostics; 

my $startTime = time;

=cut
my $refGenomeName = 'CionaInt';
my $refGenomeFile = 'AllPCGenes_CionaInt_Ens69.txt';
my $subGenomeName = 'Human';
my $subGenomeFile = 'AllPCGenes_Human_Ens69.txt';
my $orthologyFile = 'BestHits_Human-to-CionaInt.txt';
my $polyploidWindow = 250;
my $referenceWindow = 250;
my $minOrthologs  = 5;
=cut

my $refGenomeName = $ARGV[0];
my $refGenomeFile = $ARGV[1];
my $subGenomeName = $ARGV[2];
my $subGenomeFile = $ARGV[3];
my $orthologyFile = $ARGV[4];
my $referenceWindow = $ARGV[5];
my $polyploidWindow = $ARGV[6];
my $minOrthologs  = $ARGV[7];
my $outFilesDir = $ARGV[8];

my $outputBlocksFile = $refGenomeName.'-'.$subGenomeName.'_AllBlocks_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $anchorFile = $refGenomeName.'-'.$subGenomeName.'_Anchors_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $filteredBlocksFile = $refGenomeName.'-'.$subGenomeName.'_FilteredBlocks_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
#-------------------------------------------------------------------------------------------------------------------


#print "$refGenomeName\n$refGenomeFile\n$subGenomeName\n$subGenomeFile\n$orthologyFile\n$polyploidWindow\n$minOrthologs\n$outputBlocksFile\n$anchorFile\n";

print "Filtering blocks with multiple genes...";

local $/ = '>';
open FH, "$outFilesDir\/$outputBlocksFile" or die $!;
my @blocks = <FH>;
shift @blocks;
close (FH);

open OUTFILE, ">$outFilesDir\/$anchorFile" or die $!;
open OUT, ">$outFilesDir\/$filteredBlocksFile" or die $!;


foreach (@blocks){
	
	my @lines = split "\n", $_;
	pop @lines if (/\>/);
	my $anchor = shift @lines;
	#print $lines[-1],"\n"
	#print ">$anchor\n";
	my %humanGenes; my %cionaGenes; # hashes to identify unique human and ciona genes
	foreach (@lines){
		
		my @pair = split "\t", $_;
	 	map {$_=~s/\n|^\s+|\s+$//g} @pair;
		#print "$pair[1]\t$pair[3]\n";
		push @{$cionaGenes{$pair[1]}}, $pair[3]; # I am just pushing in an array. Check that there can't be a duplicate here 
		push @{$humanGenes{$pair[3]}}, $pair[1];
	}
	
	if (scalar keys %humanGenes >= $minOrthologs && scalar keys %cionaGenes >= $minOrthologs){ # Only blocks in which human/ciona genes equals or exceeds minimum ortholg required.		
				
		my $color = compressMultipleLinks(\%humanGenes, \%cionaGenes);
		
		print OUTFILE "$color|$anchor\n";
		print OUT ">$color|$anchor\n";
		print OUT join ("\n", @lines),"\n";
	}
}

close (OUTFILE);
close (OUT);

# This is to compress the multiple links
# The formula used is
# 1/Kj(Sum for all Ki(1/Ki))
#
# Where Kj = degree of polyploid gene (how many outgroup genes does it hit to)
# Ki = degree of outgroup gene (how many polyploid gene does it hit to)
#
# I have to calculate for each of the polyploid gene -> Degree of all outgroup orthologs 1/[Ki(i = 1..n)] and sum that
# Then divide the sum by Kj
#
# For best hit it is simple because there can't be any cross connections but for ensembl or other datasets where
# There can be 1-> many relationship from both sides, it will be useful.


sub compressMultipleLinks {
	
	my $hum = shift;
	my $ci = shift;
	
	my %human = %$hum;
	my %ciona = %$ci;
	
	#print "human\tciona\n";
	my $totalColor = 0;
	
	# Kj = Degree of the polyploid gene in the block. It will be always 1 for BestHits because polyploid gene can only point to 1 outgroup gene  
	# Ki = Degree of outgroup gene to which a polyploid gene points
	
	foreach (keys %human){ # Foreach polyploiud gene in the block
		
		# @{$human{$_}} has all the outgroup genes it points to
		my $Kj = scalar (@{$human{$_}}); # Get the degree of this gene i.e. how many gene it points to in outgroup (always 1 for besthit dataset)
	#	print "$_\t@{$human{$_}}\tKj = $Kj\t";
		
		# get sum of Ki
		my $sumOf1byKi = 0;
		foreach (@{$human{$_}}){ # Foreach outgroup gene it points to 
			
			my $Ki = scalar @{$ciona{$_}}; # get the degree i.e. polyploid gene it points too 
			$sumOf1byKi = $sumOf1byKi + (1/$Ki); # get the desired number
		}
		my $color = $sumOf1byKi/$Kj;
	#	print "\tSum of 1/ki =  $sumOf1byKi\tColor of this gene: $color\n";
		$totalColor = $totalColor + $color; # Sum the color
	}
	#print "Total Color = $totalColor\n";
	return($totalColor);
	
}

print "...........[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}
