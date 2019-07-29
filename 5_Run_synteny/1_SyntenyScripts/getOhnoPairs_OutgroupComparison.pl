# Generate all possible ohnolog pairs and the ohnolog families. For human I am using gene ids as of now but
# It's better to use gene ids instead since not all the protein coding gene names are unique.
#
# Input : Synteny file with ohno candidacy and position 
# Output: 1. Ohno pairs file
#         2. Gene family file

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

my $syntenyFileForPlotting = $refGenomeName.'-'.$subGenomeName.'_Synteny_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $ohnoPairFile = $refGenomeName.'-'.$subGenomeName.'_OhnoPairs_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $FamilyFile = $refGenomeName.'-'.$subGenomeName.'_GeneFamilies_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';

print "Getting pairs of ohnologs...";

#print "$syntenyFileForPlotting";

open FH, "$outFilesDir\/$syntenyFileForPlotting" or die $!;
my @file = <FH>;
shift @file;

open OUT, ">$outFilesDir\/$ohnoPairFile" or die $!;
open OUTFAM, ">$outFilesDir\/$FamilyFile" or die $!;


# Get a hash of gene symbols -- I will use it for 
open IDS, "$subGenomeFile" or die $!;
my @ids = <IDS>; close (IDS);
shift (@ids);

# This is to print gene symbols for families. This is only for humans because for many invertebrates there is no gene symbol 
my %Ids_Symbol;
foreach (@ids){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	$Ids_Symbol{$line[0]} = $line[3];
}



my %pairs; 
foreach (@file){ # Foreach line of synteny file
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[5] == 1){ # If its ohnolog
		push @{$pairs{$line[0]}{$line[1]}}, \@line; # Push line to array with this outgroup chromosome-gene as key --- i.e. all PP genes sharing same outgroup gene
		#print "$line[0]\t$line[1]\n";
	}
}

foreach my $scaffold (keys %pairs){ # For each scaffold/chromosome
	
	#if ($scaffold eq '15'){
		
		foreach my $amphiPos (keys %{$pairs{$scaffold}}){ # For each position of outgroup 
		
			
			#if ($amphiPos == 1){
			
				#print "$scaffold---$amphiPos\t=>\t";
				#print "@{$pairs{$scaffold}{$amphiPos}}\n";
				getPairs(@{$pairs{$scaffold}{$amphiPos}});   # Get pairs using this subroutein
			#}
		}
	#}
}

# SUBROUTEIN to get the pairs given genes sharing the same outgroup gene -- Array of lines of these genes is starting point
# If they lie on different chromosome, they are ohnologs to each other. If on same chromosome, decide based on positions.
sub getPairs {
	
	my @lines = @_;
	my %decision; my %probabilities;
	
	# From the ohnologs to decide pairs and families - I make this hash to separate genes on same and different chromosome.
	# genes are sorted based on 'chromosome-positon' => @(gene list), so geens on same chromosome and within wd are together in the same array. 
	foreach(@lines){		
		#print "${$_}[11]\n";
		push @{$decision{${$_}[2]."\t".${$_}[-1]}}, ${$_}[7]; # I push the ensembl ids of the genes in this hash: desision.
		#print "${$_}[2]\t${$_}[-1]\t${$_}[7]\n";
		#  2 = Human chromosome
		# -1 = position - last element
		#  7 = Id
		$probabilities{${$_}[7]} = [${$_}[10], ${$_}[11], ${$_}[12], ${$_}[13], ${$_}[14], ${$_}[15]]; # Id and P(chr), P(genome) and P1
		
		#print "${$_}[7]\t${$_}[10]\n";
	}
	
	my @toPrint; my @check;
	
	print OUTFAM scalar keys %decision,"\t"; # print gene family size
	#print scalar keys %decision,"\t"; # print gene family size
	
	foreach (keys %decision){ # For each element of decision i.e. candidates which are ohnolog with one another
		
		# get gene symbols for families
		my @geneSymbols;
		foreach (@{$decision{$_}}){
			push @geneSymbols, $Ids_Symbol{$_};
			#print "$Ids_Symbol{$_}\t";
		}
		
		# Print the family where different elements are separated by \t, and SSD on same chromosome by pipe |
		#print OUTFAM join  ('|', @geneSymbols), "\t"; # Print gene symbols for families
		print OUTFAM join  ('|', @{$decision{$_}}), "\t"; # Print gene ids for family
		
		@toPrint = (@toPrint, join '|',@{$decision{$_}}); # To separate pairs of genes and infer their ohnolog relationship, I push them in this array.
	}
	print OUTFAM "\n";
	#print "\n";
	
	my $len = scalar @toPrint;
	
	# Get the pairs with SSDs
	for (my $i = 0; $i < $len; $i++){
		
		for (my $j = $i+1; $j < $len; $j++){
			
			push @check, "$toPrint[$i]\t$toPrint[$j]";
		}
	}
	
	# Separate SSDs and generated individual gene pairs
	foreach (@check){
		
		my @line = split "\t", $_;
			
		my @array1 = split '\|', $line[0];
		my @array2 = split '\|', $line[1];
		
		for (my $i = 0; $i < scalar(@array1); $i++){
				
			for (my $j = 0; $j < scalar(@array2); $j++){
					
				print OUT "$array1[$i]\t$array2[$j]\t";
				print OUT join ("\t", @{$probabilities{$array1[$i]}}),"\t";
				print OUT join ("\t", @{$probabilities{$array2[$j]}}),"\n";
				#print OUT eval ($probabilities{$array1[$i]} * $probabilities{$array2[$j]}),"\n";
			}
		}
	}
}


print "......................[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}

