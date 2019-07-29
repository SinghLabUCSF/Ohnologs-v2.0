# To filter the anchors belonging to multiple chromosomes, I get such 
# outgroup-polyploid anchors in which outgroup chr-gene pair belongs to 
# multiple human genes. Since these are anchors- they belong to a synteny blocks, and
# if one outgroup gene is shared by at least 2 anchors with different polyploid genes, 
# it means that it belongs to two regions in polyploid genome.
# I also make sure that if the polyploid chromosomes are same, the windows having anchors do not overlap.
#
# Input : 1. Ortholog relation file
#         2. Anchor file
#
# Output : Ohno candidate file having '1' in front of gene if it is an ohnolog candidate.
#

use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
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
my $shrink = $ARGV[8];
my $probability = $ARGV[9];
my $outFilesDir = $ARGV[10];

my $orthologRelationFile = $refGenomeName.'-'.$subGenomeName.'_Ortholog_Relations.txt';
#my $anchorFile = $refGenomeName.'-'.$subGenomeName.'_Anchors_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $anchorFileWithP = $refGenomeName.'-'.$subGenomeName.'_AnchorsWithProbability_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $ohnoCandidatesFile = $refGenomeName.'-'.$subGenomeName.'_OhnoCandidates_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $polyploidChrGenes = $subGenomeName.'_Chromosomes-Genes.txt';
#----------------------------------------------------------------------------------------------------------

#print "$refGenomeName\n$refGenomeFile\n$subGenomeName\n$subGenomeFile\n$orthologyFile\n$polyploidWindow\n$minOrthologs\n$orthologRelationFile\n$anchorFile\n$ohnoCandidatesFile\n";

print "Generating anchors belonging to multiple blocks...";

open OUT, ">$outFilesDir\/$ohnoCandidatesFile" or die $!;


# Read file with polyploid chromosomes and genes on them
open PP, "$outFilesDir\/$polyploidChrGenes" or die $!;
my @ppchrgn = <PP>;
close(PP);

my %ppChrGenes;
foreach (@ppchrgn){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$ppChrGenes{$line[0]} = $line[1];
}

open FH1, "$outFilesDir\/$orthologRelationFile" or die $!;
my @allHsGenes = <FH1>;
close (FH1);
shift @allHsGenes;

open FH2, "$outFilesDir\/$anchorFileWithP" or die $!;
my @anchors = <FH2>;
close (FH2);
my %Anchors; my %ci_genes_in_anchors; my %AnchorsLessSignificant;

foreach (@anchors){
	
	my ($color, $anc) = split '\|', $_, 2;
	$anc =~s/\s*Anchor\: //gs;
	$anc =~s/\n//gs;
	#print "*$color*\t*$anc*\n";
	my @line = split "\t", $anc;
	
	# Here i define two kinds of anchor hashes - one with anchors of probability more than cutoff, and one with less than cutoff
	
	if ($line[9] <= $probability){ # $line[9] = anchor probability, more than cutoff
		$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"} = [$color, $line[8], $line[9], $line[10], $line[11], $line[12], $line[13]]; # value contains color and the 6 probabilities
		#print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[9]\n";
		push @{$ci_genes_in_anchors{$line[0]."\t".$line[1]}}, $line[2]."\t".$line[3];
	}
	else { # less significant anchors
		$AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"} = [$color, $line[8], $line[9], $line[10], $line[11], $line[12], $line[13]];
		#print "* $line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[9]\n";
	}
}

#print scalar keys %Anchors;
#$print join "\n", keys %Anchors;
#print "@{$ci_genes_in_anchors{'7_579'}}";
#print join "\n", @{$ci_genes_in_anchors{'IV_184'}};

#print OUT "$refGenomeName Chr\t$refGenomeName Pos\t$subGenomeName Chr\t$subGenomeName Pos\tR\tG\tB\tOhno Candidate\t$refGenomeName Id\t$subGenomeName Id\t$refGenomeName Orientation\t$subGenomeName Orientation\n";

print OUT "$refGenomeName Chr\t$refGenomeName Pos\t$subGenomeName Chr\t$subGenomeName Pos\tColor\tOhno Candidate\t$refGenomeName Id\t$subGenomeName Id\t$refGenomeName Orientation\t$subGenomeName Orientation\t";
print OUT "P(chr)\tP(genome)\tP1(>=k)\tP2(>=k)\tP1\tP2\n";

foreach (@allHsGenes){ # Foreach line in the ortholog relation file
	
	my @line = split "\t", $_;
	$_=~s/\n//g;
	map {$_=~s/\n//g} @line;
	
	if ($line[2] eq '' || $line[3] eq ''){ # If there are no paralogs do this
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t-1\t0\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t\t\t\t\t\t\n";
	}
	else { # if there are paralogs for this ciona gene
		
		if (exists $Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}){ # And it exists in anchor with probability more than cutoff - it is an ohnolog. As I have already filtered genes in the same window
			
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t".${$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0]."\t"; # print genes and positions for Og and Pp; and synteny support in terms of color: ${$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0]
			
			print OUT "1"; # Print that it is ohno candidate
			print OUT "\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t"; # Print rest of the things about this gene
			print OUT join ("\t" , @{$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[1..6]),"\n"; # Print the probabilities from Anchor hash
			
			
		}
		else { # Else, check if it's there is Anchors with probability less than cutoff
			
			if (exists $AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}){
				print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t".${$AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0];
				print OUT "\t0\t"; # Print that it is not ohno candidate
				print OUT "$line[4]\t$line[5]\t$line[6]\t$line[7]\t"; # Rest of the things about this gene
				print OUT join ("\t", @{$AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[1..6]),"\n"; # Print probabilities from less significant anchors
				
			}
			else { # If not exist in anchors at all, it's not an ohno candidate
				print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t0\t0\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t\t\t\t\t\t\n";
			}
		}
	}
	
}

print "[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}

