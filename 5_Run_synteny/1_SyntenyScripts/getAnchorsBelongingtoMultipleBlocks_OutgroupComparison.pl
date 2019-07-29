# To filter the anchors belonging to multiple chromosomes, I get such 
# outgroup-polyploid anchors in which outgroup chr-gene pair belongs to 
# multiple human genes. Since these are anchors- they belong to a synteny blocks, and
# if one outgroup gene is shared in at least 2 anchors with different polyploid genes, 
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
	
	if ($line[2] eq '' || $line[3] eq ''){ # If there are no orthologs do this
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t-1\t0\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t\t\t\t\t\t\n";
	}
	else { # if there are orthologs for this ciona gene
		
		if (exists $Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}){ # And it exists in anchor with probability more than cutoff
			
			#my @line = split "\t", $_; # get individual componensts of the line
			
			my %check; # Declare a hash
			#print "$line[0]_$line[1]\n"; # outgroup gene and position
			
			foreach my $pair (@{$ci_genes_in_anchors{$line[0]."\t".$line[1]}}){ # For each polyploid gene-position pair
				
				$check{$pair} = \@line; # take them into this hash
				#print ">$pair\n";
			}
			#print "----------\n";
			
			#if ($Anchors{$_} < 2){print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t9\t9\t9\t";}
			#if ($Anchors{$_} == 2){print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t1\t0.5\t1\t";}
			#if ($Anchors{$_} == 3){print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t1\t0\t1\t";}
			#if ($Anchors{$_} > 3){print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t1\t0\t0\t";}
			
			print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t".${$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0]."\t"; # print genes and positions for Og and Pp; and synteny support in terms of color: ${$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0]
			
			#if (exists $Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}){
			#	print ${$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0],"\n";
			#}
			#else {print "$line[0]\t$line[1]\t$line[2]\t$line[3]\n";}
			
			# This is the important part to decide if this is an ohnolog candidate
			if (scalar keys %check > 1){ # If this outgroup gene_position pair is shared by multiple polyploid genes. It is a likely ohno candidate.
			
				my $answer = getAnchorsOnSameChr($_, $line[0], $line[1], $polyploidWindow, $shrink); # Check if it lies on the same chromosome, within the defined window.
				print OUT "$answer"; 
				#print "*** $answer\t$_\n";
				
				
			}
			else { # If it is not shared by multiple polyploid genes, it can't be an ohnolog.
				print OUT "0";
			}
			print OUT "\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t"; # Print rest of the things about this gene
			print OUT join ("\t" , @{$Anchors{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[1..6]),"\n"; # Print the probabilities from Anchor hash
			
			#if (scalar keys %check == scalar @{$ci_genes_in_anchors{$line[0].'_'.$line[1]}}){print "kk\n";}
			
			
		}
		else { # Else, check if it's there is Anchors with probability less than cutoff
			
			if (exists $AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}){
				print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t".${$AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[0]."\t0\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t";
				print OUT join ("\t", @{$AnchorsLessSignificant{"$line[0]\t$line[1]\t$line[2]\t$line[3]"}}[1..6]),"\n"; # Print from less significant anchors
				
			}
			else {
				print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t0\t0\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t\t\t\t\t\t\n";
			}
		}
	}
	
} 


# Subroutein to check if multiple polyploid gene_chr pairs lie on the same chromosome
sub getAnchorsOnSameChr {
	
	my $anchor = shift;
	my $cionaChr = shift;
	my $cionaGene = shift;
	my $wd = shift;
	my $shr = shift;
	#print "*$anchor\n";
	#print "$cionaChr\t$cionaGene\t\t";
	#print "@{$ci_genes_in_anchors{$cionaChr.'_'.$cionaGene}}\n";
	
	# I am here making a hash in which polyploid chromosomes for this outgroup pair are key so that i get unique chromosomes.
	my %check;
	foreach (@{$ci_genes_in_anchors{$cionaChr."\t".$cionaGene}}){	
		my  ($chr, $gn) = split "\t", $_, 2;
		$check{$chr} = '';
	}
	

	if (scalar keys %check == 1){ # If this outgroup gene_chr pair just maps to one chromosome
		
		#print ">$cionaChr\t$cionaGene\n";
		my @test;
		my $flag = 0;
		my $chromosome;
		
		# If this outgroup gene_chr pair just maps to one chromosome - get all the positions of genes in anchor in polyploid chromosome
		foreach (@{$ci_genes_in_anchors{$cionaChr."\t".$cionaGene}}){ 
					
			#print "$_\n";
			my ($c, $l) = split "\t", $_, 2;
			push @test, $l;
			$chromosome = $c;
		}
		
		# Sort the positions in ascending order
		@test = sort {$a <=> $b} @test;
		
		#print join '|', @test,"\n";;
		
		# Foreach of the position do the following
		for (my $i = 0; $i < scalar @test; $i++){
			
			if ($i+1 < (scalar @test)){
				
				#print "$test[$i]\t$test[$i+1]\n";
				
				# Calculate the difference
				#my $diff = $test[$i+1] - $test[$i];  # The distance for these windows not to overlap should be > wd size 
				#print "$diff\n";
				
				my $leftWd1; my $rightWd1; my $leftWd2; my $rightWd2;
				
				if ($shr eq 'y'){
					($leftWd1, $rightWd1) = getSymmatricWindowPosition($test[$i], $ppChrGenes{$chromosome}, ($wd/2)); # window size in this function should be half
					($leftWd2, $rightWd2) = getSymmatricWindowPosition($test[$i+1], $ppChrGenes{$chromosome}, ($wd/2));
				}
				if ($shr eq 'n'){
					($leftWd1, $rightWd1) = getAsymmatricWindowPosition($test[$i], $ppChrGenes{$chromosome}, ($wd/2)); # window size in this function should be half
					($leftWd2, $rightWd2) = getAsymmatricWindowPosition($test[$i+1], $ppChrGenes{$chromosome}, ($wd/2));
				}				
				
				#print "$leftWd1|$rightWd1\t\t$leftWd2|$rightWd2\n";
				#if ($diff > $wd){ # Because the anchors lie at the centre of the windows, difference between the anchors should be greater than twice the wd size , for the windows not to overlap
				if ($leftWd2 > $rightWd1){ # if the left window boundary of i+1 does not overlap with the right wd boundary of i - it's an ohno candidate 
					$flag = 1;
				}
			}
		}
		
		
		#print "$shr\t$flag\n*$anchor\n-----------------------------------\n";
		return($flag)
	}
	else { # If it maps to multiple polyploid chromosome, it must be an ohno candidate. Even if some of them map to a same chromosome.
		foreach (@{$ci_genes_in_anchors{$cionaChr."\t".$cionaGene}}){
			
			#print "$ciGenes\t$_\n";	
		}
		return(1);
	}
}


sub getAsymmatricWindowPosition {

	my $pos = shift;
	my $N = shift;
	my $wd = shift; # This is already a half window
	
	#print "$pos\t$N\t$wd\n";
	
	my $leftWd; my $rightWd;
	my $lOffset; my $rOffset;
	
	# generate regular symmatric windows
	$leftWd = $pos - $wd;
	$rightWd = $pos + $wd;
	
	# If left window has a problem reset it and the right wd
	if ($leftWd <= 0){
		
		$lOffset = abs($leftWd) + 1;
		$rightWd = $rightWd + $lOffset;
		$leftWd = 1;
	}
	
	# Fix the right window is its greater than chromosome
	if ($rightWd > $N){
			
		$rOffset = $rightWd - $N;
		$leftWd = $leftWd - $rOffset;
		$rightWd = $N;
	}
	
	# If still they exceed the chromosome limit - due to chromosome smaller than windows - reset them to boundaries
	if ($leftWd <= 0){$leftWd = 1};
	if ($rightWd > $N){$rightWd = $N};
	
	return ($leftWd, $rightWd);
}

sub getSymmatricWindowPosition {

	my $pos = shift;
	my $N = shift;
	my $wd = shift;
	
	my $leftWd; my $rightWd;
	
	# Get the window boundaries simply by adding the half window to both sides of the gene
	# So the total number of the genes will be complete window +1 
	# eg if position is at 755 window of size 200  will be 755+100 and 755-100 = 655 to 855 (total 201 genes including 655th and 855th gene)
	# For the windows to not overlap the anchors sitting in the middle should be more than wd size apart for symmatric windows and middle segments
	
	$leftWd = $pos - $wd;
	$rightWd = $pos + $wd;
	
	# If the window boundaries exceed the segment length reset them to segment boundaries
	#if ($leftWd <= 0){$leftWd = 1};
	#if ($rightWd > $N){$rightWd = $N};
	
	return ($leftWd, $rightWd);
	
}

print "[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}

