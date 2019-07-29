# To get the positions for use in the plot programs
# The idea is to see if the ohnologs genes on chromosome lie within or outside the window
#
# Input : Ohno candidate file 
# Output: Synteny file

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
my $shrink = $ARGV[8];
my $outFilesDir = $ARGV[9];

my $ohnoCandidatesFile = $refGenomeName.'-'.$subGenomeName.'_OhnoCandidates_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $syntenyFileForPlotting = $refGenomeName.'-'.$subGenomeName.'_Synteny_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $polyploidChrGenes = $subGenomeName.'_Chromosomes-Genes.txt';
#----------------------------------------------------------------------------------------------------------

#print "$refGenomeName\n$refGenomeFile\n$subGenomeName\n$subGenomeFile\n$orthologyFile\n$polyploidWindow\n$minOrthologs\n$ohnoCandidatesFile\n$syntenyFileForPlotting\n";

print "Generating positions of ohnologs to plot...";

open OUTFILE, ">$outFilesDir\/$syntenyFileForPlotting" or die $!;



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



open FH, "$outFilesDir\/$ohnoCandidatesFile" or die $!;
my @file = <FH>;
close (FH);

my $header = shift @file;
chomp $header;

my %details; my %candidates; my @sort1;

foreach (@file){
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	push @sort1, \@line;
}

my $refS = sortKaro(\@sort1);
my @sortedOne = @$refS;

foreach (@sortedOne){

	my @line = @{$_};
	map {$_=~s/\n//g} @line;
	
	#print "$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\n";
	
	push @{$candidates{$line[0]."\t".$line[1]}}, $line[2]."\t".$line[3];
	$details{$line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]} = [$line[4], $line[5], $line[6], $line[7], $line[8], $line[9], $line[10], $line[11], $line[12], $line[13], $line[14], $line[15]];
}


my %toSort;

foreach my $ciona (keys %candidates){ # For each of ciona chromosome_gene pair
	
	#print "> Outgroup Chr_Gene Pair => Polyploid orthologs of this gene\n";
	#print "> $ciona\t=>\t";
	
	my @human = @{$candidates{$ciona}}; # Get all its human orthologs - on same or different chromosome
	
	#print join ("||", @human),"\n\n";
	
	# Take them in a hash where key is chromosome and value is all human genes on that chromosome
	# Since I am not taking care of the genes which do not have a human ortholog - for such gene it would be blank.
	my %chr;
	foreach (@human){
		
		my ($ch, $gene) = split "\t", $_, 2;
		#print "\t$ch -> $gene\n";
		push @{$chr{$ch}}, $gene;
	}
	
	#print "\tPolyChr -> Gene\n";
	foreach my $hchr (keys %chr){ # For each human chromosome. In case there is no ortholog the key would be blank and would be ignored later
		
		my @genes = @{$chr{$hchr}}; # Get all the gene on this chromosome
		my $count = 0;              # Initialize a position hash
		
		@genes = sort {$a <=> $b} @genes if ($genes[0] ne ''); # Sort the genes based on their position. Ignore the blanks i.e. genes without orthologs. If there is just one chromosome on this sort doesn't matter.
		#print "\t*$hchr*->(@genes)\n";
		
		for (my $i = 0; $i < scalar @genes; $i++){ # For each of the polyploid ortholog on this human chromosome
			
			# Only if there are multiple genes on this chromosome, the idea is to check if they lie
			# outside the window. For this i will take all sorted genes and check consequtive genes if the 
			# distance between them is greater than the overlap between windows, their position would be increased.
			# Or else it will remain the same
			if (scalar @genes > 1){
				
				#print "*$hchr*->(@genes)\n";
				my $diff = 0;          # Difference between 2 consequtive genes
				my $candidacy1 = 0;    # if 1st one is an ohnolog
				my $candidacy2 = 0;    # If 2nd one is an ohnolog
				my $leftWd1=0; my $rightWd1=0; my $leftWd2=0; my $rightWd2=0;
				
				if ($i+1 < scalar @genes){ # do this until the last gene is reached. this is to avoid the warning of array element overcounting
				
					if ($shrink eq 'y'){
						($leftWd1, $rightWd1) = getSymmatricWindowPosition($genes[$i], $ppChrGenes{$hchr}, ($polyploidWindow/2)); # Input: position, number of genes on this chromosome, half polyploid window size
						($leftWd2, $rightWd2) = getSymmatricWindowPosition($genes[$i+1], $ppChrGenes{$hchr}, ($polyploidWindow/2));
					}
					if ($shrink eq 'n'){
						($leftWd1, $rightWd1) = getAsymmatricWindowPosition($genes[$i], $ppChrGenes{$hchr}, ($polyploidWindow/2)); # window size in this function should be half
						($leftWd2, $rightWd2) = getAsymmatricWindowPosition($genes[$i+1], $ppChrGenes{$hchr}, ($polyploidWindow/2));
						#print "$genes[$i]\t$ppChrGenes{$hchr}\t", $polyploidWindow/2,"\n";
					}						
					
					$diff = $genes[$i+1] - $genes[$i];                                 # get the difference between consequtive genes. if there is just one gene the diff would be 0 as i+1 won't be less than array length
					$candidacy1 = ${$details{$ciona."\t".$hchr."\t".$genes[$i]}}[1];   # If the column definitions change the ohnocandidate column should change here
					$candidacy2 = ${$details{$ciona."\t".$hchr."\t".$genes[$i+1]}}[1]; # and here too	
					#print "$hchr\t| $genes[$i+1]| $genes[$i]| $diff|\n";
					
				}
	
				#print "\t| $genes[$i+1]| $genes[$i]| $diff|\n";
				#print "\t$candidacy1|$candidacy2\n";
				#print "$leftWd1|$rightWd1\t\t$leftWd2|$rightWd2\n";
				if ($leftWd2 > $rightWd1 && $candidacy1 == 1 && $candidacy2 == 1){ # if difference is alright and both are ohno candidate - increase the count
					#print "$ciona\t$hchr\t";
					#print "$genes[$i]\t";
					#print join "\t", (@{$details{$ciona."\t".$hchr."\t".$genes[$i]}});
					#print "\t$count\n";
	
					my ($cGene, $cId) = split "\t", $ciona, 2;
					my $str = join "\t", (@{$details{$ciona."\t".$hchr."\t".$genes[$i]}});
					push @{$toSort{$cGene}{$cId}}, "$cGene\t$cId\t$hchr\t$genes[$i]\t$str\t$count";
					$count++;
				
				}
				else {	#print "$ciona\t$hchr\t";
					#print "$genes[$i]\t";
					#print join "\t", (@{$details{$ciona."\t".$hchr."\t".$genes[$i]}});
					#print "\t$count\n";
					
					my ($cGene, $cId) = split "\t", $ciona, 2;
					my $str = join "\t", (@{$details{$ciona."\t".$hchr."\t".$genes[$i]}});
					push @{$toSort{$cGene}{$cId}}, "$cGene\t$cId\t$hchr\t$genes[$i]\t$str\t$count";
				}
			}
			else {	# If there is just a single gene push it into the array as is with position = 0
			
				my ($cGene, $cId) = split "\t", $ciona, 2;
				my $str = join "\t", (@{$details{$ciona."\t".$hchr."\t".$genes[$i]}});
				push @{$toSort{$cGene}{$cId}}, "$cGene\t$cId\t$hchr\t$genes[$i]\t$str\t$count";
			}
		}
	}
}


no warnings;
# Sort the array based on outgroup chromosome - gene and print in the file
print OUTFILE "$header\tPosition\n";
foreach my $cGn (sort {$a <=> $b} keys %toSort){
	
	my %geneHash = %{$toSort{$cGn}};
	
	foreach my $cCh (sort {$a <=> $b} keys %geneHash){
	
		print OUTFILE join ("\n", @{$geneHash{$cCh}}),"\n";
	}
}

# ---------------------
use warnings;

# Subroutein to perform a nested sort on multicolumn array based on different columns
sub sortKaro {
	
	my $ref = shift;
	my @a = @$ref;
	
	#print scalar @a;
	no warnings;	# I disable warnings in this block as I don't want this to shout at me just because I am using numeric sort on chromosome/scaffold, which can be alphanumeric also
	
	my @b = sort {   # if the column definitions change the order should change here

		$a->[0] <=> $b->[0] ||    # outgroup chromosome
		$a->[1] <=> $b->[1] ||    # outgroup position
		$a->[2] <=> $b->[2] ||    # polyploid chromosome
		$a->[5] <=> $b->[5] ||    # Ohnolog candidate 0/1
		$a->[3] <=> $b->[3]       # polyploid position

	} @a;
	
	#foreach (@b){	
	#	print join ("\t", @{$_}),"\n";
	#}
	
	return (\@b);
}

# THIS NEEDS HALF WINDOW
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

# THIS NEEDS HALF WINDOW TOO
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

print ".......[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}
