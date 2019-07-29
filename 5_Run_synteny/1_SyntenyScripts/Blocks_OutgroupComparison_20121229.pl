# Get the blocks of conserved synteny having at least 2 different genes in the defined windows
# Input : Ortholog relation file
# Output: List of blocks of conserved synteny
# 
# Input : Ortholog relation file
# Output: Block file
#
# IMPORTANT : 1. The blocks are genes where there are at least 2 genes within the defined window

use strict;
use warnings;
use Data::Dumper;
use diagnostics; 

my $startTime = time;

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

my $outputBlocksFile = $refGenomeName.'-'.$subGenomeName.'_AllBlocks_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $orthologRelationFile = $refGenomeName.'-'.$subGenomeName.'_Ortholog_Relations.txt';
my $polyploidChrGenes = $subGenomeName.'_Chromosomes-Genes.txt';
my $outgroupChrGenes = $refGenomeName.'_Chromosomes-Genes.txt';
#----------------------------------------------------------------------------------------------------------

#print "$refGenomeName\n$refGenomeFile\n$subGenomeName\n$subGenomeFile\n$orthologyFile\n$orthologRelationFile\n$polyploidWindow\n$minOrthologs\n$outputBlocksFile\n$orthologRelationFile\n";

print "Getting syntenic blocks...";

open FH, "$outFilesDir\/$orthologRelationFile" or die $!;
my @asd = <FH>;
shift @asd;
close(FH);


open OUT, ">$outFilesDir\/$outputBlocksFile" or die $!;


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
#print scalar keys %ppChrGenes;


# Read file with outgroup chromosomes and genes on them
open OG, "$outFilesDir\/$outgroupChrGenes" or die $!;
my @ogchrgn = <OG>;
close(OG);


my %ogChrGenes;
foreach (@ogchrgn){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$ogChrGenes{$line[0]} = $line[1];
}
#print scalar keys %ogChrGenes;



# make the superhash
my %superHash;

foreach (@asd){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	
	
	if ($line[0] ne '' && $line[1] ne '' && $line[2] ne '' && $line[3] ne ''){
		$superHash{$line[0]}{$line[1]}{$line[2]}{$line[3]} = '';
		#print "$line[0]\t$line[1]\t$line[2]\t$line[3]\n";
	}
}
#print Dumper %superHash;

for (my $i = 0; $i < scalar (@asd); $i++){

	my @line = split "\t", $asd[$i];
	map {$_=~s/\n//g} @line;

	if ($line[0] ne '' && $line[1] ne '' && $line[2] ne '' && $line[3] ne ''){


		my @syntenyBlocks = getBlocks(\@line, $polyploidWindow, $referenceWindow, $shrink, \%superHash, \%ogChrGenes, \%ppChrGenes);

		if (scalar (@syntenyBlocks) >= 2){ # change number of genes in a block here

			print OUT ">Anchor: ";
			print OUT join ("\t", @line),"\n";

			foreach (@syntenyBlocks){	
				print OUT join ("\t", @{$_}),"\n";
			}
		}
	}
}


sub getBlocks{
	
	my $ref = shift;
	my $pwd = shift;
	$pwd = $pwd / 2;
	my $rwd = shift;
	$rwd = $rwd / 2;
	my $shrink = shift;
	my $supRef = shift;
	my %superHash = %$supRef;
	my $ogRef = shift;
	my %ogChrGenes = %$ogRef;
	my $ppRef = shift;
	my %ppChrGenes = %$ppRef;
	
	my @anchor = @{$ref};
	my @block = ();
	
	#if ($anchor[0] == 14 && $anchor[1]== 502  && $anchor[2] == 4 && $anchor[3] == 638){
	#print "$anchor[1]\t$anchor[3]\n";
	#print "Anchor: ";
	#print join ("\t", @anchor),"\n";
##	print "$pwd|$rwd\t\t";
	
	my $leftOutgroupWd; my $rightOutgroupWd;
	my $leftPolyploidWd; my $rightPolyploidWd;
	
	if ($shrink eq 'y'){  # Get Symmatric window. It will shrink window size at the corners.	
		
		($leftOutgroupWd, $rightOutgroupWd) = getSymmatricWindowPosition($anchor[1], $ogChrGenes{$anchor[0]}, $rwd); # For outgroup genome
		($leftPolyploidWd, $rightPolyploidWd) = getSymmatricWindowPosition($anchor[3], $ppChrGenes{$anchor[2]}, $pwd); # For polyploid genome
	
	}
	if ($shrink eq 'n'){ 	# Get Asymmatric window to not to shrink window size at the corners
	
		($leftOutgroupWd, $rightOutgroupWd) = getAsymmatricWindowPosition($anchor[1], $ogChrGenes{$anchor[0]}, $rwd); # For outgroup genome
		($leftPolyploidWd, $rightPolyploidWd) = getAsymmatricWindowPosition($anchor[3], $ppChrGenes{$anchor[2]}, $pwd); # For polyploid genome
	}
	
	my @OutgroupWindow = ($leftOutgroupWd..$rightOutgroupWd); # generate a window in the outgroup genome
	
#	print "Chr = $ogChrGenes{$anchor[0]}\tBoundary = $leftOutgroupWd|$rightOutgroupWd\n";
#	print "Chr = $ppChrGenes{$anchor[2]}\tBoundary = $leftPolyploidWd|$rightPolyploidWd\n";
#	print scalar @OutgroupWindow," ***\n";
	
	#my @OutgroupWindow = ($anchor[1]-$rwd..$anchor[1]+$rwd);
	
	# get the forward and a backward boundary limit in the polyploid genome
	#my $forwardLimit = $anchor[3] + $pwd;
	#my $backwardLimit = $anchor[3] - $pwd;
	
	#! THE BACKWARD LIMIT CAN BE IN NEGATIVE AND THE FW LIMIT CAN BE GREATER THAN THE SIZE OF CHROMOSOME BUT IT DOESN'T MATTER
	#! BECAUSE OF  THE CRITERIA THAT ORTHOLOG'S POSITION SHOULD BE GREATER THAN THE BW AND SMALLER THAN THE FW LIMIT.
	
##	print "$forwardLimit|$backwardLimit\n";
	
	#my @PolyploidWindow = ($anchor[3]-$pwd..$anchor[3]+$pwd); # and here
	
	map {if ($_ <= 0){$_ = ''}} @OutgroupWindow;
	#map {if ($_ <= 0){$_ = ''}} @PolyploidWindow;
	
	#print "@OutgroupWindow\n\n";
	#print "@PolyploidWindow\n";

	#for (my $i = 0; $i < scalar (@PolyploidWindow); $i++){
	#for (my $i = 0; $i < scalar (@PolyploidWindow); $i++){
		my $blockSize = 0;
		
		for (my $j = 0; $j < scalar (@OutgroupWindow); $j++){ # Foreach gene in the outgroup windows
			
			if (exists $superHash{$anchor[0]}{$OutgroupWindow[$j]}{$anchor[2]}){ # If there are human orthologs on current human chromosome
			
				my %humanPositions = %{$superHash{$anchor[0]}{$OutgroupWindow[$j]}{$anchor[2]}}; # get positions of those genes in this hash keys
				
#				print "* $anchor[0]\t$OutgroupWindow[$j]\t$anchor[2]\t"; #Ciona gene and position
				
				foreach (keys %humanPositions){ # And for each position
	
#						print " $_";
						
						if ($_ >= $leftPolyploidWd && $_ <= $rightPolyploidWd){ # Check if they lie within the polyploid window. Because the booundary is the part of this block i have the condition as >= and not >. Earlier in older version it was > which was wrong!
#							print "*";
							push @block, [$anchor[0], $OutgroupWindow[$j], $anchor[2], $_]; # Push them to block array
							$blockSize++;
						}
				}
#				print "\n";
			}
			
			#if (exists $superHash{$anchor[0]}{$OutgroupWindow[$j]}{$anchor[2]}{$PolyploidWindow[$i]}){
			#	push @block, [$anchor[0],$OutgroupWindow[$j],$anchor[2],$PolyploidWindow[$i]];
			#}
		}
##		print "Block Size \= $blockSize\n\n";
	#}
	#}
	return (@block);
}

sub getAsymmatricWindowPosition {

	my $pos = shift; # position
	my $N = shift;   # Total number of genes
	my $wd = shift;  # This is already a half window
	
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
	if ($leftWd <= 0){$leftWd = 1};
	if ($rightWd > $N){$rightWd = $N};
	
	return ($leftWd, $rightWd);
	
}


print "........................[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}
