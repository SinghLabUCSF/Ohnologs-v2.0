# Calculate the probability that observed genes within the block occur by chance.
# I do it at 3 levels - 
# (1) Probability of finding all observed genes in given window in polyploid genome just by chance -- wrt the current chromosome
# (2) Probability of finding all observed genes in given window in polyploid genome just by chance -- wrt the entire genome
# (3) Probability of finding exactly k number of ANY observed genes within a given outgroup window
#
# To calculate first 2 - I will first calculate all windows on human chromosomes without a human ortholog of the currect outgroup
# gene, and all possible windows on human chromosome. The division is the P that no ortholog of current gene will be observed in any wd.
# 1-P then is the probability that you will observe at least 1 of the ortholog by chance. I calculate this probabilitty for each gene and multi-
# plication would be the probability that all genes in this block lie within at least one.


use strict;
use warnings;
use diagnostics;
use Math::Combinatorics;
use Math::BigFloat;

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

my $orthologRelationFile = $refGenomeName.'-'.$subGenomeName.'_Ortholog_Relations.txt';
my $filteredBlocksFile = $refGenomeName.'-'.$subGenomeName.'_FilteredBlocks_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $ppGeneCount = $subGenomeName.'_Chromosomes-Genes.txt';
my $ogGeneCount = $refGenomeName.'_Chromosomes-Genes.txt';
my $blocksWithProbability = $refGenomeName.'-'.$subGenomeName.'_FilteredBlocksWithProbability_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $anchorFileWithP = $refGenomeName.'-'.$subGenomeName.'_AnchorsWithProbability_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
#-------------------------------------------------------------------------------------------------------------------

print "Calculating probability of random occurrence...";



# open outfiles
open PROB, ">$outFilesDir\/$blocksWithProbability" or die $!;
open ANCPROB, ">$outFilesDir\/$anchorFileWithP" or die $!;

#----- get the sum for all Ks that herve wants --------   NOT USING RIGHT NOW
#open K1, '>Sum of P1 for all k.txt' or die $!;
#open K2, '>Sum of P2 for all k.txt' or die $!;


# read ortholog relation file
open OR, "$outFilesDir\/$orthologRelationFile" or die $!;
my @orf = <OR>;
shift @orf;
close(OR);

# and make the superhash. I need it to get all the orthologs not just ones in the blocks
my %superHash;

foreach (@orf){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[0] ne '' && $line[1] ne '' && $line[2] ne '' && $line[3] ne ''){
		$superHash{$line[0]}{$line[1]}{$line[2]}{$line[3]} = '';
	}
}

# Read WGD chromosome and gene count file
open PGC, "$outFilesDir\/$ppGeneCount" or die $!;
my @pgc = <PGC>;
close(PGC);

my %PolyploidGeneCount;
foreach (@pgc){
	#print "$_";
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$PolyploidGeneCount{$line[0]} = $line[1];
}



# Read outgroup chromosome and gene count file
open OGC, "$outFilesDir\/$ogGeneCount" or die $!;
my @ogc = <OGC>;
close(OGC);

my %OutgroupGeneCount;
foreach (@ogc){
	#print "$_";
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$OutgroupGeneCount{$line[0]} = $line[1];
}

#print join "\t", keys %OutgroupGeneCount;

# Read filtered blocks 
local $/ = '>';
open FH, "$outFilesDir\/$filteredBlocksFile" or die $!;
my @blocks = <FH>;
shift @blocks;
close (FH);


foreach (@blocks){
	
	my @genesInthisBlock = split "\n", $_;
	pop @genesInthisBlock if (/\>/);      # all genes in this block
	my $anchor = shift @genesInthisBlock; # anchor
	my @genesInBlockWithProbability;      # all genes in block with P appended at the end - used later
	
	#--- A few variables to store log p values for the block, and calculate final probability ----------------------#
	my $pAllGenesWithinWindow = 1;    # Multiplication of Ps for orthologs in Og-PP window for entire genome -- EXCLUDING ANCHOR
	my $pONEgenesWithinWindowChr; # Ps for orthologs in Og-PP window for only current chromosoms
	my $pAllGenesWithinWindowChr = 1; # Multiplication of Ps for orthologs in Og-PP window for only current chromosoms -- EXCLUDING ANCHOR
	my $pAnchorGenome;		  # P for all genes anchor FOR ENTIRE GENOME
	my $pAnchorChr;		  # P for all genes anchor FOR JUST CURRENT CHROMOSOME
	my $sumlogP = 0;              # Summation of log of probability values 
	my $sumlog1minusP = 0;        # Summation of log of 1 minus probability values
	my %ogGenesWithOrthInThisWd;  # Hash to identify unique Og genes having orthologs in this window ACROSS ENTIRE GENOME
	my %ogGenesWithOrthInThisWdChr;  # Hash to identify unique Og genes having orthologs in this window ONLY ON CURRENT CHROMOSOME
	my $k;                        # Total UNIQUE og genes having orthologs in current PP window -- EXCLUDING ANCHOR
	my $N;                        # Total UNIQUE og genes having orthologs anywhere in PP chromosomes -- EXCLUDING ANCHOR
	
#	print "$anchor\n";
	$anchor =~/Anchor\: (.+)\t(\d+)\t(.+)\t(\d+)\t.+\t.+\t.+\t.+/g; # get the chromosome and positions for the anchor
	my $ogAncChr = $1;
	my $ogAncPos = $2;
	my $ppAncChr = $3;
	my $ppAncPos = $4;
	
#	print "$ogAncChr $ogAncPos $ppAncChr $ppAncPos\n";
#	print "$ogAncPos, $OutgroupGeneCount{$ogAncChr}, ", $referenceWindow/2,"\n";
	
	#--- Get the window boundary for outgroup ---#
	my $leftWdPos; my $rightWdPos;
	
	if ($shrink eq 'n'){
		($leftWdPos, $rightWdPos) = getAsymmatricWindowPosition($ogAncPos, $OutgroupGeneCount{$ogAncChr}, $referenceWindow/2);
	}
	if ($shrink eq 'y'){
		($leftWdPos, $rightWdPos) = getSymmatricWindowPosition($ogAncPos, $OutgroupGeneCount{$ogAncChr}, $referenceWindow/2);
	}
	
	#--- For each outgroup gene position in window do the following ---#
	
	my $uniqueGenesWithOrth = 0; # This is total number of outgroup genes in current window having orthologs in WGD genome
	
	foreach my $position ($leftWdPos..$rightWdPos){ # Foreach position
		
		if (exists $superHash{$ogAncChr}{$position}){ # If it has orthologs
			
			$uniqueGenesWithOrth++;
			
			# Declare a hash having WGD chrs as key and empty array for positions			
			my %allPPChrPos;
			foreach (keys %PolyploidGeneCount){
				$allPPChrPos{$_} = [];
			}
#			print "$ogAncChr - $position => ";
			
			# Populate this hash with the ortholog polyploid chromosomes and the position array if any
			foreach my $ppChr (keys %{$superHash{$ogAncChr}{$position}}){
				
				$allPPChrPos{$ppChr} = [keys %{$superHash{$ogAncChr}{$position}{$ppChr}}];
			}
			
			# Now I have to calculate the probability of finding one ortholog of the given Og gene in at least one PP window
			# For that i have to get total possible windows on all PP chrs and total windows on PP chrs without any ortholog of this gene
			# And then 1 - the division of above 2 values is my Probability
			my $totalWdsWithoutAnyGenes = 0;
			my $totalWdsonALLChromosomes = 0;
			
			# to do that for each chromosome I get the possible windows and them sum them 
			foreach (keys %allPPChrPos){
				
				# pass the gene count on $_ chromosome and all genes on $_ chromosome
				# If there is no genes on any chromosome the empty array will be passed, as declared above 
				
				my $totWd; my $wdNoGene;
				
				if ($shrink eq 'y'){
					# Get the probability that this gene lies within any window by chance
					($wdNoGene, $totWd) = getSymmWindowsWithoutAnyOrtholog($PolyploidGeneCount{$_}, $polyploidWindow, \@{$allPPChrPos{$_}});  # Pass total genes on current chromosome, window size and array having all gene positions on this chromosome
				}
				if ($shrink eq 'n'){
					# Get the probability that this gene lies within any window by chance
					($wdNoGene, $totWd) = getAsymmWindowsWithoutAnyOrtholog($PolyploidGeneCount{$_}, $polyploidWindow, \@{$allPPChrPos{$_}}); # Pass total genes on current chromosome, window size and array having all gene positions on this chromosome
				}
				
				# printing outgroup chromosome|position and position of orthologs on all polyploid chromosomes if there is any ortholog, and possible windows
#				print "$ogAncChr|$position\t$_\t|@{$allPPChrPos{$_}}|\t$wdNoGene\t$totWd\n";
				
				
				# If the total possible windows for this chromosome is negative i.e. if the chr/scaff length is less than window, it will not be added to the calculations
				if ($totWd >= 0){
					$totalWdsWithoutAnyGenes = $totalWdsWithoutAnyGenes + $wdNoGene;
					$totalWdsonALLChromosomes = $totalWdsonALLChromosomes + $totWd;
				}
				
				#print "$_\t$totalWdsWithoutAnyGenes\t$totalWdsonALLChromosomes\n";
				
				# If the polyploid anchor is the current chromosome take it in a new variable 
				# to get the probability that we will find at least 1 orth ONLY ON CURRENT POLYPLOID CHROMOSOME randomly in the window
				if ($ppAncChr eq $_){
					
					#print "$_\t|@{$allPPChrPos{$_}}|\t$wdNoGene\t$totWd\n";
					
					# Total windows can be zero if the chromosome size is equal to window size. In that case there is division by 0 error
					# To avoid that if total possible windows is 0, just assign 1 to probability because all orthologs will be in that window anyway. 
					if ($totWd != 0){
						$pONEgenesWithinWindowChr = 1 - $wdNoGene/$totWd;
					}
					else {$pONEgenesWithinWindowChr = 1;}
					#print "$pONEgenesWithinWindowChr\n";
				}
				
			}
			
#			print "Total Wds without any gene = $totalWdsWithoutAnyGenes\nTotal possible windows = $totalWdsonALLChromosomes\n";
			# Get the probability of no gene within window
			my $pOFNoGeneInWindow = $totalWdsWithoutAnyGenes/$totalWdsonALLChromosomes;
#			print "p of no gene in window = $pOFNoGeneInWindow\n";
			# Get the p of at least this gene in at least one window (1 - p of no gene in window)
			my $pOfONEGeneInWindow = 1 - $pOFNoGeneInWindow;
#			print "p of at least ONE gene in window for = $pOfONEGeneInWindow\n\n";
			
			# Calculate the summation of all the probabilities and 1-p to calculate geometric mean later on --- TO EXCLUDE ANCHOR I SUBTRACT THE log(P) OF ANCHOR LATER ON FROM THESE SUM
			# Since log (0) gives error ignore zeros. Prob. can be zero if there is just 1 gene on a chromosome that is smaller than window. e.g. chromosome- Y
			# Also prob. can be zero if there are orthologs on all chromosomes and not even a single window can be placed without any gene. Ignore in that case too
			if ($pOfONEGeneInWindow > 0 && $pOfONEGeneInWindow < 1){
				$sumlogP = $sumlogP + log($pOfONEGeneInWindow);
				$sumlog1minusP = $sumlog1minusP + log(1 - $pOfONEGeneInWindow);
			}
			
			# Now I have to check which one of them exist in the block and push the probability at the end.
			# If there is no ortholog on any PP chr i.e. $_, i figured perl won't go in this loop. So only chromosomes having an ortholog are listed here.
			# If there are multiple polyploid orthologs on this gene, it doesn't matter as the probability will be same for them because they share same outgroup.
			foreach my $pChr (keys %allPPChrPos){
				
				foreach my $ppPos (@{$allPPChrPos{$pChr}}){
				
					#print "*$ogAncChr\t$position\t$pChr\t$ppPos\n";
					foreach my $blockLine (@genesInthisBlock){
				
						#print "$blockLine\n";
						if ($blockLine eq "$ogAncChr\t$position\t$pChr\t$ppPos"){
							#print "->$ogAncChr\t$position\t$pChr\t$ppPos\t$pOfONEGeneInWindow\t*$pONEgenesWithinWindowChr*\n";
							push @genesInBlockWithProbability, "$ogAncChr\t$position\t$pChr\t$ppPos\t$pONEgenesWithinWindowChr\t$pOfONEGeneInWindow";
							
							$ogGenesWithOrthInThisWd{$position} = $pOfONEGeneInWindow; # To get unique outgroup genes having orthologs ON ENTIRE PP GENOME in this window - push the og position in key.
							$ogGenesWithOrthInThisWdChr{$position} = $pONEgenesWithinWindowChr; # To get unique outgroup genes having orthologs ON THIS CHROMOSOME OF PP GENOME in this window - push the og position in key.
							
							if ($anchor =~ /$ogAncChr\t$position\t$pChr\t$ppPos/){
								$pAnchorGenome = $pOfONEGeneInWindow;
								$pAnchorChr = $pONEgenesWithinWindowChr;
						  	#	print "*$ogAncChr\t$position\t$pChr\t$ppPos\t$pOfONEGeneInWindow\n";
							}
							
						}
					}
				}
			}
#			print "\n----------------------------\n";
		}
		#print join "\n", @genesInthisBlock;
		#print "$uniqueGenesWithOrth\n";
	}
	
	#-- Calculate probability of finding all the unique Og genes in this window randomly ACROSS ENTIRE GENOME
	foreach (keys %ogGenesWithOrthInThisWd){
		# Here ignore the 0 values for individual genees. It can happen if there is just one ortholog on a chr/scaff which is smaller than the window size e.g. Y in human
		#print "$_ => $ogGenesWithOrthInThisWd{$_}\n";
		if ($ogGenesWithOrthInThisWd{$_} > 0){
			$pAllGenesWithinWindow = $pAllGenesWithinWindow * $ogGenesWithOrthInThisWd{$_};
		}
	}
	
	#-- Calculate probability of finding all the unique Og genes in this window randomly ON THIS CHROMOSOME ONLY
	foreach (keys %ogGenesWithOrthInThisWdChr){
		if ($ogGenesWithOrthInThisWdChr{$_} > 0){
		#print "===> $_ => $ogGenesWithOrthInThisWdChr{$_}\n";
			$pAllGenesWithinWindowChr = $pAllGenesWithinWindowChr * $ogGenesWithOrthInThisWdChr{$_};
		}
	}
	if ($pAnchorChr > 0){ # Remove the probability for anchor 
		$pAllGenesWithinWindowChr = $pAllGenesWithinWindowChr/$pAnchorChr;
	}
	
	
	# To exclude P of anchor divide by anchor probability for joint probability where i multiply
	# And subtract of log values where I add.
	# Now this could be problematic as p for anchor can be zero if anchor is on chr smaller than PP wd
	# Also if there are so many genes that its not possible to find a window without any gene ******************* THINK ABOUT IT ALSO IF I NEED TO IGNORE THAT ???
	# So i only remove ancchor if its P is greater than 0. Else i dont remove it.
	if ($pAnchorGenome > 0 && $pAnchorGenome < 1){
		$pAllGenesWithinWindow = $pAllGenesWithinWindow/$pAnchorGenome; # Since I have to EXCLUDE ANCHOR, I divide by anchor probability for join probability calculation
		$sumlogP = $sumlogP - log($pAnchorGenome);   # I deduct the log(P) for anchor from this sum 
		$sumlog1minusP = $sumlog1minusP - log(1 - $pAnchorGenome); # Here too I deduct the 1 - log(P) for anchor from this sum 
	}
	
#	print "* P of all genes within window: $pAllGenesWithinWindow\n";
	#my $pAllGenesWithinWindow = sprintf("%.3E", $pAllGenesWithinWindow);   # change to scientific notation

#	print "P of Anchor = $pAnchorGenome\n";
	$N = $uniqueGenesWithOrth - 1;         # Total outgroup genes having orthologs in PP genome-- EXCLUDING ANCHOR
	$k = scalar (keys %ogGenesWithOrthInThisWd) - 1;   # Total outgroup genes in this window having orthologs in PP window-- EXCLUDING ANCHOR
	my $LogpBar = $sumlogP/$N;             # Geomatric mean of log P for entire window
	my $P1 = exp ($LogpBar);               
	my $LogOneMinusPBar = $sumlog1minusP/$N; # Geomatric mean of log 1-P for entire window
	my $P2 = 1 - exp($LogOneMinusPBar);
		
	
		
	# Calculate Probability of observing any k number of genes within this window -- EXCLUDING ANCHOR	
	my $pAny_k_GenesWithinWindow_1; # From log P
	my $pAny_k_GenesWithinWindow_2; # From log 1-P
	
	# P from log approximation which Hervé told me to do earlier ---------- OLDER
	#$pAny_k_GenesWithinWindow_1 = (-1 * $k * (log($k/$N) - log($P1))) - (($N-$k) * (log(($N-$k)/$N) - log(1-$P1)));
	#$pAny_k_GenesWithinWindow_2 = (-1 * $k * (log($k/$N) - log($P2))) - (($N-$k) * (log(($N-$k)/$N) - log(1-$P2)));


	#--- Print probability in new anchor file, rest of the probabilities will be printed below when i evaluate the P for all possible K's
	print ANCPROB "$anchor\t$pAllGenesWithinWindowChr\t";
	print ANCPROB "$pAllGenesWithinWindow\t";
	
#	print ANCPROB exp($pAny_k_GenesWithinWindow_1),"\t";
#	print ANCPROB exp($pAny_k_GenesWithinWindow_2),"\n";
	
	#--- Print probability of block in new block file, sameway rest of the probabilities will be printed below when i evaluate the P for all possible K's
	print PROB ">$anchor\t$pAllGenesWithinWindowChr\t";
	print PROB "$pAllGenesWithinWindow\t";
#	print PROB exp($pAny_k_GenesWithinWindow_1),"\t";
#	print PROB exp($pAny_k_GenesWithinWindow_2),"\n";
	
#	foreach (@genesInBlockWithProbability){print PROB "$_\n";}

=cut
	print "N = $N\nk = $k\n";
	print "Sum log(p) = ", $sumlogP,"\n";
	print "Log P bar = Sum_log(P)/N = ",$LogpBar,"\n";
	print "P1 = exp (Log P bar) = ", $P1,"\n----------\n";
	print "Sum log (1-p) = ",$sumlog1minusP,"\n";
	print "Log 1 - Pbar = Sum log(1-p)/N = ",$LogOneMinusPBar,"\n";
	print "Exp (Log 1 - Pbar) = ",exp($LogOneMinusPBar),"\n";
	print "P2 = exp (1 - log(1 - Pbar)) = ", $P2, "\n";
=cut

	my $sumkP1 = 0; my $sumkP2 = 0; # Sum of P1 and P2 for all Ks
	my $P_MoreThanEqualK_1 = 0; my $P_MoreThanEqualK_2 = 0; # P of finding any K or more orthologs randomly
	
	foreach (0..$k){					
		
		my $n = Math::BigFloat->new($N); # Get the big float from N
		#my $n = $N;
		my $kLocal = $_;
		
		my $N_C_k = $n->bnok($kLocal); # Calculate "N choose k" from this function from library of Math::bigFloat
#		print "$kLocal\t$N\t$N_C_k\n";
	
		my $P1k = $P1 ** $kLocal;
		my $P2k = $P2 ** $kLocal;
#		print "\t$P1k\t$P2k\n";
		
		my $OneMinusP1k = (1-$P1) ** ($N-$kLocal);
		my $OneMinusP2k = (1-$P2) ** ($N-$kLocal);	
#		print "\t$OneMinusP1k\t$OneMinusP2k\n";
		
		# Convert it to float
		my $x1 = Math::BigFloat->new($N_C_k);
		my $x2 = Math::BigFloat->new($P1k);
		my $x3 = Math::BigFloat->new($OneMinusP1k);
		my $x4 = Math::BigFloat->new($P2k);
		my $x5 = Math::BigFloat->new($OneMinusP2k);
=cut		
		my $x1 = $N_C_k;
		my $x2 = $P1k;
		my $x3 = $OneMinusP1k;
		my $x4 = $P2k;
		my $x5 = $OneMinusP2k;
=cut

		#Math::BigFloat->accuracy(5);
		#print "\n->$x1\n$x2\n$x3\n";
		my $testP1 = $x1 * $x2 * $x3; # That's the P for current k
		my $testP2 = $x1 * $x4 * $x5;
		
#		print "* $testP1\n* $testP2\n";
		#$testP1 = sprintf("%.3E", $testP1);
		#$testP2 = sprintf("%.3E", $testP2);
		#print "-> $testP1\t$testP2\n";
		
		# Sum P1 and P2's
		$sumkP1 = $sumkP1 + $testP1;
		$sumkP2 = $sumkP2 + $testP2;
		
		if ($kLocal == ($k-1)){ # For all P's and Sum of P below actual K
			
			$P_MoreThanEqualK_1 = 1 - $sumkP1;
			$P_MoreThanEqualK_2 = 1 - $sumkP2;
			#print "$kLocal\n$sumkP1\n$sumkP2\n\n";
			print ANCPROB "$P_MoreThanEqualK_1\t$P_MoreThanEqualK_2\t"; # Prinf P of finding any K or more genes in the block
			print PROB "$P_MoreThanEqualK_1\t$P_MoreThanEqualK_2\t"; # Prinf P of finding any K or more genes in the block
		}
		
		if ($k == $kLocal){
		
			print ANCPROB "$testP1\t$testP2\n"; # P of finding exactly any K gene randomly in the block
			print PROB "$testP1\t$testP2\n";    # P of finding exactly any K gene randomly in the block
			foreach (@genesInBlockWithProbability){print PROB "$_\n";}
			#print "-> $testP1\t$testP2\n";
		}
		#print "$kLocal\t$testP1\t$testP2\t$sumkP1\t$sumkP2\n";
		

	}
	
	#my $sumkP1_Formatted = sprintf("%.3E", $sumkP1);
	#my $sumkP2_Formatted = sprintf("%.3E", $sumkP2);
	
	#print K1 "$anchor\t$sumkP1_Formatted\t$sumkP1\n";
	#print K2 "$anchor\t$sumkP2_Formatted\t$sumkP2\n";
		
}


# In case if the window is symmatric do the following
sub getSymmWindowsWithoutAnyOrtholog {

	my $totalGenes = shift;
	my $wdSize = shift;
	my $ref = shift;
	my @genes = @$ref;

	@genes = sort {$a <=> $b} @genes;
	
	
#	print "total genes = $totalGenes\twindow = $wdSize\tTotal orth on this chr = ";
#	print scalar keys @genes, "\tOrthologs = ";
	
	# Since I have to calculate the segments in between these genes from start of the 
	# chromosome till end - i am appending start position = 0 and end position = $totalGenes+1
	# at the beginning and end of this array genes, only if 1st and last element is not already 1 or length of gene.
	push @genes, ($totalGenes + 1);
	

#	print join ('|', @genes),"\n";
	
	my $start = 0; # first position of gene
	my $totalPossibleWindows = 0; # total possible windows without any ortholog of this gene
	for (my $i = 0; $i < scalar @genes; $i++){
		
		my $segment = $genes[$i] - $start - 1;
#		print "Start = $start\tGenePos = $genes[$i]\tSegment = $segment\n";
		if ($start == 0 || $genes[$i] == ($totalGenes + 1)){  # If it is the first or last end segment. Window will shrink and the would be half for symmatric window
			
			if ($segment > ($wdSize/2)){
				$totalPossibleWindows = $totalPossibleWindows + ($segment - ($wdSize/2));
#				print "* Windows = $totalPossibleWindows\n";
			}
			
		}
		if ($start != 0 && $genes[$i] != ($totalGenes + 1)) {  # If it is middle segment the window will not shrink
			
			if ($segment > $wdSize){
				$totalPossibleWindows = $totalPossibleWindows + ($segment - $wdSize);
#				print " -> Windows = $totalPossibleWindows\n";
			}
		}
		
		
		$start = $genes[$i];
		
	}
#	print "Total Windows for this gene = $totalPossibleWindows\n";
	my $probabilityOf_NONE_OfGenesWithinWindow = $totalPossibleWindows/$totalGenes;
	my $probabilityOf_ONE_OfGeneWithinWindow = 1 - $probabilityOf_NONE_OfGenesWithinWindow;
	
#	print "$probabilityOf_NONE_OfGenesWithinWindow\t$probabilityOf_ONE_OfGeneWithinWindow\n";
	
	return($probabilityOf_ONE_OfGeneWithinWindow);
	
}

# In case of Asymmatric window, do the following
sub getAsymmWindowsWithoutAnyOrtholog {

	my $totalGenes = shift;
	my $wdSize = shift;
	my $ref = shift;
	my @genes = @$ref;

	@genes = sort {$a <=> $b} @genes;
	
	
#	print "total genes = $totalGenes\twindow = $wdSize\tTotal orth on this chr = ";
#	print scalar keys @genes, "\tOrthologs = ";
	
	# Since I have to calculate the segments in between these genes from start of the 
	# chromosome till end - i am appending start position = 0 and end position = $totalGenes + 1
	# at the beginning and end of this array genes
	push @genes, ($totalGenes + 1);
	
	#print ">";
	#print join ('|', @genes),"\n";
	
	my $start = 0; # first position of gene
	my $totalPossibleWindows = 0; # total possible windows WITHOUT ANY ORTHOLOG of this gene
	for (my $i = 0; $i < scalar @genes; $i++){
		
		my $segment = $genes[$i] - $start - 1;
#		print "Start = $start\tGenePos = $genes[$i]\tSegment = $segment\n";

		if ($segment > $wdSize){ # For asymmatric window the window size will always be same so no shrinking
			$totalPossibleWindows = $totalPossibleWindows + ($segment - $wdSize);
#			print "Windows = $totalPossibleWindows\n";
		}
		
		$start = $genes[$i];
	}
#	print "Total Windows for this gene = $totalPossibleWindows\n";
	
#	print "$probabilityOf_NONE_OfGenesWithinWindow\t$probabilityOf_ONE_OfGeneWithinWindow\n";
	
	return($totalPossibleWindows, ($totalGenes-$wdSize));
	
}


# The same subrouteins to get the position of the asymmatric and symmatric window to get total og genes in the window having orthologs anywhere on the genome
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

print "...[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}
