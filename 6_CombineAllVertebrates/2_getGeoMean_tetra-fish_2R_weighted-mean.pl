# This calculates the phylogenetically biased q-score and final ohnolog pairs for filtering.
# The output is all the candidate ohnologs combined from all outgroups, windos, and vertebrates.
# I calculate 2 q-scores for each probability -> a global with all vertebrates and one for each of the 3 major groups. We are using glocal q-score for now.
# ****************************************
use strict;
use warnings;
use Statistics::Descriptive;

local $| = 1;

# --------------------------------------------------------------------------------------------------------------------
# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 8: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'oanatinus', 'sscrofa', 'cfamiliaris', 'ggallus' # REMOE OANATINUS

my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

# Group of species for geometric mean
my @fishetc = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'loculatus'); # Fish etc
my @mammals = ('fcatus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'csabaeus'); # mammals
my @sauropsida = ('acarolinensis', 'ggallus', 'mgallopavo', 'tguttata');  # Lizard and birds
#my @primates = ('hsapiens', 'ptroglodytes', 'ggorilla', 'mmulatta', 'cjacchus', 'panubis', 'pabelii', 'csabaeus'); # primates. NOT USING IT ANY MORE
my $wgdtype = '_2R'; # _2R THIS IS ONLY FOR @R FOR 3R USE ANOTHER SCRIPT
# --------------------------------------------------------------------------------------------------------------------

foreach my $organism (@allorgs){
	
#if ($organism eq 'olatipes'){ # test for one organism

	print "> Combining $organism\n";
	my $infile = "$organism/$organism\_allVertebrates$wgdtype\.txt";
	my $outfile = "$organism/$organism\_allVertebrates$wgdtype\_withGeoMean.txt";
	print "  INPUT:  $infile\n";
	print "  OUTPUT: $outfile\n";

	open CONSENSUS, $infile or die $!; # Open Consensus file 
	open OUT, ">$outfile" or die $!; # Open outfile
	
	# ----------------------------------------------------------------------------------------
	# Arrays to hold the index of the groups --> This will be decided based on the header lines and then applied to rest of the lines
	my %fishindex; my %mamindex; my %birdindex; my %primateindex; 
	my %globalindex;

	# ----------------------------------------------------------------------------------------
	# Foreach gene
	my $linecount = 0; # This is just to print the stuff once inside the loop below
	while (<CONSENSUS>){
	
		my @line = split "\t", $_;
		map {$_=~s/\t|\n//g} @line;
		$linecount++;

		# ----------------------------------------------------------------------------------------
		# Print the header and get the index of the self and outgroup q-score for all organisms that needs to be averaged				
		if ($_ !~/^ENS.+/){
		
			# Check the order in the header file

			print OUT join ("\t", @line[0..6]),"\t";
			print OUT "$organism mult. P-og(>=k)\t$organism P-self(>=k)\t"; 
			print OUT "Og weighted global mean\tSelf weighted global mean\t"; # Global q-score
			print OUT join ("\t", @line[8..12]),"\t"; 
			print OUT join ("\t", @line[14..$#line]),"\t";
			print OUT "Og weighted group mean\tSelf weighted group mean\n"; # Group q-score
			
			# To get the index of the matches for different groups -- this index will be used to get the group specific geo-means
			my $refF = getIndex(\@line, \@fishetc); # FISH
			%fishindex = %{$refF};

			my $refM = getIndex(\@line, \@mammals); # MAMMALS
			%mamindex = %{$refM};

			my $refB = getIndex(\@line, \@sauropsida); # BIRTS LIZARDS ETC
			%birdindex = %{$refB};

			#my $refP = getIndex(\@line, \@primates); # PRIMATES, NOT USING IT ANY MORE
			#%primateindex = %{$refP};		
			
			# Global index with all the organisms
			my $refG = getIndex(\@line, \@allorgs); # ALL ORGANISMS
			%globalindex = %{$refG};
			
			next;
			
		};

		# Print some of the values, I'll print rest of the ones after global index below
		print OUT join ("\t", @line[0..7]),"\t";
		print OUT "$line[13]\t";

		
		# ----------------------------------------------------------------------------------------
		# Decide which group to use --> for fish group will be fish. For mammals -> mammals. For others -> birds
		
		my %groupindex; # this will have the index for the group and organism names as key
		my $weightfile = ''; # this will have the name of the weight file
		
		if ($organism ~~ @fishetc){
			%groupindex = %fishindex;
			$weightfile = 'all_fish_2R_filtered.txt';
		}
		#elsif ($organism ~~ @primates){ # NOT USING IT ANY MORE
		#	%groupindex = %primateindex;
		#	$weightfile = 'primates_2R_filtered.txt';
		#}
		elsif ($organism ~~ @mammals){
			%groupindex = %mamindex;
			$weightfile = 'all_mammals_2R_filtered.txt';
		}
		elsif ($organism ~~ @sauropsida){
			%groupindex = %birdindex;
			$weightfile = 'sauropsid_2R_filtered.txt';
		}
		else {die "Unidentified organism name!!\n";}
		
		# Print which group it is just once
		if ($linecount == 2){		
			if ($organism ~~ @fishetc){print "  -> Combining $organism with other FISH.\n";}
			#elsif ($organism ~~ @primates){print "  -> Combining $organism with other PRIMEATES\n";}
			elsif ($organism ~~ @mammals){print "  -> Combining $organism with other MAMMALS.\n";}
			elsif ($organism ~~ @sauropsida){print " -> Combining $organism with other SAUROPSIDAS.\n";}
		}			
		
		# ----------------------------------------------------------------------------------------
		# Open the weight file and read in the weights in GroupWeights hash
		my %GroupWeights;  # group weights

		open WF, "Analysis_for_q-score\/$weightfile" or die $!;
		foreach (<WF>){
			my @wghts = split "\t", $_;
			map {$_=~s/\t|\n//g} @wghts;
			#print "$wghts[1]\t$wghts[3]\n";
			$GroupWeights{$wghts[1]} = $wghts[3]

		}
		close (WF);
		
		# Open the weight file and read in the weights in GlobalWeights hash
		my %GlobalWeights; # global weights
		open GWF, "Analysis_for_q-score\/all_vertebrates_2R_filtered.txt" or die $!;
		foreach (<GWF>){
			my @wghts = split "\t", $_;
			map {$_=~s/\t|\n//g} @wghts;
			#print "$wghts[1]\t$wghts[3]\n";
			$GlobalWeights{$wghts[1]} = $wghts[3]

		}
		close (GWF);
				
		# ----------------------------------------------------------------------------------------		
		# Now add the values for the current organism that could not be added in the subroutein. See NOTE in the subroutein.
		$groupindex{$organism}{'og'} = 7; # Index for outgroup for self organism, check!
		$groupindex{$organism}{'self'} = 13; # Index for self for self organism, check!
		$globalindex{$organism}{'og'} = 7; # Index for outgroup for self organism, check!
		$globalindex{$organism}{'self'} = 13; # Index for self for self organism, check!


		# ----------------------------------------------------------------------------------------		
		# Calculate global q-score
		# Q-score formula is: qscore^k = Exp[ Sum_{species i in group k or global} p^k_i log(species i qscore)] (see email by HI)
		# For global q score the weights will be from all-vertebrates file
		my $sumogGlobal = 0; my $sumselfGlobal = 0;
		#print "$line[2] - $line[3]\tGLOBAL\n";
		# Foreach organism, weighted sum of global q-scores
		foreach my $organ (keys %globalindex){
						
			#print "$organ\t", $globalindex{$organ}{'og'},"\t", $globalindex{$organ}{'self'}, "\n";
			
			# Calculate for outgroups ---------------------			
			if ($line[$globalindex{$organ}{'og'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$globalindex{$organ}{'og'}], "\t", log($line[$globalindex{$organ}{'og'}]), "\t", $GlobalWeights{$organ}, "\t", (log($line[$globalindex{$organ}{'og'}]) * $GlobalWeights{$organ}), "\n";
				$sumogGlobal = $sumogGlobal + (log($line[$globalindex{$organ}{'og'}]) * $GlobalWeights{$organ})
			}

			# Calculate for self ---------------------		
			if ($line[$globalindex{$organ}{'self'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$globalindex{$organ}{'self'}], "\t", log($line[$globalindex{$organ}{'self'}]), "\t", $GlobalWeights{$organ}, "\t", (log($line[$globalindex{$organ}{'self'}]) * $GlobalWeights{$organ}), "\n";
				$sumselfGlobal = $sumselfGlobal + (log($line[$globalindex{$organ}{'self'}]) * $GlobalWeights{$organ})
			}			
		}
		
		# print global averages and rest of the lines
		print OUT exp($sumogGlobal), "\t";
		print OUT exp($sumselfGlobal), "\t";
		
		print OUT join ("\t", @line[8..12]),"\t";
		print OUT join ("\t", @line[14..$#line]),"\t";

		# test print global averages
		#print $sumogGlobal, "\t", exp($sumogGlobal), "\n";
		#print $sumselfGlobal, "\t", exp($sumselfGlobal), "\n";
		
		

		# ----------------------------------------------------------------------------------------		
		# Calculate group q-score
		# Q-score formula is: qscore^k = Exp[ Sum_{species i in group k or global} p^k_i log(species i qscore)] (see email by HI)
		my $sumog = 0; my $sumself = 0;
		#print "$line[2] - $line[3]\t GROUP\n";
		# Foreach organism, weighted sum group q-scores
		foreach my $organ (keys %groupindex){
						
			#print "$organ\t", $groupindex{$organ}{'og'},"\t", $groupindex{$organ}{'self'}, "\n";
			
			# Calculate for outgroups ---------------------			
			if ($line[$groupindex{$organ}{'og'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$groupindex{$organ}{'og'}], "\t", log($line[$groupindex{$organ}{'og'}]), "\t", $GroupWeights{$organ}, "\t", (log($line[$groupindex{$organ}{'og'}]) * $GroupWeights{$organ}), "\n";
				$sumog = $sumog + (log($line[$groupindex{$organ}{'og'}]) * $GroupWeights{$organ})
			}

			# Calculate for self ---------------------		
			if ($line[$groupindex{$organ}{'self'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$groupindex{$organ}{'self'}], "\t", log($line[$groupindex{$organ}{'self'}]), "\t", $GroupWeights{$organ}, "\t", (log($line[$groupindex{$organ}{'self'}]) * $GroupWeights{$organ}), "\n";
				$sumself = $sumself + (log($line[$groupindex{$organ}{'self'}]) * $GroupWeights{$organ})
			}
			
		}
		
		#print $sumog, "\t", exp($sumog), "\n";
		print OUT exp($sumog), "\t";
		#print $sumself, "\t", exp($sumself), "\n";
		print OUT exp($sumself), "\n";


	}
	close (CONSENSUS);

#} # end test for one organism
}



sub getIndex {
	
	my $lineref = shift; 
	my @line = @$lineref;
	my $orgref = shift;
	my @orgs = @$orgref;
	
	my @ind; my %index;
	foreach my $org (@orgs){
		
#		print "$org\t";
		@ind = grep { $line[$_] =~ /$org/ } 0..$#line; # The @ind here will have 3 columns that have the organims names. 0: gene symbols; 1: p_og; 2: p_self. I will extract the p_og and self later
#		print "*@ind*\n";
		$index{$org}{'og'} = $ind[1];
		$index{$org}{'self'} = $ind[2];
		#print "@index\n";
	}	
	return (\%index);
	# NOTE: Remember that the values for self organism here are blank. I will add these values later.
}


