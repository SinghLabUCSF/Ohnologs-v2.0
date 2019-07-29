# This calculates the phylogenetically biased q-score and final ohnolog pairs for filtering.
# The output is all the candidate ohnologs combined from all outgroups, windos, and vertebrates.
# Remember that duplication time in some of the ohnologs that are frm other vertebrates may be different than WGD nodes, I am keeping them for now.
# There are some warnings that correspond to weird lines in the input file, maybe a filesystem issue. Doesn't affect the stats, so I'll look in detail later.
# ****************************************
use strict;
use warnings;
use Statistics::Descriptive;

# --------------------------------------------------------------------------------------------------------------------
# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 8: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'oanatinus', 'sscrofa', 'cfamiliaris', 'ggallus'

my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

my $wgdtype = '_3R'; # _2R or _3R
if ($wgdtype ne '_3R'){die "Incorrect WGD type!\n";}

# --------------------------------------------------------------------------------------------------------------------
# Calculate weighted mean for each organism
# Remember that group and global mean is the same for this. Or no group mean as last 2 columns.
foreach my $organism (@fish){

#if ($organism eq 'olatipes'){ # test for one organism
	
	print "> Combining $organism\n";
	my $infile = "$organism/$organism\_allVertebrates$wgdtype\.txt";
	my $outfile = "$organism/$organism\_allVertebrates$wgdtype\_withGeoMean.txt";
	print "INPUT:  $infile\n";
	print "OUTPUT: $outfile\n\n";

	open CONSENSUS, $infile or die $!; # Open Consensus file 
	open OUT, ">$outfile" or die $!; # Open outfile
	
	# Arrays to hold the index of all fish --> This will be decided based on the header lines and then applied to rest of the lines
	my %fishindex; 
	
	# ----------------------------------------------------------------------------------------
	# Foreach gene	
	while (<CONSENSUS>){
	
	
		my @line = split "\t", $_;
		map {$_=~s/\t|\n//g} @line;

		if ($_ !~/^ENS.+/){
		
			# Check the order in the header file

			print OUT join ("\t", @line[0..6]),"\t";
			print OUT "$organism mult. P-og(>=k)\t$organism P-self(>=k)\t";
			print OUT "Og weighted global mean\tSelf weighted global mean\t";
			print OUT join ("\t", @line[8..14]),"\t";
			print OUT join ("\t", @line[16..$#line]),"\n";
			
			# Fish index with all the organisms
			my $refF = getIndex(\@line, \@fish); # ALL ORGANISMS
			%fishindex = %{$refF};

			next;
		};


		# ----------------------------------------------------------------------------------------
		# Open the weight file and read in the weights in GroupWeights hash
		my %GroupWeights;  # group weights. Again, global and group is the same here.

		open WF, "Analysis_for_q-score\/teleost_fish_3R_filtered.txt" or die $!;
		foreach (<WF>){
			my @wghts = split "\t", $_;
			map {$_=~s/\t|\n//g} @wghts;
			#print "$wghts[1]\t$wghts[3]\n";
			$GroupWeights{$wghts[1]} = $wghts[3]

		}
		close (WF);


		# ----------------------------------------------------------------------------------------		
		# Now add the values for the current organism that could not be added in the subroutein. See NOTE in the subroutein.
		$fishindex{$organism}{'og'} = 7; # Index for outgroup for self organism, check!
		$fishindex{$organism}{'self'} = 15; # Index for self for self organism, check!
		

		# ----------------------------------------------------------------------------------------		
		# Calculate global/group averaged q-score
		# Q-score formula is: qscore^k = Exp[ Sum_{species i in group k or global} p^k_i log(species i qscore)] (see email by HI)
		# The weights will be from all-fish file

		my $sumog = 0; my $sumself = 0;
		#print "$line[2] - $line[3]\t GROUP\n";
		# Foreach organism, weighted sum group q-scores

		foreach my $organ (keys %fishindex){
						
			#print "$organ\t", $fishindex{$organ}{'og'},"\t", $fishindex{$organ}{'self'}, "\n";
			
			# Calculate for outgroups ---------------------			
			if ($line[$fishindex{$organ}{'og'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$fishindex{$organ}{'og'}], "\t", log($line[$fishindex{$organ}{'og'}]), "\t", $GroupWeights{$organ}, "\t", (log($line[$fishindex{$organ}{'og'}]) * $GroupWeights{$organ}), "\n";
				if ($line[$fishindex{$organ}{'og'}] <= 0){$line[$fishindex{$organ}{'og'}] = 1e-16}; # Assign a very low probability value if it's 0 or negative
				$sumog = $sumog + (log($line[$fishindex{$organ}{'og'}]) * $GroupWeights{$organ})
			}

			# Calculate for self ---------------------		
			if ($line[$fishindex{$organ}{'self'}] ne ''){ # If the value exists for this outgroup
				#print "$organ\t", $line[$fishindex{$organ}{'self'}], "\t", log($line[$fishindex{$organ}{'self'}]), "\t", $GroupWeights{$organ}, "\t", (log($line[$fishindex{$organ}{'self'}]) * $GroupWeights{$organ}), "\n";

				if ($line[$fishindex{$organ}{'self'}] <= 0){$line[$fishindex{$organ}{'self'}] = 1e-16;} # Assign a very low probability value if it's 0 or negative
				$sumself = $sumself + (log($line[$fishindex{$organ}{'self'}]) * $GroupWeights{$organ})
			}
			
			
			
		}

		print OUT join ("\t", @line[0..7]),"\t";
		print OUT "$line[15]\t";
		print OUT exp($sumog), "\t";
		print OUT exp($sumself), "\t";
		print OUT join ("\t", @line[8..14]),"\t";
		print OUT join ("\t", @line[16..$#line]),"\n";


=cut
		# Declare to hold array for the 2 P values
		my @pOg; my @pSelf;
	
		# Push the two p-values for the 1st organism, as there are extra columns for this one. ONLY IF THEY ARE NOT BLANK
		if ($line[7] ne ''){push @pOg, $line[7];}
		if ($line[15] ne ''){push @pSelf, $line[15];}
	
		# push the outgroup P-value in hash. ONLY IF THIS IS NOT BLANK
		for (my $ogStart = 17; $ogStart < scalar (@line); $ogStart = $ogStart + 3){
			#print "$ogStart\t";
			if ($line[$ogStart] ne ''){push @pOg, $line[$ogStart];}
		}
		#print "Og: @pOg\n";
	
		# calculate geometric mean, if there is at least one P
		my $ogGM = '';
		if ((scalar @pOg) > 0){
		
			foreach (@pOg) {if ($_ <= 0){$_ = 1e-16};} # Assign a very low probability value if it's 0
		
			my $ogStat = Statistics::Descriptive::Full->new();
			$ogStat->add_data(@pOg);
			$ogGM = $ogStat->geometric_mean();
		}
		#print "*$ogGM*\n";
	
		# push the SELF P-value in hash. ONLY IF THIS IS NOT BLANK
		for (my $selfStart = 18; $selfStart < scalar (@line); $selfStart = $selfStart + 3){
			#print "$selfStart\t";
			if ($line[$selfStart] ne ''){push @pSelf, $line[$selfStart];}
		}
		#print "Self: @pSelf\t";
	
		# calculate geometric mean, if there is at least one P
		my $selfGM = '';
		if ((scalar @pSelf) > 0){

			foreach (@pOg) {if ($_ <= 0){$_ = 1e-16};} # Assign a very low probability value if it's 0

			my $selfStat = Statistics::Descriptive::Full->new();
			$selfStat->add_data(@pSelf);
			$selfGM = $selfStat->geometric_mean();
		}
		#print "-> $selfGM <-\n";
	
		print OUT join ("\t", @line[0..7]),"\t";
		print OUT "$line[15]\t";
		print OUT "$ogGM\t$selfGM\t";
		print OUT join ("\t", @line[8..14]),"\t";
		print OUT join ("\t", @line[16..$#line]),"\n";
=cut	
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

