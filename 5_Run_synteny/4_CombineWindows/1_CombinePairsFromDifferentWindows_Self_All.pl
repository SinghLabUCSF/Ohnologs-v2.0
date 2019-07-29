# Read ohno pair files for different windows as input and combine probability for anchors from different windows.
# 
# The difference in self combination is that I have pairs directly from blocks so
# in one line there is just one set of probability for both the pairs.
# Secondly, each pair would be twice because same region would be compared with self 2 times.
# So I calculata and take a geometric mean for both A->B and B->A comparison.
#

use strict;
use warnings;  # There will be a warning for tetrapods for $wgd, but it can be ignored

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 7: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus' (exclude 'oanatinus')

my %combinedHash; # Contains ids as keys and the all the probabilities string as value for that pair.


############ WGD for fish #############
my $wgdtype = '_3R';
# for tetrapods this will be blank '', use '_2R' or '_3R' only for the fish
#--------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri');

# If the WGD is _3R -> run it only for the fish
if ($wgdtype eq '_3R'){
	@allorgs = @fish;
}

# For each organism ********************************************************************************************
foreach my $polyploidName (@allorgs){
	
	my $wgd = ''; # for tetrapods
	if ($polyploidName ~~ @fish){
		$wgd = $wgdtype;
	}

if ($polyploidName eq 'nfurzeri'){ # Test for one organism

	print "> Processing: $polyploidName\t$wgd\n";

	open OUT, ">$polyploidName\/$polyploidName-$polyploidName\_CombinedPairs_AllWindows$wgd\.txt" or die $!;

	# Foreach ohno pair file for self
	foreach (<..\/3_SyntenyOutputFiles\/$polyploidName\/self$wgd\/$polyploidName*-$polyploidName\_OhnoPairs*>){ # Folder name and path

		print "  $_\n";
	
		# Read the ohno pair file 
		readOhnoPairFile($_);
		#last();
	}


	foreach (keys %combinedHash){
	
		if (scalar @{$combinedHash{$_}} > 1){ # If the gene exists in multiple windows - there will be more lines
		
			print OUT "$_\t"; # Print the ids
			#print scalar @{$combinedHash{$_}},"\t";
			#print @{$combinedHash{$_}};
		
			# Declare a loop of 12 iterations from 2 to 13. Because i have to calculate geometric mean of 12 probabilities - 6 for each gene in pair
			# The goal is to calculate geo mean for respective probability for all the windows
			for (my $i = 2; $i <=7; $i++){
			
				my @prob; # this array will hold a particular probability one at a time, for each window
				
				foreach (@{$combinedHash{$_}}){ # Now get the probability starting from Pchr for gene 1 in pair
			
					my @pLine = split "\t", $_;
					map {$_=~s/\n//g} @pLine;

					push @prob, "$pLine[$i]"; # push probability in @prob for geomean calculation
					#print "$pLine[$i]\t";	
				}
				#calculate geomean
				my $geoMean = geoMean(\@prob);
			
				print OUT "$geoMean\t";
				#print "\n";
			}
			print OUT "\n";
		}
		else { # If the gene exists only in one window
			#print OUT join "\t", @{$combinedHash{$_}};
		}
	}

	# close the output filehandle and call the program to filter one-way relationship
	#print "\nFiltering two-way relations by taking maximum\n";
	close (OUT);
	print `perl FilterTwoWayRelationships_CombinedPairs_Self.pl $polyploidName $wgd`;
	
	%combinedHash = (); # Empty the hash ********** VERY IMPORTANT*************

}# End test for one organism	

} # End the organism array *************************************************************************************




sub readOhnoPairFile {
	
	my $name = shift;
	open FH, "$name" or die $!;	
	#open FH, 'CionaInt-Human_OhnoPairs_Window-P100R100_Orthologs-2.txt' or die $!;
	#my %hash;
	
	foreach (<FH>){
		
		#print "$_";
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
		
		push @{$combinedHash{$line[0]."\t".$line[1]}}, $_;
=cut		
		# These conditions will check the duplicate pairs --- NOT NEEDED FOR SELF
		if (exists $combinedHash{$line[0]."\t".$line[1]}){
		
			push @{$combinedHash{$line[0]."\t".$line[1]}}, $_;
		}
		if (exists $combinedHash{$line[1]."\t".$line[0]}){ # If pair exists in reverse order - I also exchange the order of probability for easy calculation later on.
			
			#my $exchange = "$line[0]\t$line[1]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\n";
			#print "$exchange\n";
			#print "$_\n";
			push @{$combinedHash{$line[1]."\t".$line[0]}}, $_;
		}
		if ((not exists $combinedHash{$line[0]."\t".$line[1]}) && (not exists $combinedHash{$line[1]."\t".$line[0]})) {
			push @{$combinedHash{$line[0]."\t".$line[1]}}, $_;
		}
=cut
	}
	#print scalar keys %combinedHash, "\n";
}


sub geoMean {
	
	my $ref = shift;
	my @probabilities = @$ref; # array to hold the probability values
	
	# Remove blank values where P is not there for any window. If the ortholog pair does not belong to anchor.
	@probabilities = grep { $_ ne '' } @probabilities;
	#print join "\n", @probabilities;
	
	# Because of rounding error and some very low probability values some probabilities may be negative.
	# In that case I assign a very low arbritary probability of 1e-16. Basically these are very low probability of random occurance, the actual value doesn't matter much below a reasonably threshold.
	foreach (@probabilities){if ($_ <=0){$_ = 1e-16}};
	#print join "\n", @probabilities;
	
	# To take geometric mean calculate multiplication of all Ps
	my $pMult = 1; my $geoMean = 0;
	#print "\n->$pMult\n";
	
	map { $pMult = $pMult * $_} @probabilities;
	
	# Calculate geometric mean only if ortholog pair belongs to window for at least one window.	
	if (scalar @probabilities > 0){
		$geoMean = $pMult ** (1/(scalar @probabilities));
	}

	# If the probability is 0 assign a very low value of 1e-16. There are negligibal chances that this is a false positive, so a low value is assigned for ease of calculations.
	if ($geoMean == 0){$geoMean = 1e-16;}	
	return ($geoMean);
}
