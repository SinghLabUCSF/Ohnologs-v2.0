# This is same as combining the outgroups but I have to just print the multiplication of
# probability for all the outgroups instead of geometric mean.
# I select the maximum parobability - only P1(>=k) - from the pair for geometric mean.
# And so no need to bother about the A-B or B-A orientation for the pairs.

# ********* MAKE SURE THAT THE FILE NAMES ARE OKAY AND INCLUDE ONLY THE FILES YOU WANT TO COMBINE *************

use strict;
use warnings;

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 7: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus' (exclude 'oanatinus')

############ WGD for fish #############
my $wgdtype = '_3R';
# for tetrapods this will be blank '', use '_2R' or '_3R' only for the fish
#--------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri');
my @outgroups2r = ('bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster');
my @outgroups3r = ('lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus');

my %combinedHash; # Contains ids as keys and the all the probabilities string as value for that pair.
my %outgroups;    # this hash will hold the names and values to be printed

# Decide coorect organisms for WGD
if ($wgdtype eq '_3R'){
	@allorgs = @fish;
}
elsif ($wgdtype eq '' || $wgdtype eq '_2R'){
	# Outgroups is outgroups
}
else {die "Outgroup name is not correct!\n";}

# For each organism ********************************************************************************************
foreach my $polyploidName (@allorgs){
	
	# Decide wgd ----------------------
	my $wgd = ''; # for tetrapods
	
	if ($polyploidName ~~ @fish){
		$wgd = $wgdtype;
	}

	# Decide outgroup array ----------------------
	my @outgroups;
	if ($wgdtype eq '' || $wgdtype eq '_2R'){@outgroups = @outgroups2r}
	elsif ($wgdtype eq '_3R'){
		if ($polyploidName ~~ @fish){
			@outgroups = @outgroups3r;
		}
		else {die "_3R has not been found in tetrapods yet! Change it to _2R.\n"}
	}
	else {die 'Unidentified WGD string!';}
		

if ($polyploidName eq 'nfurzeri'){ # test for one organism

	open OUT, ">$polyploidName\/$polyploidName\_Ohnologs_CombinedFromAllOutgroups$wgd\.txt" or die $!;
	
	foreach (@outgroups){
		
		my $ogfile = "$polyploidName\/$_\-$polyploidName\_CombinedPairs_AllWindows$wgd\.txt";

		print "$ogfile\n";
		readOhnoPairFile($ogfile, $polyploidName); # Read the ohno pair file
	}
	print "\n";
#	print "\n*** WARNING ***\nCheck the file names above to make sure it includes all outgroups and no extra/missing file!\n";

	print OUT "Ohno1\tOhno2\tOutgroup Support\tMultiplication for P1(>=k)\t";
	#print scalar keys %outgroups;
	print OUT join ("\t", sort keys %outgroups), "\n";

	# Process the combined hash - It has a pair as key and all the lines from all outgroups in a value array
	foreach (keys %combinedHash){
	
		print OUT "$_\t";					# Print the ids for the pairs (key)
		print OUT scalar @{$combinedHash{$_}},"\t";		# Print the number of outgroups (from value array)
#		print @{$combinedHash{$_}},"\n------------\n";
	
		# declare and pupulate a probability array having max p value for each gene in the pair
		my @prob;

		foreach (@{$combinedHash{$_}}){ # For each of the outgroup
		
			my @pLine = split "\t", $_;
			map {$_=~s/\n//g} @pLine;
		
			# Select maximum out of P1(>=k) for the two genes and print it for each outgroup
			if ($pLine[4] > $pLine[10]){
				#print "$pLine[4]\t";
				push @prob, $pLine[4];
				$outgroups{$pLine[-1]} = $pLine[4]; # $pLine[-1] is the organism name that i appended at the end
			}
			else {
				#print  "$pLine[10]\t";
				push @prob, $pLine[10];
				$outgroups{$pLine[-1]} = $pLine[10];
			}
		}
		# multiply the probabilities
		my $mult = 1;
		map {$mult = $mult * $_;} @prob;
		print OUT "$mult\t";
	
		# Print max P1>=k for all outgroups
		my $counter = 1;
		foreach (sort keys %outgroups){
		
			# This condition is just to avoid printig the extra \t at the end
			if ($counter < (scalar keys %outgroups)){print OUT "$outgroups{$_}\t";}
			else {print OUT "$outgroups{$_}";}
		
			$outgroups{$_} = '';
			$counter++;
		}
		print OUT "\n";
	
	}
	%combinedHash = (); # Empty the hash ********** VERY IMPORTANT*************
	%outgroups = ();

} # end test for 1 organims


} # End organism loop *********************************************************************************************		
	






sub readOhnoPairFile {
	
	my $name = shift;
	my $polyploidName = shift;
	open FH, "$name" or die $!;
	
	my ($og, $rest) = split '-', $name, 2;
	$og =~s/$polyploidName\///g; # extract the outgroup name
	#print "$og\n";
	$outgroups{$og} = ''; # Initialize the outgroup hash with names
	
	# For each outgroup, here I push the entire line in the combinedHash, and then append the organism name at the end
	# I will process this hash at the end by multiplying these
	foreach (<FH>){
		
		#print "$_";
		my @line = split "\t", $_; # JUST BE SURE THAT THERE IS AN EXTRA TAB AT THE END OF SOME OF THE LINES. IT DOESN"T MATTER AT THE END AS I AM APPENDING THE ORGANISM AME AND GETTING IT BY INDEX -1. BUT BE CAUTIOUS
		
		# These conditions will check the duplicate pairs
		if (exists $combinedHash{$line[0]."\t".$line[1]}){
		
			push @{$combinedHash{$line[0]."\t".$line[1]}}, $_."\t$og";
		}
		if (exists $combinedHash{$line[1]."\t".$line[0]}){ # Since I just select the maximum probability for anchor, I don't need to correct the orientation like in combine P for windows program.
				
			push @{$combinedHash{$line[1]."\t".$line[0]}}, $_."\t$og";
		}
		if ((not exists $combinedHash{$line[0]."\t".$line[1]}) && (not exists $combinedHash{$line[1]."\t".$line[0]})) {
			push @{$combinedHash{$line[0]."\t".$line[1]}}, $_."\t$og";
		}
	}	
}

