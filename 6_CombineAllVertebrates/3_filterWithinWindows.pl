# Due to combining all vertebrates some candidates may be within the minimum window of 100 genes
# Here I filter those candidates
use strict;
use warnings;

# --------------------------------------------------------------------------------------------------------------------

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 7: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus' (exclude 'oanatinus')

my $wgd = '_2R'; # WGD type: for tetrapods '_2R' or for fish '_3R'
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

# If the WGD is _3R -> run it only for the fish
if ($wgd eq '_3R'){
	@allorgs = @fish;
}

# --------------------------------------------------------------------------------------------------------------------
foreach my $organism (@allorgs){ # Process each organism one by one ***************************************************************
	
	my $subGenomeFile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$organism.'_Ens84.txt';
	my $pairFile = "$organism/$organism\_allVertebrates$wgd\_withGeoMean.txt";
	my $outfile = "$organism/$organism\_allVertebrates$wgd\_withGeoMean_filtered.txt";

	open OUT, ">$outfile" or die $!;
	
	print "> $organism\n  $pairFile\n  $outfile\n  $subGenomeFile\n\n";

	# Read organism genome file and make hash having Id and chromosome position ---------------
	open WGDGENOME, $subGenomeFile or die $!;
	my @wgdGenome = <WGDGENOME>; 
	close (WGDGENOME);

	my %Positions;

	my $wgdGenomeCount = 1;
	for (my $i = 0; $i < scalar (@wgdGenome); $i++){

		# Read two lines of the file at a time
		my @currentLine = split "\t", $wgdGenome[$i];		
		my @previousLine = split "\t", $wgdGenome[$i-1];
		map {$_=~s/\n|^\s+|\s+$//g} @currentLine;
		map {$_=~s/\n|^\s+|\s+$//g} @previousLine;

		#print "$currentLine[2]\t$previousLine[2]\t$wgdGenomeCount\n";

		# If the chromosome number is different then reset the count to 1
		if ($currentLine[2] ne $previousLine[2]){
			$wgdGenomeCount = 1;
		}

		$Positions{$currentLine[0]} = [$currentLine[2], $wgdGenomeCount]; # hash holding Id and its chr, position
		$wgdGenomeCount++;
	}

	#print "@{$Positions{'ENSGACG00000012007'}}\t-\t";
	#print "@{$Positions{'ENSGACG00000010551'}}\n";

	# -------------------
	open PAIR, $pairFile or die $!;
	my @pairs = <PAIR>;
	my $header = shift @pairs;

	print OUT $header;
	foreach (@pairs){
	
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
	
		# If both the genes are on same chromosome
		if (${$Positions{$line[0]}}[0] eq ${$Positions{$line[1]}}[0]){
		
			my @positions = (${$Positions{$line[0]}}[1], ${$Positions{$line[1]}}[1]);
			@positions = sort {$a <=> $b} @positions;
		 
			# print only if the two genes are more than 100 genes apart
			if ($positions[1] > $positions[0]+100){
				#print "@positions\n";
				#print "$line[0]\t${$Positions{$line[0]}}[0] - ${$Positions{$line[0]}}[1]\t";
				#print "$line[1]\t${$Positions{$line[1]}}[0] - ${$Positions{$line[1]}}[1]\n";
				print OUT "$_";
			}
			else { # pairs that lie witing the windows and should be removed
				print "$line[0]\t${$Positions{$line[0]}}[0] - ${$Positions{$line[0]}}[1]\t";
				print "$line[1]\t${$Positions{$line[1]}}[0] - ${$Positions{$line[1]}}[1]\n";		
			}
		}
		else { # else print it
			print OUT "$_";
		}
	}
}
