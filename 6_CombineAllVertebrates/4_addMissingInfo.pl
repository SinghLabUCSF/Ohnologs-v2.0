# This is to add the gene name and duplication type for some of the genes that are added later
# I can add filters latr here
local $| = 1;

use strict;
use warnings;

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 7: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus' (exclude 'oanatinus')

# --------------------------------------------------------------------------------------------------------------------

my $wgd = '_2R'; # _3R and _2R
#my $paralog_path = '../2_Paralogs/No-Sarco-Neo/4_filter_multi_copy_paralogs_No-Sarco-Neo'; # Make sure that the paralog files are correct --> with or without sarco-neo
my $paralog_path = '../2_Paralogs/4_filter_multi_copy_paralogs'; # Make sure that the paralog files are correct --> with or without sarco-neo

my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

# If the WGD is _3R -> run it only for the fish
if ($wgd eq '_3R'){
	@allorgs = @fish;
}

# --------------------------------------------------------------------------------------------------------------------
foreach my $organism (@allorgs){ # Process each organism one by one **************************************************

#	if ($organism eq 'drerio'){ # tst for 1 organism
	
	print "> Processing $organism $wgd\n";
	# --------------------------------------------------------------------------------------------------------------------
	# Read paralog file for the base organism having reconciled duplication time
	my $paralogFile = "$organism\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt"; # All paralogs in the base organism
	print "  Paralogs: $paralogFile\n";
	
	open PARALOGS, "$paralog_path/$paralogFile" or die $!;

	my %Paralogs;

	while (<PARALOGS>){

		my @line = split "\t", $_;
		map {$_=~s/\t|\n//g} @line;
		#print "$line[0]\t$line[2]\n";
	
		# Both sided hash so I have to check only in one direction later	
		$Paralogs{"$line[0]\t$line[1]"} = $line[5]; # Id1 Id2 as key and duplication node as value is the gene type
	}
	close (PARALOGS);
	#print join "\n", @{$Paralogs{"ENSG00000051596	ENSG00000101849"}},"\n";

	# --------------------------------------------------------------------------------------------------------------------
	
	# Open outfile
	my $infile = "$organism/$organism\_allVertebrates$wgd\_withGeoMean_filtered.txt";
	my $outfile = "$organism/$organism\_allVertebrates$wgd\_withGeoMean_filtered_final.txt";
	
	getGeneName($infile, $outfile, \%Paralogs);
	
#	} # end test for 1 organism
} # close the organism array ****************************************************************************************************



sub getGeneName {
	
	my $in = shift;
	my $out = shift;
	my $ref = shift;
	my %paralogs = %$ref;
	
	print "  Outfile: $out\n";
	open OUT, ">$out" or die $!;
		
	# Open consensus file
	print "  Infile: $in\n\n";
	open CONSENSUS, "$in" or die $!;
	
	while (<CONSENSUS>){
		
		if ($_ !~/^ENS.+/g){print OUT $_; next;} # Print header and go to next
	
		my @line = split "\t", $_;
		map {$_=~s/\t|\n//g} @line;
		
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t";
		
		if ($line[4] eq ''){ # If the gene type is missing -> add it
				
				if (exists $paralogs{"$line[0]\t$line[1]"}){print OUT $paralogs{"$line[0]\t$line[1]"},"\t";}
				elsif (exists $paralogs{"$line[1]\t$line[0]"}){print OUT $paralogs{"$line[1]\t$line[0]"},"\t";}
				else {print OUT "$line[4]\t";}
		}
		else {print OUT "$line[4]\t";}
	
		print OUT join "\t", @line[5..$#line];
		print OUT "\n";
	}
	close (CONSENSUS);
}
