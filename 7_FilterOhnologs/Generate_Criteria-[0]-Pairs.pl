# Filter ohnologs for diferent criteria
# No warnings allowed.
use strict;
use warnings;

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 8: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'oanatinus', 'sscrofa', 'cfamiliaris', 'ggallus'

# --------------------------------------------------------------------------------------------------------------------
my $wgdtype = '_2R'; # WGD type '_2R' or '_3R' 
my $score_type = 'group'; # 'group' or 'global'

# Columns with geo. mean of the final OUTGROUP and SELF q scores
my $og = 10000; # A number more than columns, so we get a warning
my $self = 10000;

# For global q score means across all vertebrates.
if ($score_type eq 'global'){
	$og = 9;
	$self = 10;
}

# For group-wise means --------------------------------------
if ($score_type eq 'group'){
	$og = 94;
	$self = 95; 
}
# --------------------------------------------------------------------------------------------------------------------

my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

# If the WGD is _3R -> run it only for the fish
if ($wgdtype eq '_3R'){
	@allorgs = @fish;
	
	# Reset for _3R
	$og = 9;
	$self = 10;
}


foreach my $organism (@allorgs){ # Process each organism one by one 
	
#	if ($organism eq 'hsapiens'){ # test for 1 org
	
	print "> Processing $organism $wgdtype\n";	
    if ($og == 9 && $self == 10){print "  Filtering based on *** GLOBAL *** q-score\n"}
    elsif ($og == 94 && $self == 95){print "  Filtering based on *** GROUP *** q-score\n"}
    else {print "Check the column number for q-score\n"; exit;}
    print "  Index for outgroup q-score =  $og & self q-score = $self. Check again (Colum = Index + 1)!\n";
	print "  INPUT: ../6_CombineAllVertebrates/$organism\/$organism\_allVertebrates$wgdtype\_withGeoMean_filtered_final.txt\n";

	open FH, "../6_CombineAllVertebrates/$organism\/$organism\_allVertebrates$wgdtype\_withGeoMean_filtered_final.txt" or die $!;
	
	# remove the WGD string for tetrapod output file to keep it compatible with the server and old code
	my $wgd = $wgdtype;
	if (($wgdtype eq '_2R') && ($organism ~~ @tetrapods)){
		$wgd = '';
	}
	
	# Open outfile	
	my $outfile = "$organism\/".$organism.'_Criteria-[0]-Pairs'.$wgd.'.txt';
	print "  OUTFILE: $outfile\n\n";	
	open OUT, ">$outfile" or die $!;
	
	while (<FH>){

		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;		
	
	   if ($line[0] !~/^ENS.+/){print OUT $_; next;} # Print and skip to next for header
		
		# Criteria [AA]
		if ($line[$og] ne '' && $line[$self] ne '' && $line[$og] < 0.001 && $line[$self] < 0.001){
			print OUT "$_";
		}
	}
#	} # end test for 1 org
}

