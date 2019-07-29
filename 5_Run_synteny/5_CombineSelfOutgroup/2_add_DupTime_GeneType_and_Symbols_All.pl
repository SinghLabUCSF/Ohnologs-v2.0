# Add gene symbols and type from "All Genes" file and duplication time from paralog file

use strict;
use warnings;

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'

############ WGD for fish #############
my $wgdtype = '_2R';
# for tetrapods this will be blank '', use '_2R' or '_3R' only for the fish

#my $paralog_file = '../../2_Paralogs/No-Sarco-Neo/4_filter_multi_copy_paralogs_No-Sarco-Neo'; # For No-sarco-neo nodes
my $paralog_file = '../../2_Paralogs/4_filter_multi_copy_paralogs';

#-------------------------------------- Removed  'trubripes'
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');

# If the WGD is _3R -> run it only for the fish
if ($wgdtype eq '_3R'){
	@allorgs = @fish;
}

# For each organism ********************************************************************************************
foreach my $organism (@allorgs){ 

	# Decide wgd ----------------------
	my $wgd = ''; # for tetrapods
	
	if ($organism ~~ @fish){
		$wgd = $wgdtype;
	}
		
if ($organism eq 'nfurzeri'){ # Test for 1 org

	print "Getting degtails for $organism for WGD: $wgdtype\n";
	#print "** Warning: Remember I have multiple duplication time, so make sure the correct file path is given in the script.\n\n";

	open OUT, ">$organism\/$organism"."Ohno_Self+Outgp$wgd\_withDetails.txt" or die $!; # Outfile
	open DT, "$paralog_file/$organism\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt" or die $!; # Read the duplication time from filtered file
	open ALL, "../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$organism\_Ens84.txt" or die $!; # All genes and symbol file
	open OHNO, "$organism\/$organism"."Ohno_Self+Outgp$wgd\.txt" or die "$organism $!"; # Ohno file to add duplication time
	
	# Read the duplication time from filtered file -------------------------------------------------------------------------------------
	my %DupTime;

	foreach (<DT>){
	
		my @line = split "\t", $_;
		map {$_ =~s/\n//g} @line;
	
		$DupTime{$line[0]}{$line[1]} = "$line[2]";
		$DupTime{$line[1]}{$line[0]} = "$line[2]";	
	}

	# All genes file ------------------------------------------------------------------------------------------------------------------
	my %GeneInfo;

	foreach (<ALL>){
	
		my @line = split "\t", $_;
		map {$_ =~s/\n//g} @line;
	
		$GeneInfo{$line[0]} = [$line[3], $line[7]]; # This will have the gene symbol and type
	}

	# Read ohno file to add the dup time ----------------------------------------------------------------------------------------------
	my @ohno = <OHNO>;
	my $header = shift @ohno;
	my @header = split "\t", $header;
	map {$_ =~ s/\n//g} @header;

	print OUT "$header[0]	$header[1]	";
	print OUT "Symbol1	Symbol2	Gene type	Duplication time	";
	print OUT join ("\t", @header[2..$#header]), "\n";

	foreach (@ohno){
	
		my @line = split "\t", $_;
		map {$_ =~s/\n//g} @line;
	
		if (exists $DupTime{$line[0]}{$line[1]}){
		
			print OUT "$line[0]	$line[1]	${$GeneInfo{$line[0]}}[0]	${$GeneInfo{$line[1]}}[0]	${$GeneInfo{$line[1]}}[1]	$DupTime{$line[0]}{$line[1]}	"; # print ohno symbol, gene symbol for both the gene and the type for only one of them as it will be same for both paralogs
			print OUT join ("\t", @line[2..$#line]), "\n";
		}
		else {
			print OUT "$line[0]	$line[1]	${$GeneInfo{$line[0]}}[0]	${$GeneInfo{$line[1]}}[0]	${$GeneInfo{$line[1]}}[1]		";
			print OUT join ("\t", @line[2..$#line]), "\n";
		}
	}

} # end test for 1 org
} # End loop for each organism ********************************************************************************************


print `date`;
