# get the gene symbols
#
use strict;
use warnings;

# --------------------------------------------------------------------------------------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');

# get command line arguments
my $organism = $ARGV[0];
my $criteria = $ARGV[1];
my $wgd = $ARGV[2]; # _3R or _2R for fish

if ((scalar @ARGV) != 3){
	
	print "Script takes 3 arguments:
	1st: organism name e.g. hsapiens
	2nd criteria (A, C or E)
	3rd WGD type _2R or _3R
	
	Example usage: 1_DepthFirstSearchOhnolgFamilies_cl.pl drerio A _3R\n\n";
	print "Check parameters. Criteria can only be A, C or E. WGD can only be _2R and _3R and organism name must be one of the following\n";
	print "@allorgs\n\n";
	exit;
}

# Remove WGD if it's not fish
if ($organism ~~ @tetrapods){$wgd = '';} # Remove the WGD type for tetrapods

print "> Getting symbols for $organism criteria $criteria $wgd";

# -------------- INPUT & OUTPUT FILES ------------------------------------------------------------------------------
my $allpcFile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$organism.'_Ens84.txt';
my $inFile = "$organism/Families_Criteria-[$criteria]_$organism\_ProcessedFinal"."$wgd\.txt";
my $outFile = "$organism/Families_Criteria-[$criteria]_$organism\_ProcessedFinal-Symb"."$wgd\.txt";

# Open FILES ------------------------------------------------------------------------------------------------
open FH1, "$allpcFile" or die $!;
open FH, "$inFile" or die $!;
open OUT, ">$outFile" or die $!;
# ------------------------------------------------------------------------------------------------------------
my @allpc = <FH1>;
shift (@allpc);
close (FH1);

my %Symbols;

foreach (@allpc){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	$Symbols{$line[0]} = $line[3];
}


foreach (<FH>){
	
	my @line = split "\t", $_;
	# replace the brackets from the families
	$_=~s/\(|\)//g;
	
	$_ =~s/(ENS.{0,3}G\d{11})/$Symbols{$1}/g;
	
#	map {$_=~s/\n//g} @line;	
#	while ($_ =~/(ENS.{0,3}G\d{11})/g){ # symbol start
#		my $k = $1;
#		
#		if ((exists $Symbols{$k}) && ($Symbols{$k} ne '')){ # if exists in hash and symbol is not null then only convert
#			$_=~s/$1/$Symbols{$k}/g;
#		}
#	}

	print OUT $_;
}
print " ...... DONE\n\n";