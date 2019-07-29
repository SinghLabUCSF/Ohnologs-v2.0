# As I combine ohnologs from all vertebrates, it can happen that some ohnolog pairs still are 
# within the smallest window (because they are outside the smallest window in other vertebrates).
# Such ohnologs are listed as size 1 family and I remove them here.
# This will generate a family with the same name.
use strict;
use warnings;

# --------------------------------------------------------------------------------------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

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
# --------------------------------------------------------------------------------------------------------------------

my $familyFile = "$organism\/Families_Criteria-[$criteria]_$organism"."$wgd\.txt";

print "> Processing $organism criteria $criteria $wgd\n";
print "  Family file: $familyFile\n";
print "  Removed families: ";

open FAM, "$familyFile" or die "$! $familyFile\n";
my @family = <FAM>;
close (FAM);

# open outfile with the same name
open OUT, ">$familyFile" or die "$!. Could not open $familyFile for writing\n";
my $famcount = 0;
foreach (@family){
	
	my ($size, $fam) = split ("\t", $_, 2);
	#print "$fam";
	if ($size == 1){
		print "  - $size\t$fam";
		$famcount++;
	}
	else {
		print OUT "$size\t$fam";	
	}
}
print "total_removed = $famcount\n\n";
