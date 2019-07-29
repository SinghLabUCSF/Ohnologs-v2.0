# Run all the organisms
use strict;
use warnings;
local $| = 1;

# --------------------------------------------------------------------------------------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio');

my @Wgd = ('_2R', '_3R');

#my @Criteria = ('E', 'C', 'A'); # Old criteria also had E
my @Criteria = ('0', 'C', 'A');


foreach my $Org (@allorgs){	

	foreach my $Crit (@Criteria){

		foreach my $WGD (@Wgd){
			
			if (($Org ~~ @tetrapods) && ($WGD eq '_3R')){next;} # No 3R WGD for tetrapods
			
			# Run the 4 scripts one by one
			
#			if ($Org eq 'hsapiens' && $Crit eq 'A'){ # test for 1 organism
#			if ($Org eq 'ptroglodytes'){ # test for 1 organism

			print "\n--------------------------------------------------------------------------------------------------------------------\n";
			print "PROCESSING: $Org\t$Crit\t$WGD\n";
			print "--------------------------------------------------------------------------------------------------------------------\n";
			
 			print `perl 1_DepthFirstSearchOhnolgFamilies_cl.pl $Org $Crit $WGD`;
 			print `perl 1b_filterSize1families_cl.pl $Org $Crit $WGD`;
 			print `perl 2_filterAncientDuplicates_cl.pl $Org $Crit $WGD`;
 			print `perl 3_filterRecentDuplicates_cl.pl $Org $Crit $WGD`;
			#print `perl 3b_getCliqueFamilies_cl.pl $Org $Crit $WGD`; # This is to get the clique families. Change the input name in the next script if you use it to get symbols.
			print `perl 4_getGeneSymols_cl.pl $Org $Crit $WGD`;
						
#			} # end test for 1 organism
		}		
	}
}
