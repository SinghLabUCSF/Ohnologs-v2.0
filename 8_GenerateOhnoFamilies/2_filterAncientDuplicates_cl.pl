# This is to discriminate between SSDs that are in the same window which duplicated after WGDs i.e. they are young duplicates.
# Recent/young duplicates would be separated by comma and potential old duplicates by a pipe. In this first one - I just mark the recent duplicates within a family. In the next one I will also check if there are potential SSDs among the different ohno families.
# Paralog dataset used is reconsiled paralogs from all versions.
# If a duplicate is not there in paralog dataset it would be considered as old duplicate.

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

print "> Processing $organism criteria $criteria $wgd\n";
print "  **** Check the nodes to mark as old or comtemporary in the script ****\n";

# -------------- INPUT & OUTPUT FILES -------------------------------------------------------------------------------
my $familyFile = "$organism/Families_Criteria-[$criteria]_$organism"."$wgd\.txt";

#my $paralogFile = "..\/2_Paralogs\/No-Sarco-Neo\/4_filter_multi_copy_paralogs_No-Sarco-Neo\/".$organism.'_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt';
#my $paralogFile = '../OhnologsForNfur/2_Paralogs/2_Nfur_paralog_annotation/8_final_nfur_duplicates/Final_nfurzeri_paralogs.txt'; # For Nfurzeri 
my $paralogFile = "..\/2_Paralogs\/4_filter_multi_copy_paralogs\/".$organism.'_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt';

my $outFile = "$organism/Families_Criteria-[$criteria]_$organism\_Processed$wgd\.txt";

print "  Input   : $familyFile\n";
print "  Paralog : $paralogFile\n";
print "  Outfile : $outFile\n\n";      

# OPEN FILES ---------------------------------------------------------------------------------------------------------
open FH, "$paralogFile" or die $!;
open FH2, "$familyFile" or die $!;
open OUT, ">$outFile" or die $!;
# --------------------------------------------------------------------------------------------------------------------

my %ssdFamily;
my %Explored;
my %ssdArray;
my %youngDuplicates;
#my %allDuplicates;

# Read reconsiled paralog file and make a two way hash for recent duplicates
foreach (<FH>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($wgd eq '' || $wgd eq '_2R'){
	
		# These are Ids of old or contemporary duplication nodes from 2R
		if ($line[2] ne '' && 
		    $line[2] ne 'Euteleostomi' && 
		    $line[2] ne 'Vertebrata' && 
		    $line[2] ne 'Coelomata' && # This is not in the current version -> but if it were, this would be an old node
		    $line[2] ne 'Opisthokonta' && 
		    $line[2] ne 'Chordata' && 
		    $line[2] ne 'Bilateria' && 
		    $line[2] ne 'Older' &&
		    $line[2] ne 'Sarcopterygii' && 
		    $line[2] ne 'Neopterygii'){

			$youngDuplicates{$line[0]}{$line[1]} = $line[2];
			$youngDuplicates{$line[1]}{$line[0]} = $line[2];
		}
	}
	if ($wgd eq '_3R'){
	
		# If this is a 3R WGD, also includes a few more nodes as old node
		if ($line[2] ne '' && 
		    $line[2] ne 'Euteleostomi' && 
		    $line[2] ne 'Vertebrata' && 
		    $line[2] ne 'Coelomata' && 
		    $line[2] ne 'Opisthokonta' && 
		    $line[2] ne 'Chordata' && 
		    $line[2] ne 'Bilateria' && 
		    $line[2] ne 'Older' && 
		    $line[2] ne 'Clupeocephala' &&
		    $line[2] ne 'FishWGD' &&
		    $line[2] ne 'Teleosts' &&
		    $line[2] ne 'Acanthomorphata'){
	
			$youngDuplicates{$line[0]}{$line[1]} = $line[2];
			$youngDuplicates{$line[1]}{$line[0]} = $line[2];
		}
	}
}

#print scalar keys %youngDuplicates;

# --------------------------------------------------------------------


foreach (<FH2>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	my $size = shift @line;
	
	# Final string to hold processed family
	my $finalString = '';
	
#	print join "\t", @line,"\n\n";
#	print "$size\t";
	
	$finalString = $finalString. "$size\t";
	
	foreach my $ohno (@line){ # For each ohnolog 
		
		my @ssd = split '\|', $ohno; # get the SSDs
		if (scalar @ssd > 1){ # if there are ssds
			
			# get all possible pairs by a nested hash and push it in a two way SSD hash having each gene and all its ssd partners
			for (my $i = 0; $i < scalar @ssd; $i++){
				
				for (my $j = $i+1; $j < scalar @ssd; $j++){
										
					push @{$ssdArray{$ssd[$i]}}, $ssd[$j];
					push @{$ssdArray{$ssd[$j]}}, $ssd[$i];
					#print "$ssd[$i]\t$ssd[$j]\t$allDuplicates{$ssd[$i]}{$ssd[$j]}\n";
				}
			}
		
		#print join ("\n", keys %ssdArray),"\n---\n";
		
		# Foreach ssd gene
		foreach (keys %ssdArray){
		
#			print "Start Node: $_\n";
			%ssdFamily = getConnections(\%ssdArray, $_); # Mark all recent duplications using depth first search and get recent duplicated ssd family
			
			if (scalar keys %ssdFamily > 1){ # If there are more than one gene
#				print join ('|', keys %ssdFamily), "_";
				$finalString = $finalString.join (',', keys %ssdFamily)."|";
			}
			else {
				if ((not exists $Explored{$_})){ # else if there is just one gene which is not explored, its late duplicate
#					print "$_"."_";
					$finalString = $finalString."$_|";
				}
			}
			%ssdFamily = ();
		}
				
		%ssdFamily = ();
		%Explored = ();
		%ssdArray = ();
#		print "\t";
		$finalString = $finalString."\t";
		}
		else {
#			print "$ohno\t";
			$finalString = $finalString."$ohno\t";
		}
	}
#	print "\n";
	$finalString = $finalString."\n";
	
	$finalString =~s/\|\t/\t/g;
	$finalString =~s/\t\n/\n/g;
	
	print OUT "$finalString";
}


# For each ssd gene, it will get recursively all the duplicated partners of all genes, until no new partners are found
# When there are no new partners it will return the family which has all recent duplicates
# I print that fimily and separate it by others with pipe, and then go to next gene
sub getConnections {
	
	my $graphRef = shift;
	my %Graph = %$graphRef; # Has all genes and ohnolog relations
	my $node = shift;       # The curent node
	#my @family; # Ohnolog family identified by depth first search
	
#	print "\tCurrent Node: $node\n";
	
	#if (not exists $Explored{$node}){push @ssdFamily, $node;} # If the node has not been explored, start family using this node
	#$Explored{$node} = 1;                                     # and mark it explored
#	print "@{$Graph{$_}}\n";
	
	foreach (@{$Graph{$_}}){    # For each of its ohnolog partners in graph. $_ here is $node. Check it properly later !!!
		
		my $nbor = $_;      # Take one partner or ohnolog neighbour at a time
#		print "$node   Neighbour   $nbor\n";
		
		if ((not exists $Explored{$nbor}) && (exists $youngDuplicates{$nbor}{$node})){ # if this one is not explored
			
			$Explored{$nbor} = 1; # mark it as explored
			$Explored{$node} = 1; # mark it as explored
#			print "* $nbor $node\n";
			$ssdFamily{$nbor} = '';  # put it in family
			$ssdFamily{$node} = '';
			getConnections(\%ssdArray, $nbor); # and recursively call the function to get all it's relations. It will keep on doing it as long as there are no nodes unexplored.
		}
	}
	return %ssdFamily; # Return this family
}

