# This is to further refine families
#
# 1. It will check if any of the ohnolog partners (separated by tab) have any genes that are recent SSD. If yes, then merge them as a single family.
# 2. For the ohnolog regions that have multiple partners at different genomic locations, it gets the one with highest synteny support, and puts the rest in a bracket.
#
# WE REMOVE THE BRACKETS FROM THE FAMILIES ON THE SERVER TO KEEP THEM SIMPLE.
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

print "> Processing $organism criteria $criteria $wgd\n";
# --------------------------------------------------------------------------------------------------------------------

my %Explored;
my @family;

# -------------- INPUT & OUTPUT FILES -------------------------------------------------------
my $allPCfile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$organism.'_Ens84.txt';
#my $allPCfile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$organism.'_NCBI.txt'; # For Nfurzeri

#my $paralogFile = "..\/2_Paralogs\/No-Sarco-Neo\/4_filter_multi_copy_paralogs_No-Sarco-Neo\/".$organism.'_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt';
my $paralogFile = "..\/2_Paralogs\/4_filter_multi_copy_paralogs\/".$organism.'_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt';
#my $paralogFile = '../OhnologsForNfur/2_Paralogs/2_Nfur_paralog_annotation/8_final_nfur_duplicates/Final_nfurzeri_paralogs.txt'; # For Nfurzeri 


my $selfSyntenyFile = "..\/5_Run_synteny\/3_SyntenyOutputFiles/$organism\/self$wgd\/$organism\_500-$organism\_OhnoCandidates_Window-P500R500_Orthologs-2.txt";
my $familyFile = "$organism/Families_Criteria-[$criteria]_$organism"."_Processed$wgd\.txt";
my $outFile = "$organism/Families_Criteria-[$criteria]_$organism\_ProcessedFinal"."$wgd\.txt";
my $logFile = "$organism/$organism\_$criteria\_$wgd.log";
# ----------------------------------------------------------------------------------

print "  All genes   : $allPCfile\n";
print "  Input       : $familyFile\n";
print "  Paralog     : $paralogFile\n";
print "  Self synteny: $selfSyntenyFile\n";
print "  Outfile     : $outFile\n";
print "  -> Warning: Double check that the correct nodes are in place (including or exclusing Sarcopterygii and Neopterygii).\n\n";

# open outfile 
open OUT, ">$outFile" or die $!;
open LOG, ">$logFile" or die $!;

# All PC gene file to get names ------------------------------------------------------------------------------------------------
open FH1, "$allPCfile" or die "$! $allPCfile";
my @allpc = <FH1>;
shift (@allpc); 
close (FH1);

my %Symbols;

foreach (@allpc){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Symbols{$line[0]} = $line[3];
	#$Symbols{$line[3]} = $line[0]; # This is to test for symbols
}

#print scalar keys %Symbols;
# ------------------------------------------------------------------------------------------------------------------------------

# Read the duplicatinon timing file
open FH2, "$paralogFile" or die $!;
my @duplicates = <FH2>;
#shift (@duplicates);
close (FH2);

my %Timing;

foreach (@duplicates){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# These are Ids of old or contemporary duplication nodes from 2R
	if ($wgd eq '' || $wgd eq '_2R'){
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
		    	
				if ((exists $Symbols{$line[0]}) && (exists $Symbols{$line[1]})){
					$Timing{$line[0]}{$line[1]} = $line[2];
				}
		}
	}

	if ($wgd eq '_3R'){ # If this is a 3R WGD, also includes a few more nodes as old node
		if ($line[2] ne '' && 
		    $line[2] ne 'Euteleostomi' && 
		    $line[2] ne 'Vertebrata' && 
		    $line[2] ne 'Coelomata' && 
		    $line[2] ne 'Opisthokonta' && 
		    $line[2] ne 'Chordata' && 
		    $line[2] ne 'Bilateria' && 
		    $line[2] ne 'Older' && 
		    $line[2] ne 'Clupeocephala' &&
		    $line[2] ne 'FishWGD' && # I can keep fishWGD. As this is not in the No SarcoNeo version, this will not affect things.
		    $line[2] ne 'Teleosts' && # Only for Nfurzeri
		    $line[2] ne 'Acanthomorphata'){
  
				if ((exists $Symbols{$line[0]}) && (exists $Symbols{$line[1]})){
					$Timing{$line[0]}{$line[1]} = $line[2];
				}
		}
	}
}
#print scalar keys %Timing;

# ------------------------------------------------------------------------------------------------------------------------------
# This will hold human paralogs on synteny and the color/statistical support
open FH3, "$selfSyntenyFile" or die $!;
my @humanWGD = <FH3>;
shift @humanWGD;
close (FH3);

my %humanWGD;

foreach (@humanWGD){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[5] == 1){
		
		#print "$line[6]\t$line[7]\n";
		#print "$Symbols{$line[6]}\t$Symbols{$line[7]}\t$line[4]\n";
		$humanWGD{$line[6]}{$line[7]} = $line[4];
		$humanWGD{$line[7]}{$line[6]} = $line[4];
	}
}

#print scalar keys %humanWGD;

# ------------------------------------------------------------------------------------------------------------------------------


# Read combined family file
open FH4, "$familyFile" or die $!;
my @outgroupFam = <FH4>;
close (FH4);
	
	
foreach (@outgroupFam){ # read each line of family file
	
	$_=~s/\t+$|\n//g;
	#print "$_\n\n";
	my @line = split "\t", $_; # @line holds ohno pairs with ssds separated with pipe or comma
	shift (@line);             # remove the first element of line that is the family size
	map {$_=~s/\n//g} @line;

	my %recentDuplicates; my %ohnoPairs; my %allOhno;
	
	# Foreach combination of ohno pairs
	for (my $i = 0; $i < scalar(@line); $i++){

		for (my $j = $i+1; $j < scalar(@line); $j++){
			
			#print "$line[$i]\t$line[$j]\n";
			
			my @array1 = split /[\|\,]/, $line[$i]; # separate them in 2 arrays having SSD genes also if any
			my @array2 = split /[\|\,]/, $line[$j];
			
			for (my $k = 0; $k < scalar(@array1); $k++){ # generate individual pairs by another nested loop
			
				for (my $l = 0; $l < scalar(@array2); $l++){
					
					#print "$array1[$k]\t$array2[$l]\n";
					my $id1 = $array1[$k]; my $id2 = $array2[$l];
					
					if ((exists $Timing{$id1}{$id2}) || (exists $Timing{$id2}{$id1})){ # see if they have been duplicated later than WGD i.e. young duplicates
							
						#print "$id1\t$id2\t$Timing{$id2}{$id1}\n";
						#print "$line[$i]\t$line[$j]\n";
						push @{$recentDuplicates{$line[$i]}}, $line[$j];  # Have a hash of late duplicates. THIS HAS SSD ALSO - NOT INDIVIDUAL IDS
						push @{$recentDuplicates{$line[$j]}}, $line[$i];
					}
					else { # If no then make 2 hashes - one having just the ohno ids and one having ohno pairs in both directions
						#print "* $line[$i]\t$line[$j]\n";
						$allOhno{$line[$i]} = $line[$j];  # THIS HAS SSD ALSO - NOT INDIVIDUAL IDS
						$allOhno{$line[$j]} = $line[$i];
						
						#print "* $id1\t$id2\n";
						$ohnoPairs{$id1}{$id2} = '';
						$ohnoPairs{$id2}{$id1} = '';
					}
				}
			}
		}
	}
# This will merge 2 families if any of the genes is a recent dupliacte
my %lateDuplicates;
foreach (keys %recentDuplicates){ # For each ohnolog gene
	
	if (not exists $Explored{$_}){                       # If it has not been explored already
		
		my @fam = getConnections(\%recentDuplicates, $_);   # get the family from this subroutein. Pass the current node and all relations
		#print "@fam\n";
		my $str = join '~', @fam;
		$lateDuplicates{$str} = '';
		@family = (); # Empty the family array as it's a global one
	}

}

	#print join "*\t", keys %lateDuplicates,"\n";
	#print join "+\t", keys %ohnoPairs,"\n";
	#print join "-\t", keys %allOhno,"\n";
	
	# In cases where there are multiple ohnolog partners on diferent regions, this gets the one with most statistical support from the self synteny file.
	my $processedLateDuplicates = getProperOhnolog(\%lateDuplicates, \%ohnoPairs); 
	#my %LateDuplicatesWithColor = %$ldRef;
	
	#my $processedLateDuplicates = getSyntenicLateDuplicate(\%LateDuplicatesWithColor);
	
	#chop $processedLateDuplicates; # remove comma at the end of processed family
	
	my $FinalProcessedFamily = '';	
	#print "$processedLateDuplicates" if ($processedLateDuplicates ne ''); # If their are no late duplicates - this string will be blank
	$FinalProcessedFamily = $FinalProcessedFamily.$processedLateDuplicates if ($processedLateDuplicates ne '');
	
	foreach (keys %allOhno){
		
		if (not exists $recentDuplicates{$_}){
			#print "\t$_";
			$FinalProcessedFamily = $FinalProcessedFamily."\t$_";
		}
	}
	#print "\n";
	
	my @size = split "\t", $FinalProcessedFamily;
	#print "*$size[0]*";
	if ($size[0] eq ''){shift @size}; # remove an element because of extra tab if there is any
	
	# Convert symbols to ids here
	foreach (@size){
		
		while ($_ =~/(ENS.{0,3}\d{11})/g){ # symbol start for Human
			my $k = $1;
			
			if (exists $Symbols{$k} && $Symbols{$k} ne ''){ # if exists in hash and symbol is not null then only convert
				#$_=~s/$1/$Symbols{$k}/g;
			}
		}
	}
	
	if (scalar (@size) == 1){
		
		print LOG scalar (@size),"\t";
		print LOG join ("\t", @size),"\n"; # print gene symbols
	}
	else {
		
		print OUT scalar (@size),"\t";
		print OUT join ("\t", @size),"\n"; # print gene symbols
	}
}

# This one gets the strength of synteny between groups of genes from the self synteny file of largest window
# This is used to decide the proper ohnolog, in case theer are more than one at different locations
sub getProperOhnolog {
	
	my $LDref = shift;
	my $OPref = shift;
	
	# Get late duplicates and ohno pair hash dereferenced
	my %LateDuplicates = %{$LDref};
	my %OhnoPairs = %{$OPref};
	my $finalFamly = '';
	
	foreach my $regions (keys %LateDuplicates){      # For each bunch of genes in late duplicates
		
		#print "$regions\n";
		my @cluster = split '~', $regions;
		my %LateDuplicatedRegions;
		
		foreach my $geneGroup (@cluster){
			
			#print "$geneGroup\n";
			my @ssd = split /[\|\,]/, $geneGroup;          # Get individual gene
			my $sum = 0; my $avg = 0; my $count = 0;

			foreach my $ssd (@ssd){  		   # For each of the gene in this gene group

	#			print "$geneGroup\t$ssd->\t";
	#			print join " ", keys %{$OhnoPairs{$ssd}},"\n";


				foreach (keys %{$OhnoPairs{$ssd}}){ # get it's ohno pairs

	#				print "* $geneGroup\t$ssd ->\t$_\n";

					if (exists $humanWGD{$ssd}{$_}){
	#					print "Group: $geneGroup\tGene: $ssd\tOhno: $_\tColor: $humanWGD{$ssd}{$_}\n";
						$sum = $sum + $humanWGD{$ssd}{$_};
						$count++;
					}
					else {
	#					print "->Error: $ssd or $_ do not exist in human paralogs\n";
					}
				}
			}
			if ($count !=0){$avg = $sum / $count;} # $count will be zero if none of the genes duplicated later exist in human paralog dataset
			#print "Sum = $sum\tCount: $count\tAverage = $avg\n";
			$LateDuplicates{$geneGroup} = $avg;
			$LateDuplicatedRegions{$geneGroup} = $avg;
			#print "*$LateDuplicates{$geneGroup}*\n";
		}
		#foreach (keys %LateDuplicatedRegions){
		#	print "$_\t$LateDuplicatedRegions{$_}\n";
		#}
		#print "---------------\n";
		my $processedLateDuplicates = getSyntenicLateDuplicate(\%LateDuplicatedRegions);
		chop $processedLateDuplicates;
		#print "$processedLateDuplicates\t";
		#print "---------------\n";
		$finalFamly = $finalFamly."\t".$processedLateDuplicates;
	}
	return ($finalFamly);
}



sub getSyntenicLateDuplicate {
	
	my $LDref = shift;
	my %ColoredLateDuplicates = %{$LDref};
	my $duplicateString = '';
	
	my @sorted = sort { $b <=> $a } values %ColoredLateDuplicates;  # Sort late duplicates based on their synteny in humans - in DESCENDING order
	#print join "\t", @sorted,"\n";
	
	foreach (keys %ColoredLateDuplicates){
		
		if ($ColoredLateDuplicates{$_} == $sorted[0]){ # If the region is highest synteny value (teher can be multiple regions with highest synteny value) -> print it outside the bracket
			#print "[$_]-\n";
			$duplicateString = $duplicateString."$_".',';
			#print "$duplicateString\n";
		}
		else { # Else print it in the bracket
			#print "$_-\n";
			$duplicateString = $duplicateString.'('."$_".')'.',';
			#print "*$duplicateString\n";
		}
	}
	#print "*** $duplicateString\n";
	
	return ($duplicateString) 
}


sub getConnections {
	
	my $graphRef = shift;
	my %Graph = %$graphRef; # Has all genes and ohnolog relations
	my $node = shift;       # The curent node
	#my @family; # Ohnolog family identified by depth first search
	
	#print "**$node** : ";
	
	if (not exists $Explored{$node}){push @family, $node;} # If the node has not been explored, start family using this node
	$Explored{$node} = 1;                                  # and mark it explored
	#print "@{$Graph{$_}}\n\n";
	foreach (@{$Graph{$_}}){    # For each of its ohnolog partners in graph. $_ here is $node. Check it properly later !!!
		
		my $nbor = $_;      # Take one partner or ohnolog neighbour at a time

		if ((not exists $Explored{$nbor})){ # if this one is not explored
			
			$Explored{$nbor} = 1; # mark it as explored
			
			push @family, $nbor;  # put it in family
	#		print "*\t$node\t$nbor\t$ogsup\n";
	#		print "Family: @family\n\n";
			getConnections(\%Graph, $nbor); # and recursively call the function to get all it's relations. It will keep on doing it as long as there are no nodes unexplored.
		}
	}
	#print "-------------\n< @family >\n------------\n";
	return @family; # Return this family
}
