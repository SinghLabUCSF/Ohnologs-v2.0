# This is to reconcile the duplication node from different ensembl versions
# I consider 7 (then) latest versions v80 to 86
# I decide final version from majority rule.
# Run it one by one for each organism and make sure there are no errors and warnings.
#
use strict;
use warnings;

local $| = 1;

my %Symbols;
my %Duplicates;
my %AllNodes;

my $organism = 'fcatus';

open OUT, ">2_Combined_nodes\/$organism\_CombinedDuplicates_Ens80-86.txt" or die $!;
open NODES, ">2_Combined_nodes\/AllNodes_Ens80-86_$organism\.txt" or die $!;
print OUT "Id1	Id2	Symbol1	Symbol2	Reconciled node	v80	v81	v82	v83	v84*	v85	v86	Total	Older	Vertebrates	Younger\n";

my $fileCount = 0;

foreach (<1_Paralogs_from_7_latest_versions/$organism*.txt>) { # Check that the versions are proper, they dont seem proper
	
	print "$_\n";
	
	open FH,"$_" or die $!;
	my @file = <FH>;
	shift @file;
	close ($_);
		
	foreach (@file){

		my @line = split "\t" or die $!;
		map {$_=~s/\n//g} @line;
		
		#$Symbols{$line[0]} = $line[1];

		if ($line[0] ne '' && $line[1] ne '' &&  $line[2] ne '' && $line[0] ne 'NA' && $line[1] ne 'NA' &&  $line[2] ne 'NA'){ # Its duplicate timing and id should not be null or NA
			
			if ((not exists $Duplicates{$line[0]."\t".$line[1]}) && (not exists $Duplicates{$line[1]."\t".$line[0]})){ # check if tha pair exists in any direction
				$Duplicates{$line[0]."\t".$line[1]} = ['', '', '', '', '', '', '']; # initialize the hash having an array of 7 for each ensembl version from 80 to 86
				${$Duplicates{$line[0]."\t".$line[1]}}[$fileCount] = $line[2];
			}
			else {
				# Here if a pair exists in both directions, the node will be replaced with the one that appear later. Although it shouldn't matter as
				# in the same version a pair cannot be associated with 2 nodes. Both must be same.
				if (exists $Duplicates{$line[0]."\t".$line[1]}){${$Duplicates{$line[0]."\t".$line[1]}}[$fileCount] = $line[2];} # push the value in the appropriate array
				if (exists $Duplicates{$line[1]."\t".$line[0]}){${$Duplicates{$line[1]."\t".$line[0]}}[$fileCount] = $line[2];} # push the value in the appropriate array
			}
			
			# push the node in nodes array
			$AllNodes{$line[2]} += 1;
		}
	}
	$fileCount++;		
}


# For all duplication nodes - find out the final class here based on majority rule
foreach (keys %Duplicates){
	
	# Print ids
	print OUT "$_\t";
	my ($id1, $id2) = split "\t", $_;
	
	# Print symbols if exists
	if (exists $Symbols{$id1}){print OUT "$Symbols{$id1}\t";}
	else {print OUT "\t";}
	
	if (exists $Symbols{$id2}){print OUT "$Symbols{$id2}\t";}
	else {print OUT "\t";}
	
	my $count = 0;
	my $nonBlank = 0;
	
	# initialize variables to decide if the variable is old, intermediate or recent
	my %category;
	$category{'old'} = 0;
	$category{'intermediate'} = 0;
	$category{'fishwgd'} = 0;
	$category{'young'} = 0;
	
	# Populate hash with counts
	foreach (@{$Duplicates{$_}}){
		
		if ($_ ne ''){
			
			if ($_ eq 'Opisthokonta' || $_ eq 'Bilateria' || $_ eq 'Coelomata'){ # Also check in the all nodes file that there is not any older node
				$category{'old'} += 1;
			}
			elsif ($_ eq 'Vertebrata' || $_ eq 'Euteleostomi' || $_ eq 'Chordata' || $_ eq 'Sarcopterygii' || $_ eq 'Neopterygii'){ # 
				$category{'intermediate'} += 1;
			}
			elsif ($_ eq 'Clupeocephala' || $_ eq 'Acanthomorphata'){
				$category{'fishwgd'} += 1;
			}
			else {
				$category{'young'} += 1;
			}
			$nonBlank++;
		}
	}
	
		
	# Decide for different versions if the duplicate can be classified properly or not - uncertain
	if ($nonBlank == 6 || $nonBlank == 7){
		
		if ($category{'old'}             >= 4){print OUT "Older than vertebrates";}
		elsif ($category{'intermediate'} >= 4){print OUT "Vertebrates";}
		elsif ($category{'fishwgd'}      >= 4){print OUT "FishWGD";}
		elsif ($category{'young'}        >= 4){print OUT "Younger than vertebrates";}
		else {print OUT "Uncertain";}
	}
	if ($nonBlank == 4 || $nonBlank == 5){
			
		if ($category{'old'}             >= 3){print OUT "Older than vertebrates";}
		elsif ($category{'intermediate'} >= 3){print OUT "Vertebrates";}
		elsif ($category{'fishwgd'}      >= 3){print OUT "FishWGD";}
		elsif ($category{'young'}        >= 3){print OUT "Younger than vertebrates";}
		else {print OUT "Uncertain";}
	}
	if ($nonBlank == 2 || $nonBlank == 3){
				
		if ($category{'old'}             >= 2){print OUT "Older than vertebrates";}
		elsif ($category{'intermediate'} >= 2){print OUT "Vertebrates";}
		elsif ($category{'fishwgd'}      >= 2){print OUT "FishWGD";}
		elsif ($category{'young'}        >= 2){print OUT "Younger than vertebrates";}
		else {print OUT "Uncertain";}
	}
	if ($nonBlank == 1){
				
		if ($category{'old'}             >= 1){print OUT "Older than vertebrates";}
		elsif ($category{'intermediate'} >= 1){print OUT "Vertebrates";}
		elsif ($category{'fishwgd'}      >= 1){print OUT "FishWGD";}
		elsif ($category{'young'}        >= 1){print OUT "Younger than vertebrates";}
		else {print OUT "Uncertain";}
	}
	
	print OUT "\t";
	# Print the duplicate timing for different versions
	print OUT (join "\t", @{$Duplicates{$_}}),"\t";
	
	# print number of times a pair exists in different versions and if it is intermediate, old or young duplicate count
	print OUT "$nonBlank\t".$category{'old'}."\t".$category{'intermediate'}."\t".$category{'young'};
	
	print OUT "\n";
}


# print the nodes in the nodes file -- check this file properly and see if all the nodes are covered in this script
foreach (keys %AllNodes){
	
	print NODES "$_\t$AllNodes{$_}\n";
}

print `date`;
