# Combine the ohnologs from all vertebrates together and get their q-scores to take geometric mean at the next step
# 
# Here is how I do it -> Read ohno file of the base organism, and then for all other organisms.
#                     -> Get the orthologous ohno pairs in each organism one at a time
#                     -> If there is a single organism - print probability values
#                     -> If there are more than 1 ortholog pairs, get pair with best probabilities and print
#                     -> For all the ohnologs in the other organism, do above steps again
#
# MODIFIED FROM THE TETRAPOD SCRIPT
# --------------------------------------------------------------------------------------------------------------------

local $| = 1;

use strict;
use warnings;

# FISH: 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes', 'nfurzeri'
# TETRAPODS: 'acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus'
# Outgroups for 2R = 5: 'bfloridae', 'celegans', 'cintestinalis', 'csavignyi', 'dmelanogaster'
# Outgroups for 3R = 7: 'lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'sscrofa', 'cfamiliaris', 'ggallus' (exclude 'oanatinus')

# --------------------------------------------------------------------------------------------------------------------
# *************** CHECK THAT THE PATH IS FOR CORRECT FILES FOR DUPLICATEION NODES ***************
my $wgdtype = '_3R'; # WGD type. Use '_3R' for the fish for this script.
#my $paralog_path = '../2_Paralogs/No-Sarco-Neo/4_filter_multi_copy_paralogs_No-Sarco-Neo'; # Make sure that the paralog files are correct --> with or without sarco-neo
my $paralog_path = '../2_Paralogs/4_filter_multi_copy_paralogs'; # Make sure that the paralog files are correct --> with or without sarco-neo
my $ohnofilepath = '../5_Run_synteny/5_CombineSelfOutgroup';
my $totalOrganisms = 3; # This should be total organisms to be combined - 1 (as the starting organism will be the first one)
my @allorgs = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio'); # All organism array is fish
if ($wgdtype ne '_3R'){die "WGS type is not _3R\n";} # Check that WGS is proper

# --------------------------------------------------------------------------------------------------------------------
foreach my $organism (@allorgs){ # Process each organism one by one **************************************************
# $organism is the base organism 

print "*** Processing $organism\n";
print "*** WGD type: $wgdtype\n";

# --------------------------------------------------------------------------------------------------------------------

my $combinedOhnoFile = "$organism"."Ohno_Self+Outgp$wgdtype\_withDetails.txt"; # Decide the input file name for tetrapods or fish
my @organisms = ('tnigroviridis', 'gaculeatus', 'olatipes', 'drerio'); # Decide the organism names and order in which to add other vertebrates
my $allPCFile = "AllGenes_$organism\_Ens84.txt"; # all genes file for the base organism
my $paralogFile = "$organism\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt"; # All paralogs in the base organism
my $orgCount = 1; # Organism counts for proper printing

# --------------------------------------------------------------------------------------------------------------------
# All PC gene file for the start/base organism 
open ALLPC, "../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/$allPCFile" or die $!;

my %AllPC;
while (<ALLPC>){

	my @line = split "\t", $_;
	map {$_=~s/\t|\n//g} @line;
	
	$AllPC{$line[0]} = $line[3];
}
close (ALLPC);
#print scalar keys %AllPC;

# --------------------------------------------------------------------------------------------------------------------
# Read paralog file for the base organism having reconciled duplication time

open PARALOGS, "$paralog_path/$paralogFile" or die $!;

my %Paralogs;

while (<PARALOGS>){

	my @line = split "\t", $_;
	map {$_=~s/\t|\n//g} @line;
	#print "$line[0]\t$line[2]\n";
	
	# Both sided hash so I have to check only in one direction later	
	$Paralogs{"$line[0]\t$line[1]"} = $line[2]; # Id1 Id2 as key and duplication node as value
}
close (PARALOGS);
#print scalar keys %Paralogs;

# --------------------------------------------------------------------------------------------------------------------

my $currentFile = $combinedOhnoFile; # This is the file under consideration. To start with I assign it to the Ohno file for base organism, and for following iterations, it is assigned to the outfile generated from that iteration. So basically stuff is added to this file name iteratively, and its name also changes.
my $done = ''; # Something to hold a name tag for what's beed already added

print "*** Organisms to be combined: @organisms\n\n";
foreach my $org2add (@organisms){
	
	# $org2add will be the other vertebrate organisms that will be added one by one to the base organism
	if ($org2add ne $organism){ # If it's not the same organism
		
		print "-> Adding $org2add to $currentFile\n";
		# ---------------------------------------------------------------------------------------------------
		if ($orgCount == 1){ # If first organism -> get the final file 
			open CONSENSUS, "$ohnofilepath/$organism/$currentFile" or die $!;
		}
		else { # else use the combined file
			open CONSENSUS, "$organism/$currentFile" or die $!;
		}
		
		# ---------------------------------------------------------------------------------------------------
		# OUTFILE
		my $tag = substr $org2add, 0, 2;
		$done = $done."+$tag";
		
		# If it's the final file with with $orgcount == 26, rename it final file. Otherwise name it as an intermediate
		my $outFile;
		if ($orgCount == $totalOrganisms){
			$outFile = "$organism"."_allVertebrates$wgdtype\.txt";
		}
		else {$outFile = "$organism"."$done$wgdtype\.txt";}
		
		open OUT, ">$organism/$outFile" or die $!;

		print "-> Opening $outFile for writing\n";
		$currentFile = $outFile; # Assign current file to outfile for the next iteration
		
		# ---------------------------------------------------------------------------------------------------
		# The other organism that I wanna add to the above base organism
		my $otherVfile = "$org2add"."Ohno_Self+Outgp$wgdtype\_withDetails.txt";
		print "-> Reading $otherVfile\n";
		
		open VERT2, "$ohnofilepath/$org2add/$otherVfile" or die $!;
		my %vertebrate2;
		while (<VERT2>){
			
			my @line = split "\t", $_;
			map {$_=~s/\t|\n//g} @line;
			#print "$line[0]\t$line[1]\n";
			
			$vertebrate2{"$line[0]\t$line[1]"} = $_;	
		}
		close (VERT2);
		
		# ---------------------------------------------------------------------------------------------------
		# Ortholog file between the two organisms. There can be 2 possible file names, and I check for both		
		my $filename1 = "../3_Orthologs/Vertebrate_Orthologs_filtered/$organism\_$org2add\_orthologs_filtered_Ens84.txt";
		my $filename2 = "../3_Orthologs/Vertebrate_Orthologs_filtered/$org2add\_$organism\_orthologs_filtered_Ens84.txt";
		
		if (-e $filename1){
			print "-> Ortholog $filename1\n" or die $!;
			open ORTH, $filename1 or die $!;
		}
		elsif (-e $filename2){
			print "-> Ortholog $filename2\n" or die $!;
			open ORTH, $filename2 or die $!;
		}
		else {
			die "No ortholog file for $org2add - $organism\n\n";
		}
		
		# Create ortholog hash - both sided
		my %Orthologs;
		while (<ORTH>){
			
			my @line = split "\t", $_;
			map {$_=~s/\t|\n//g} @line;	
			
			$Orthologs{$line[0]}{$line[1]} = '';
			$Orthologs{$line[1]}{$line[0]} = '';
		}
		close (ORTH);
		# --------------------------------------------------------------------------------------------------------------------

		print "-> Adding the orthologous ohnologs\n";
		my %checkedVertIds; # Variable to hold checked ids for the organism that we are adding

		while (<CONSENSUS>){

			my @line = split "\t", $_;
			map {$_=~s/\t|\n//g} @line;
			
			# This is just to print proper header. I am using <CONSENSUS>, otherwise it takes too much memory.
			if (/^ENS.+/){
				print OUT join ("\t", @line),"\t";
			}
			else {
				my $head = $_;
				chomp $head;
				print OUT "$head\t$org2add ohno\t$org2add P-Og\t$org2add P-Self\n";
				next;
			}
		
			
			# Check if the the current ohno pair have one or more orthologous ohno pairs between these vertebrates. And push all ohno pairs in the array @orthPairs.
			if ((exists $Orthologs{$line[0]}) && (exists $Orthologs{$line[1]})){

				my @orthPairs;
				foreach my $gid1 (keys %{$Orthologs{$line[0]}}){

					foreach my $gid2 (keys %{$Orthologs{$line[1]}}){

						if (exists $vertebrate2{"$gid1\t$gid2"}){

							push @orthPairs, $vertebrate2{"$gid1\t$gid2"};
							$checkedVertIds{"$gid1\t$gid2"} = '';
						}
						elsif (exists $vertebrate2{"$gid2\t$gid1"}){

							push @orthPairs,  $vertebrate2{"$gid2\t$gid1"};
							$checkedVertIds{"$gid2\t$gid1"} = '';
						}
						else {
							# No pairs, no problem.
						}
					}
				}
				
				# Get the probability values for printing. If there is none or just 1 orthologous pair, just print it.
				# ****** @orthPairs have the lines from other vertabrate final ohno file.
				# If there are more than 1 -> Get the best one and print that.
				if (scalar @orthPairs <= 1){ # If there is just one line, or no line -> print it

					if (scalar @orthPairs == 0){
						print OUT "\t\t";
					}
					else {				
						my @otLine = split "\t", $orthPairs[0]; # if there is just one line in this array, split it
						map {$_=~s/\n//} @otLine;

						print OUT "$otLine[2],$otLine[3]\t$otLine[7]\t$otLine[15]"; # 2,3 gene symbols. 7: p-outgroup. 13: p-self.
					}
				}
				else { # if there are multiple lines -> Get the best one

					my $ref = SortProbabilities(\@orthPairs);
					my @sortedVline = @$ref;
					map {$_=~s/\n//} @sortedVline;
					print OUT "$sortedVline[2],$sortedVline[3]\t$sortedVline[7]\t$sortedVline[15]";
				}
			}
			else {print OUT "\t\t";}
			print OUT "\n";

		}

		print "-> Adding the orthologous ohnologs in reverse direction\n";
		# --------------------------------------------------------------------------------------------------------------------
		# ************* Now check the reverse direction to see that the ohnolog pairs in the vertebrate file to be added have any that are missed in the comparison above.
		# This reverse comparison is important because there may be some discrepancy in orthologs in different vertebrates. Often for an ohno pair in the other vertebrate, their orthologs are not ohnologs due to duplication time etc. This will get all those also, which I can decide to add or not add later.

		my %ReverseComparison;

			foreach my $gid (keys %vertebrate2){

				if (not exists $checkedVertIds{$gid}){ # I dont need to check here in both the directions here as this is the same vertebarte file and the order would be always 1 and 2

					my ($gid1, $gid2) = split "\t", $gid, 2;

					if ((exists $Orthologs{$gid1}) && (exists $Orthologs{$gid2})){ # Both the genes must have an ortholog
						
						foreach my $id1 (keys %{$Orthologs{$gid1}}){ # These nested loops basically generate all the possible ohnolog pairs for the current gene and push the lines in above array.
							
							foreach my $id2 (keys %{$Orthologs{$gid2}}){
								
								if ((exists $AllPC{$id1}) && (exists $AllPC{$id2}) && ((exists $Paralogs{"$id1\t$id2"}) || (exists $Paralogs{"$id2\t$id1"}))){ # I am only taking the ones that are in ALL PC file and have a paralog in the base organism. 
									# ALL GENES WITHOUT A PARALOG WILL BE IGNORED HERE. I CAN MODIFY THIS LATER IF NEEDED.								
									push @{$ReverseComparison{"$id1\t$id2"}}, $vertebrate2{"$gid1\t$gid2"};
									#print "$id1\t$id2\n";
								}
							}
						}
					}
				}
			}

		foreach (keys %ReverseComparison){
			
			# Print the base file information
			my ($orth1, $orth2) = split "\t", $_, 2;
			print OUT "$orth1\t$orth2\t$AllPC{$orth1}\t$AllPC{$orth2}\t\t";
			if (exists $Paralogs{"$orth1\t$orth2"}){print OUT $Paralogs{"$orth1\t$orth2"};}
			elsif (exists $Paralogs{"$orth2\t$orth1"}){print OUT $Paralogs{"$orth2\t$orth1"};}
			else {die "Check why paralog cannot be found!\n";}
			print OUT "\t\t\t\t\t\t\t\t\t\t\t";

			if ($orgCount > 1){
				for (1..(($orgCount - 1)*3)){
					print OUT "\t";
				}
			}
			
			# Foreach of the ortholog in the other vertebarte file, get the best one and print 
			
			if (scalar @{$ReverseComparison{$_}} <= 1){ # If there is just one line, or no line -> print it

				if (scalar @{$ReverseComparison{$_}} == 0){
					die "No ohno line found\n\n";
				}
				else {				
					my @otLine = split "\t", ${$ReverseComparison{$_}}[0]; # if there is just one line in this array, split it
					map {$_=~s/\n//} @otLine;

					print OUT "$otLine[2],$otLine[3]\t$otLine[7]\t$otLine[15]"; # 2,3 gene symbols. 7: p-outgroup. 13: p-self.
				}
			}
			else { # if there are multiple lines -> Get the best one

				my $ref = SortProbabilities(\@{$ReverseComparison{$_}});
				my @sortedVline = @$ref;
				map {$_=~s/\n//} @sortedVline;
				print OUT "$sortedVline[2],$sortedVline[3]\t$sortedVline[7]\t$sortedVline[15]";
			}
			print OUT "\n";	
		}
		
		$orgCount++;
		print "-> DONE. $orgCount organisms combined\n\n";
	}
}

} # close the organism array ****************************************************************************************************

# Here if there are multiple pairs corresponding to one in the vertebrate under consideration I filter them based on Probabilities as follows:
# Sort based on - P for outgroup (minimum), then P for self comparison (minimum), and then Outgroup support (maximum).
# So basically then we seelct the best ohno pair in that organism
sub SortProbabilities {
	
	my @p = @{$_[0]};
	
	my @temp; # A temporary array having array of arrays
	foreach (@p){
		
		$_ = "$_\n"; # Split wont work properly if there is nothing at the end. There must be at least \n at the end, so i put it there.
		my @line = split "\t", $_;
		
		map {$_=~s/\n//g} @line;
		push @temp, \@line;
	}
	
	no warnings; # no warnings because there may not be ps in some cases - otherwise it prints a lot of junk
	### ********** IMPORTANT: check the column numbers for the probabilities and the outgroup support
	my @sorted = sort {$a->[7] cmp $b->[7] || $a->[15] cmp $b->[15] || $b->[6] <=> $a->[6]} @temp; # sort based on 7:P-og then 13:P-self then 6:Og sup
	
	return (\@{$sorted[0]}); # return the sorted 1st line
}
