# Get a reconciled node from all 7 versions based on majority rule.
# If the consensus points to multiple node, I take prefer versions.
# FishWGD, is a consensus duplication node corresponding to teh base of teleost fish.
#

use strict;
use warnings;

# Process all the organisms -- must be no warnings and errors
foreach (<2_Combined_nodes/*.txt>){

		if ($_ =~/2_Combined_nodes\/(.+)_CombinedDuplicates_Ens80-86\.txt/g){
			#print "$1\n";
			my $org = $1;
			
#			if ($org eq 'hsapiens'){ # tets for one organims	

			print "processing... $org\n";
			Process($org); # Do all the processing
		
#			}
		}
}


# Subroutein 
sub Process {

	my $organism = shift;

	open FH, "2_Combined_nodes\/$organism\_CombinedDuplicates_Ens80-86.txt" or die $!;
	open OUT, ">3_Reconciled_nodes\/$organism\_ReconsiledParalogs_Ens80-86.txt" or die $!;

	foreach (<FH>){

		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
	
		# get the count of nodes for each version
		my %counts;
		foreach (@line[5..11]){
					
			if ($_ ne ''){
				$counts{$_} +=1;
			}
		}
		# If the reconciled node is vertebrate print it
		if ($line[4] eq 'Vertebrates'){
				
			print OUT "$line[0]\t$line[1]\tVertebrata\n";
		}
		# If it is older print it
		if ($line[4] eq 'Older than vertebrates'){
				
			print OUT "$line[0]\t$line[1]\tOlder\n";
		}
		if ($line[4] eq 'FishWGD'){
			
			print OUT "$line[0]\t$line[1]\tFishWGD\n";
		}
		if ($line[4] eq 'Younger than vertebrates'){
				
			print OUT "$line[0]\t$line[1]\tYounger\n";
		}
		# If the node could not be decided
		if ($line[4] eq 'Uncertain'){
			
			my @kys = sort { $counts{$b} <=> $counts{$a} } keys %counts; # Take sorted values in array
			
		
			if ($counts{$kys[0]} > $counts{$kys[1]}){ # If there is a highest value
			
				print OUT "$line[0]\t$line[1]\t$kys[0]\n"; # print it
			}
			elsif ($counts{$kys[0]} == $counts{$kys[1]}){ # if there is no highest value
			
				# print node of the most recent version
				if ($line[11] ne ''){print OUT "$line[0]\t$line[1]\t$line[11]\n";}
				elsif ($line[10] ne ''){print OUT "$line[0]\t$line[1]\t$line[10]\n";}
				elsif ($line[9] ne ''){print OUT "$line[0]\t$line[1]\t$line[9]\n";}
				elsif ($line[8] ne ''){print OUT "$line[0]\t$line[1]\t$line[8]\n";}
				elsif ($line[7] ne ''){print OUT "$line[0]\t$line[1]\t$line[7]\n";}
				elsif ($line[6] ne ''){print OUT "$line[0]\t$line[1]\t$line[6]\n";}
				else {print OUT "$line[0]\t$line[1]\t$line[5]\n";}
			}	
		}	
	}
}

print `date`;
