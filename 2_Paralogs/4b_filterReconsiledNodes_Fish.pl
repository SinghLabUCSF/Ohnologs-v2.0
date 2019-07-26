# Filter the fish paralogs for WGD
# The nodes selected for 2R are: Vertebrata, Euteleostomi, Chordata, Sarcopterygii (for tetrapods), Neopterygii (for fish)
# Nodes for 3R WGD: FishWGD, Clupeocephala, Acanthomorphata
use strict;
use warnings;

# Define a fish name hash
my %Fish = (
	'drerio' => 'Zebrafish',
	'nfurzeri' => 'Turquoise killifish',
	'olatipes' => 'Medaka',
	'tnigroviridis' => 'Tetraodon' ,
	'trubripes' => 'Fugu',
	'gaculeatus' => 'Stickleback'
);

# Process all the organisms -- must be no warnings and errors
foreach (<4_filter_multi_copy_paralogs/*.txt>){

		if ($_ =~/4_filter_multi_copy_paralogs\/(.+)_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs\.txt/g){
	
			#print "$1\n";
			my $org = $1;
			
#			if ($org eq 'hsapiens'){ # tets for one organims
			
			if (exists $Fish{$org}){					
				print "processing ... $org\n";
				Process($org); # Do all the processing
			}
			else {print "Skipped --- $org\n";}		
#			}
		}
};


sub Process {

	my $organism = shift;	
	
	open FH, "4_filter_multi_copy_paralogs\/$organism\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt" or die $!;
	open OUT2R, ">4_filtered_paralogs\/$organism\_ReconsiledParalogs_Filtered_2R.txt" or die $!;
	open OUT3R, ">4_filtered_paralogs\/$organism\_ReconsiledParalogs_Filtered_3R.txt" or die $!;

	print OUT2R "Id1\tId2\n";
	print OUT3R "Id1\tId2\n";

	foreach (<FH>){
	
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
		
		# If the reconciled node is vertebrate print it
		if ($line[2] eq 'Vertebrata' || $line[2] eq 'Euteleostomi' || $line[2] eq 'Chordata' || $line[2] eq 'Sarcopterygii' || $line[2] eq 'Neopterygii'){
			print OUT2R "$line[0]\t$line[1]\n";
		}
		if ($line[2] eq 'FishWGD' || $line[2] eq 'Clupeocephala' || $line[2] eq 'Acanthomorphata'){
			print OUT3R "$line[0]\t$line[1]\n";
		}
	}
}

print `date`;
