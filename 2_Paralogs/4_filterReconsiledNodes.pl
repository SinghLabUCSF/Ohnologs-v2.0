# Filter reconciled nodes for 2R WGD organisms -> tetrapods excluding fish
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
	
			print "$1\n";
			my $org = $1;
			
			if ($org eq 'fcatus'){ # test for one organims
			
			if (not exists $Fish{$org}){					
				print "processing ... $org\n";
				Process($org); # Do all the processing
			}
			else {print "Skipped --- $org\n";}		
			}
		}
};

# Subroutein --- just to simplify the reading --------------------------------------------
sub Process {

	my $organism = shift;

	open FH, "4_filter_multi_copy_paralogs\/$organism\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt" or die $!;
	open OUT, ">4_filtered_paralogs\/$organism\_ReconsiledParalogs_Filtered.txt" or die $!;
	
	print OUT "Id1\tId2\n";
		
	foreach (<FH>){
		
		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
		
		# If the reconciled node is vertebrate print it
		if ($line[2] eq 'Vertebrata' || $line[2] eq 'Euteleostomi' || $line[2] eq 'Chordata' || $line[2] eq 'Sarcopterygii' || $line[2] eq 'Neopterygii'){
			print OUT "$line[0]\t$line[1]\n";
		}	
	}
}

print `date`;