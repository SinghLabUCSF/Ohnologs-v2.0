# Create the output directory structure for fish
#

use strict;
use warnings;

# Variables to make directories and the files
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes');
my @wgd = ('2R', '3R');

foreach my $fish (@fish){
	
	# Make directory structuure for fish and tetrapods
	print `mkdir $fish`;
	print `mkdir $fish\/outgroup_2R`;
	print `mkdir $fish\/outgroup_3R`;
	print `mkdir $fish\/self_2R`;
	print `mkdir $fish\/self_3R`;
	
}

print `date`;

