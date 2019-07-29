# Create the output directory structure for fish and tetrapods
#

use strict;
use warnings;

# Variables to make directories and the files
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');

# Make directory structuure for fish and tetrapods

foreach my $tetra (@tetrapods){	
	
	if (!-e $tetra){
		print `mkdir $tetra`;
	}
	else {
		print "Directory $tetra already exists, check and make sure all is well!\n";
	}
}

foreach my $fish (@fish){	
	
	if (!-e $fish){
		print `mkdir $fish`;
	}
	else {
		print "Directory $fish already exists, check and make sure all is well!\n";
	}
}


print `date`;


