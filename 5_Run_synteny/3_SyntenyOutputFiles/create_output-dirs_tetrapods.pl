# Create the output directory structure for tetrapods
#
use strict;
use warnings;

# Variables to make directories and the files
my @tetrapods = ('acarolinensis', 'fcatus', 'ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
#my @tetrapods = ('hsapiens', 'mmusculus'); # test

foreach my $tetra (@tetrapods){
	
	# Make directory structuure for fish and tetrapods
	print `mkdir $tetra`;
	print `mkdir $tetra\/outgroup`;
	print `mkdir $tetra\/self`;
	
}

print `date`;
