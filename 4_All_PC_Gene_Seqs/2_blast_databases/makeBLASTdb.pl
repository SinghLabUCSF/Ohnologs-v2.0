# Make blast database for each organism
use strict;
use warnings;

foreach (<..\/1_get_longest_seqs\/*.txt>){ # Foreach organism's fasta file
	
	# get the first and last name
	$_=~/.+\/(.+)_LongestProteins.txt/g;
	my $name = "$1";
	print "$1\n";
	
#	if ($1 eq 'celegans'){
	print `makeblastdb -in $_ -dbtype prot -title $1 -out $1 >> makeblastdb.log.txt 2>&1`; # 2>&1 is to redirect both stderr and stdout to the same file
#	}
}
