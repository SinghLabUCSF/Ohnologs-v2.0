# Add the amphioxus genes in here from the next step to 2R orthologs

use strict;
use warnings;

foreach (<..\/4_All_PC_Gene_Seqs\/3_run_blast_amphioxus\/BestHits_*.txt>){
	
	#print "$_\n";
	$_=~/.+\/BestHits_(.+)-to-bfloridae.txt/g;
	print "$1\n";
	my $name = $1;
	
#	if ($name eq 'acarolinensis'){ # test for one gene
	
	# open an outfile
	open OUT, ">2R_Orthologs\/$name\_bfloridae_orthologs_BestHits.txt" or die $!;	
	print OUT "ensembl_gene_id	bfloridae_homolog_id	bfloridae_homolog_orthology_confidence\n";

	# open input file	
	open IN, "$_" or die $!;
	
	foreach (<IN>){

		my @line = split "\t", $_;
		map {$_=~s/\n//g} @line;
		
		$line[0] =~s/\|.+//g;
		$line[1] =~s/\|.+//g;
		
		print OUT "$line[0]	$line[1]	$line[2]\n";
		
	}
	close (OUT);
	
#	} # end test for one gene
}

print `date`;