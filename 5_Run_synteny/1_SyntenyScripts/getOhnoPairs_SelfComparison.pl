# Generate all possible ohnolog pairs for self comparison
# In self comparison the blocks themselves are ohnologs, so I just need to remove the unwanted columns and lines

use strict;
use warnings;
use diagnostics; 

my $startTime = time;

=cut
my $refGenomeName = 'CionaInt';
my $refGenomeFile = 'AllPCGenes_CionaInt_Ens69.txt';
my $subGenomeName = 'Human';
my $subGenomeFile = 'AllPCGenes_Human_Ens69.txt';
my $orthologyFile = 'BestHits_Human-to-CionaInt.txt';
my $polyploidWindow = 250;
my $referenceWindow = 250;
my $minOrthologs  = 5;
=cut

my $refGenomeName = $ARGV[0];
my $refGenomeFile = $ARGV[1];
my $subGenomeName = $ARGV[2];
my $subGenomeFile = $ARGV[3];
my $orthologyFile = $ARGV[4];
my $referenceWindow = $ARGV[5];
my $polyploidWindow = $ARGV[6];
my $minOrthologs  = $ARGV[7];
my $outFilesDir = $ARGV[8];

my $ohnoCandidateFile = $refGenomeName.'-'.$subGenomeName.'_OhnoCandidates_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';
my $ohnoPairFile = $refGenomeName.'-'.$subGenomeName.'_OhnoPairs_Window-P'.$polyploidWindow.'R'.$referenceWindow.'_Orthologs-'.$minOrthologs.'.txt';

print "Getting pairs of ohnologs...";

#print "$ohnoCandidateFile";

open FH, "$outFilesDir\/$ohnoCandidateFile" or die $!;
my @file = <FH>;
shift @file;

open OUT, ">$outFilesDir\/$ohnoPairFile" or die $!;

foreach (@file){
	
	#print "$_";
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if ($line[5] == 1){
		
		print OUT "$line[6]\t$line[7]\t";
		print OUT join ("\t", @line[10..15]), "\n";
	}
}


print "......................[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}

