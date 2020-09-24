# 20160928: In this version I include all the genes that are also on unplaced scaffolds in Ensembl.
#
# to process the Biomart downloaded file and make a final file for 
# OHNOLOGS server and runs.
use strict;
use warnings;

local $| = 1;

# This is the hash for pattern of chromosomes to be removed. 
# I just exclude MT and some scaffolds that are haplotypes
my %Chrs = (
	     'hsapiens'		=> '^MT$|^CHR_.+$',
	     'mmusculus'	=> '^MT$|^CHR_.+$',
	     'rnorvegicus'	=> '^MT$',
	     'sscrofa'		=> '^MT$',
	     'cfamiliaris'	=> '^MT$',
         'ggallus'  	=> '^MT$',
         'olatipes'		=> '^MT$',
         'drerio' 		=> '^MT$',
         'acarolinensis'	=> '^MT$',
         'fcatus'		=> '^MT$',
         'ptroglodytes'	=> '^MT$',
         'btaurus'		=> '^MT$',
         'ggorilla'		=> '^MT$',
         'ecaballus'	=> '^MT$',
         'mmulatta'		=> '^MT$',
         'cjacchus'		=> '^MT$',
         'panubis'		=> '^MT$',
         'mdomestica'	=> '^MT$',
         'pabelii'		=> '^MT$|^.+Hap.+',
         'oanatinus'	=> '^MT$',
         'ocuniculus'	=> '^MT$',
         'oaries'		=> '^MT$',
         'mgallopavo'	=> '^MT$',
         'csabaeus'		=> '^MT$',
         'tguttata'		=> '^MT$',
         'loculatus'	=> '^MT$',
         'gaculeatus'	=> '^MT$',
         'tnigroviridis'	=> '^MT$',
         'trubripes'	=> '^MT$',
         'cintestinalis'	=> '^MT$',
         'csavignyi'	=> '^MT$',
         'lchalumnae' => '^MT$',
         'dmelanogaster'	=> 'dmel_mitochondrion_genome',
         'celegans'		=> '^MT$'
);


foreach (<1_BioMart_gene_attributes\/*.txt>){

	my $infile = $_;
	#print "$_\n";
	
	$_=~/\/(.+)_gene_ensembl_biomaRt_v84\.txt/g;
	my $organism = $1;
		
#	if ($organism eq 'cjacchus'){ # test for one organism
	
	print "$organism\n";
	
	# output directory
	my $directory = "3_Prepare_final_gene_files_all-scaffolds";

	# created if not exists
	unless(-e $directory or mkdir $directory) {
        die "Unable to create $directory\n";
    }

	my $outfile = $directory."\/AllGenes_".$organism."_Ens84.txt";
	
	# Open the organism file downloaded from bio_mart
	open FH, $infile or die $!;
	my @file = <FH>;
	close (FH);
	my $header = shift @file;
		
	# open outfile
	open OUT, ">$outfile" or die $!;

	# process file to filter genes and make a hash of filtered genes
	my %AllPC;

	foreach (@file){
	
		my @line = split "\t", $_;
		map {$_=~s/\n|\r|^\s+|\s+$//g} @line;
		
		# Filter the line if the gene is on one of the chromosomes in the karyotypes AND is one of the pre-defined gene_type
		# I manually selected these gene types as they have ortholog and paralog information in Ensembl
		if (($line[2] !~/$Chrs{$organism}/gi) && ($line[7] =~/miRNA|mirna|misc_RNA|protein_coding|rRNA|snRNA|snoRNA/gi)){
			
			# replace symbol with ids if its blank
			if ($line[3] eq ''){$line[3] = $line[0];}
			
			# remove extra stuff from the description column
			$line[8] =~s/\s{0,2}\[.+\].*//g;
			
			# push proper ids in AllPC hash
			$AllPC{$line[0]} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
			
			#print OUT join ("\t", @line),"\n";
		}
	}
	
	
	# Open the GO file for this organism
	open GO, "2_BioMart_GO_attributes\/$organism\_gene_ensembl_GO_v84.txt" or die $!;
	my @gofile = <GO>;
	close (GO);
	my $goheader = shift @file;
	
	# Hashes to hold GO
	my %GOId;
	my %GOTerm;
	
	# get GO terms and Ids in a hash to be appended to the all genes file
	foreach (@gofile){
	
		my @line = split "\t", $_;
		map {$_=~s/\n|\r+|^\s+|\s+$//g} @line;
		
		#print "$line[0]\t$line[2]\t$line[3]\n";
		
		# Create two hashes for GOID and Terms
		if ($line[2] ne ''){push @{$GOId{$line[0]}}, $line[2]};
		if ($line[3] ne ''){push @{$GOTerm{$line[0]}}, $line[3]};	
	
	}	

	# Proecss to add GO terms and Ids
	foreach (keys %AllPC){
	
		my @line = split "\t", $AllPC{$_};
		map {$_=~s/\n|^\s+|\s+$//g} @line;
			
		print OUT join ("\t", @line[0..8]),"\t";
		
		if (exists $GOId{$line[0]}){
			print OUT join (',', @{$GOId{$line[0]}}),"\t";
		}
		else {print OUT "\t";}
	
		if (exists $GOTerm{$line[0]}){
			print OUT join (',', @{$GOTerm{$line[0]}});
		}
		else {print OUT "";}
	
		print OUT "\n";		
	
	}
	close (OUT);

	# open the file again to sort based on chromosome and start position
	open FH2, $outfile or die $!;
	my @allpc; # 2D hash to sort
	foreach (<FH2>){

		my @line = split "\t", $_;
		map {$_=~s/\n|^\s+|\s+$//g} @line;
		push @allpc, \@line;
	}
	close (FH2);

	# delete it to sort and print
	unlink "$outfile";

	# Sort gives warnings because some chromosomes are not numbers, so I ignore them for this line
	no warnings;
	#my @sortedallpc = sort { $a->[2] <=> $b->[2] || $a->[4] <=> $b->[4]} @allpc;
	my @sortedallpc = sort { $a->[2] cmp $b->[2] || $a->[4] <=> $b->[4]} @allpc;
	use warnings;
	
	# open another outfile
	open OUT2, ">$outfile" or die $!;
	
	# Print header in the file
	$header =~s/\n//g;
	print OUT2 $header,"\tgo_id	name_1006\n";
	foreach (@sortedallpc){

		print OUT2 join "\t", @{$_};
		print OUT2 "\n"
	}
	#} # End bracket for one organism
} 

print "done";
