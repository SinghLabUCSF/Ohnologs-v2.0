# Read the filteretd Id files and isolate the longest sequences in a new fasta file
use strict;
use warnings;

foreach (<..\/fasta\/*>){ # Foreach organism's fasta file
	
	# get the first and last name
	$_=~/.+\/(.{1}).+_(.+)/g;
	my $name = "$1$2";
	#print "$name\n";
	
	if ($name eq 'hsapiens'){ # test for one gene
	
	# open the gene list file -------------------------------------------------------------------------------------
	open FH, "..\/..\/1_All_Genes\/3_Prepare_final_gene_files_all-scaffolds\/AllGenes_$name\_Ens84.txt" or die $!;
	my @allPC = <FH>;
	shift @allPC;
	
	# make hash for all ensembl ids as keys
	my %filteredIdsList;
	foreach (@allPC){
		
		my @line = split "\t",$_;
		$filteredIdsList{$line[0]} = \@line;
	}
	
	
	# open the sequence file --------------------------------------------------------------------------------------
	foreach (<$_\/pep\/*.fa>){ # foeach here doesn't matter as there is just one fasta file in each folder
	
		print "$_\n";
		local $/ = '>';
		open FH, "$_" or die $!;
		my @seqs = <FH>;
		shift(@seqs);
		
		# make sequence hash with Ensembl ids as keys and value contains an array with all proteins
		# I am appending the protein id to the sequence and later on I will filter the longest sequence
		my %seqHash;
		foreach (@seqs){
			
			my @lines = split "\n", $_;
			if ($lines[-1] eq '>'){pop @lines};
			my $header = shift (@lines);
			
			# get ensembl id and protein ids
			$header =~/^(.+?)\s+.+gene\:(.+?)\s+/g;
			my $protId = $1; my $geneId = $2;
			
			$protId =~s/\.\d+//g;
			$geneId =~s/\.\d+//g;
			#print "$protId\t$geneId\n";
			
			# This is to test if there are some issues with the header -- if this prints something - figure out what's going on
			print "**** \n$header\n****\n" if (not defined $protId);
			print "**** \n$header\n****\n" if (not defined $geneId);
			
			my $seq = join '', @lines;
			
			# append protein id to sequence
			$seq = $protId.'~~'.$seq; # Use ~~ and not _ because some Ids may have _ in them e.g. celegans
			
			
			# push all the proteins for this ensembl id into an array
			push @{$seqHash{$geneId}}, $seq;
		
		}
		
		# Open a new file, filter the ------------ 
		open FH3, ">$name\_LongestProteins.txt" or die $!;
		
		foreach my $ens (keys %seqHash){
			
			# if it exists in all protein coding genes file
			if (exists $filteredIdsList{$ens}){
			
				# sort @{$seqHash{$ens}} (which contains all the proteins for one genes) on the basis of their length
				@{$seqHash{$ens}} = sort {length($a) <=> length($b)} @{$seqHash{$ens}};
			
				# split the protein id and sequence for the longest protein. If two proteins would have same lengths, the last one would be picked.
				my ($prot, $seq) = split '~~', ${$seqHash{$ens}}[-1];
			
				# print Ensembl and protein id as header and sequence
				print FH3 ">$ens|$prot\n$seq\n";
			}
		}


	}
	} # end test for one gene
}

