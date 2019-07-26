# To Filter the genes that have more than 30 paralogs.
# This is because some non-coding classes have thousands of duplicates and runnig with everything will take months.
# NO WARNINGS ALLOWED -> ONLY GENE NAMES SHOULD BE PRINTED ON COMMAND LINE
# Create a directory named 4_filter_multi_copy_paralogs where you run it
#
use strict;
use warnings;

foreach (<3_Reconciled_nodes/*.txt>){
	
	print "$_\n";
	$_=~/.+\/(.+)_ReconsiledParalogs_Ens80-86.txt/g;
	#print "$1\n";
	my $org = $1;
	
if ($org eq 'fcatus'){
	FilterSSD($org);
}	
}

sub FilterSSD {

my $org = shift;

# Open all PC file ------------------------------------------------------------------------------------------------
open ALLGENES, "../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$org\_Ens84.txt" or die $!;
my @allgenes = <ALLGENES>;
shift @allgenes;
my %AllGenes;

foreach (@allgenes){

	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	#if ($line[7] ne 'protein_coding'){
		$AllGenes{$line[0]} = [$line[3], $line[7]];
	#}
	
	#print "$line[0]\t$line[3]\n";
	#last();
};
#print scalar keys %AllGenes;

# Open reconciled paralog files -----------------------------------------------------------------------------------
open PARA, "3_Reconciled_nodes/$org\_ReconsiledParalogs_Ens80-86.txt" or die $!;
my %Paralogs;          # 2D hash paralog pairs => duplication time 
my %ParalogCounts;     # Gene => Total paralogs
my %Duplicates;        # 2D hash paralog pair => '' to determine if it's been checked
my %GeneTypesInParalogs; # Gene type => Total genes with paralogs before filtering
my %GeneTypesInParalogsF; # Gene type => Total genes with paralogs after filtering

foreach (<PARA>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# If the pair has not been checked
	if ((not exists $Duplicates{$line[0]}{$line[1]}) && (not exists $Duplicates{$line[1]}{$line[0]}) &&
	     (exists $AllGenes{$line[0]}) && (exists $AllGenes{$line[1]})){
		
		$Paralogs{$line[0]}{$line[1]} = $line[2];
		$ParalogCounts{$line[0]} += 1;
		$ParalogCounts{$line[1]} += 1;
		#print "$line[0]\t$line[1]\n";
		
		$Duplicates{$line[0]}{$line[1]} = '';
		$Duplicates{$line[1]}{$line[0]} = '';
	}
 }

# Print paralog counts in a file -------------------------------------------------------------------------------
open OUT, ">4_filter_multi_copy_paralogs/$org\_paralog-counts.txt" or die $!;
print OUT "Id	Total paralog pairs	Symbol	Gene type\n";
foreach (keys %ParalogCounts){	
	if (exists $AllGenes{$_}){
		print OUT "$_\t$ParalogCounts{$_}\t";
		print OUT join "\t", @{$AllGenes{$_}}, "\n";
	}
}

# Print all the pairs with number of paralog info -------------------------------------------------------------
open OUTPARA, ">4_filter_multi_copy_paralogs/$org\_ReconsiledParalogs_Ens80-86_RemovedPotentialSSDs.txt" or die $!;
open OUTPARA2, ">4_filter_multi_copy_paralogs/$org\_ReconsiledParalogs_Ens80-86_AllParalogs.txt" or die $!;

print OUTPARA "id1	id2	node	symb1	symb2	type	# paralogs1	# paralogs2\n" or die $!;
print OUTPARA2 "id1	id2	node	symb1	symb2	type	# paralogs1	# paralogs2\n" or die $!;

foreach my $key1 (keys %Paralogs){
	foreach my $key2 (keys %{$Paralogs{$key1}}){
		
		if (${$AllGenes{$key2}}[1] eq 'protein_coding'){
		
			print OUTPARA "$key1\t$key2\t$Paralogs{$key1}{$key2}\t${$AllGenes{$key1}}[0]\t${$AllGenes{$key2}}[0]\t${$AllGenes{$key2}}[1]\t$ParalogCounts{$key1}\t$ParalogCounts{$key2}\n";

			$GeneTypesInParalogsF{${$AllGenes{$key1}}[1]}{$key1} = '';
			$GeneTypesInParalogsF{${$AllGenes{$key2}}[1]}{$key2} = '';

		}
		elsif ($ParalogCounts{$key1} < 30 && $ParalogCounts{$key2} < 30){
			print OUTPARA "$key1\t$key2\t$Paralogs{$key1}{$key2}\t${$AllGenes{$key1}}[0]\t${$AllGenes{$key2}}[0]\t${$AllGenes{$key2}}[1]\t$ParalogCounts{$key1}\t$ParalogCounts{$key2}\n";

			$GeneTypesInParalogsF{${$AllGenes{$key1}}[1]}{$key1} = '';
			$GeneTypesInParalogsF{${$AllGenes{$key2}}[1]}{$key2} = '';

		}
		
		# Print all the pars with additional information
		print OUTPARA2 "$key1\t$key2\t$Paralogs{$key1}{$key2}\t${$AllGenes{$key1}}[0]\t${$AllGenes{$key2}}[0]\t${$AllGenes{$key2}}[1]\t$ParalogCounts{$key1}\t$ParalogCounts{$key2}\n";
		$GeneTypesInParalogs{${$AllGenes{$key1}}[1]}{$key1} = '';
		$GeneTypesInParalogs{${$AllGenes{$key2}}[1]}{$key2} = '';

	}
}


# Print counts of genes in different categories with at least one paralog
open OUTG, ">4_filter_multi_copy_paralogs/$org\_gene-type_counts.txt" or die $!;

print OUTG "Numebr of genes with at least one paralog BEFORE filtering\n";
foreach (keys %GeneTypesInParalogs){	
	print OUTG "$_\t";
	print OUTG scalar keys %{$GeneTypesInParalogs{$_}},"\n";
}

print OUTG "\n\n\n";

print OUTG "Numebr of genes with at least one paralog AFTER filtering the genes with >30 paralogs\n";
foreach (keys %GeneTypesInParalogsF){	
	print OUTG "$_\t";
	print OUTG scalar keys %{$GeneTypesInParalogsF{$_}},"\n";
}

} # end subroutein

print `date`;







