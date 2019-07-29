# Create the parameter file to run the script
# This is specific to the current directory structure but can be adapted for other contexts also. To adapt it, change the relative drectory paths and try to run manually for an organism.
# This is for the FISH genomes -> since there needs to be two set of ohnologs for fish, 2R and 3R
#

use strict;
use warnings;

# Variables to make directories and the files
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'trubripes');
my @wgd = ('2R', '3R');
my @outgroups_2R = ('dmelanogaster', 'cintestinalis', 'csavignyi', 'celegans', 'bfloridae');
my @outgroups_3R = ('lchalumnae', 'loculatus', 'hsapiens', 'mmusculus', 'oanatinus', 'sscrofa', 'cfamiliaris', 'ggallus');


foreach my $fish (@fish){
	
	# Make directory structuure for fish and tetrapods
	print `mkdir $fish`;
	print `mkdir $fish\/outgroup_2R`;
	print `mkdir $fish\/outgroup_3R`;
	print `mkdir $fish\/self_2R`;
	print `mkdir $fish\/self_3R`;

	# Polyploid genomes: Same for 2R and 3R --------------------------------------------------
	open PP2R, ">$fish/outgroup_2R/PolyploidGenomes.in" or die $!;
	print PP2R "$fish	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt";
	close (PP2R);
	
	open PP3R, ">$fish/outgroup_3R/PolyploidGenomes.in" or die $!;
	print PP3R "$fish	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt";
	close (PP3R);

	open SPP2R, ">$fish/self_2R/PolyploidGenomes.in" or die $!;
	print SPP2R "$fish	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt";
	close (SPP2R);

	open SPP3R, ">$fish/self_3R/PolyploidGenomes.in" or die $!;
	print SPP3R "$fish	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt";
	close (SPP3R);

	
	# Outgroup genomes --------------------------------------------------
	open OG2R, ">$fish/outgroup_2R/OutgroupGenomes.in";
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print OG2R "$og\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$og\_Ens84.txt\n";
		}
	}
	close (OG2R);

	open OG3R, ">$fish/outgroup_3R/OutgroupGenomes.in";
	foreach my $og (@outgroups_3R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print OG3R "$og\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$og\_Ens84.txt\n";
		}
	}
	close (OG3R);	

	open SOG2R, ">$fish/self_2R/OutgroupGenomes.in";
	open SOG3R, ">$fish/self_3R/OutgroupGenomes.in";
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
	 	print SOG2R "$fish\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt\n";
	 	print SOG3R "$fish\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$fish\_Ens84.txt\n";
	}
	close (SOG2R);	
 	close (SOG3R);
	
	# Run parameters -------------------------------------------------	
	open RP2R, ">$fish/outgroup_2R/RunParameters.in" or die $!;
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
			print RP2R "$og\_$_	$fish	$_	$_	2\n";
		}
	}
	close (RP2R);

	open RP3R, ">$fish/outgroup_3R/RunParameters.in" or die $!;
	foreach my $og (@outgroups_3R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
			print RP3R "$og\_$_	$fish	$_	$_	2\n";
		}
	}
	close (RP3R);
	
	open SRP2R, ">$fish/self_2R/RunParameters.in" or die $!;
	open SRP3R, ">$fish/self_3R/RunParameters.in" or die $!;
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
		print SRP2R "$fish\_$_	$fish	$_	$_	2\n";
		print SRP3R "$fish\_$_	$fish	$_	$_	2\n";
	}
	close (SRP2R);
	close (SRP3R);
	
	# Ortholog files --------------------------------------------------
	open ORTH2R, ">$fish/outgroup_2R/OrthologFiles.in" or die $!;
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print ORTH2R "$og\_$_	$fish	../../3_Orthologs/2R_Orthologs_filtered/$fish\_$og\_orthologs_filtered_Ens84.txt\n"; 
		}
	}
	close (ORTH2R);

	open ORTH3R, ">$fish/outgroup_3R/OrthologFiles.in" or die $!;
	foreach my $og (@outgroups_3R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print ORTH3R "$og\_$_	$fish	../../3_Orthologs/3R_Orthologs_filtered/$fish\_$og\_orthologs_filtered_Ens84.txt\n"; 
		}
	}
	close (ORTH3R);

	open SORTH2R, ">$fish/self_2R/OrthologFiles.in" or die $!;
	open SORTH3R, ">$fish/self_3R/OrthologFiles.in" or die $!;
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
	 	print SORTH2R "$fish\_$_	$fish	../../2_Paralogs/4_filtered_paralogs/$fish\_ReconsiledParalogs_Filtered_2R.txt\n"; 
	 	print SORTH3R "$fish\_$_	$fish	../../2_Paralogs/4_filtered_paralogs/$fish\_ReconsiledParalogs_Filtered_3R.txt\n"; 
	}
	close (SORTH2R);
	close (SORTH3R);
	
}


