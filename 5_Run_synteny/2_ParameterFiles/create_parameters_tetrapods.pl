# Create the parameter file to run the script
# This is specific to the current directory structure but can be adapted for other contexts also. To adapt it, change the relative drectory paths and try to run manually for an organism.
# This is only for tetrapods
#
use strict;
use warnings;

# Variables to make directories and the files
my @tetrapods = ('acarolinensis', 'fcatus', 'ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'oanatinus', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
#my @tetrapods = ('hsapiens', 'mmusculus'); # test
my @outgroups_2R = ('dmelanogaster', 'cintestinalis', 'csavignyi', 'celegans', 'bfloridae');


foreach my $tetra (@tetrapods){
	
	# Make directory structuure for fish and tetrapods
	print `mkdir $tetra`;
	print `mkdir $tetra\/outgroup`;
	print `mkdir $tetra\/self`;

	# Polyploid genomes --------------------------------------------------
	open PP2R, ">$tetra/outgroup/PolyploidGenomes.in" or die $!;
	print PP2R "$tetra	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$tetra\_Ens84.txt";
	close (PP2R);
	
	open SPP2R, ">$tetra/self/PolyploidGenomes.in" or die $!;
	print SPP2R "$tetra	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$tetra\_Ens84.txt";
	close (SPP2R);

	# Outgroup genomes --------------------------------------------------
	open OG2R, ">$tetra/outgroup/OutgroupGenomes.in";
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print OG2R "$og\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$og\_Ens84.txt\n";
		}
	}
	close (OG2R);

	open SOG2R, ">$tetra/self/OutgroupGenomes.in";
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
	 	print SOG2R "$tetra\_$_	../../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_$tetra\_Ens84.txt\n";
	}
	close (SOG2R);	
	
	# Run parameters -------------------------------------------------	
	open RP2R, ">$tetra/outgroup/RunParameters.in" or die $!;
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
			print RP2R "$og\_$_	$tetra	$_	$_	2\n";
		}
	}
	close (RP2R);

	open SRP2R, ">$tetra/self/RunParameters.in" or die $!;
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
		print SRP2R "$tetra\_$_	$tetra	$_	$_	2\n";
	}
	close (SRP2R);
	
	# Ortholog files --------------------------------------------------
	open ORTH2R, ">$tetra/outgroup/OrthologFiles.in" or die $!;
	foreach my $og (@outgroups_2R){
		
		for (1..5){ # Window sized in the multiples of 100
			$_ = $_*100;
 		 	print ORTH2R "$og\_$_	$tetra	../../3_Orthologs/2R_Orthologs_filtered/$tetra\_$og\_orthologs_filtered_Ens84.txt\n"; 
		}
	}
	close (ORTH2R);

	open SORTH2R, ">$tetra/self/OrthologFiles.in" or die $!;
	for (1..5){ # Window sized in the multiples of 100
		$_ = $_*100;
	 	print SORTH2R "$tetra\_$_	$tetra	../../2_Paralogs/4_filtered_paralogs/$tetra\_ReconsiledParalogs_Filtered.txt\n"; 
	}
	close (SORTH2R);
	
}


