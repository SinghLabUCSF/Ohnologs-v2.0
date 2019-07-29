# Get the families that are cliques i.e. each member is ohnolog with every other member. 
# I am not using them for now.
#
use strict;
use warnings;
use Algorithm::Combinatorics qw/combinations permutations/;
local $| = 1;
# --------------------------------------------------------------------------------------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');

# get command line arguments
my $organism = $ARGV[0];
my $criteria = $ARGV[1];
my $wgd = $ARGV[2]; # _3R or _2R for fish

#my $organism = "ptroglodytes";
#my $criteria = "E";
#my $wgd = "_2R"; # _3R or _2R for fish


if ((scalar @ARGV) != 3){
	
	print "Script takes 3 arguments:
	1st: organism name e.g. hsapiens
	2nd criteria (A, C or E)
	3rd WGD type _2R or _3R
	
	Example usage: 1_DepthFirstSearchOhnolgFamilies_cl.pl drerio A _3R\n\n";
	print "Check parameters. Criteria can only be A, C or E. WGD can only be _2R and _3R and organism name must be one of the following\n";
	print "@allorgs\n\n";
	exit;
}

# Remove WGD if it's not fish
if ($organism ~~ @tetrapods){$wgd = '';} # Remove the WGD type for tetrapods

print "> Processing $organism criteria $criteria $wgd\n";

# -------------- INPUT & OUTPUT FILES -------------------------------------------------------
my $pairFile = "../7_FilterOhnologs/$organism/$organism\_Criteria-\[$criteria\]-Pairs$wgd\.txt";
my $familyFile = "$organism/Families\_Criteria-\[$criteria\]_$organism\_ProcessedFinal$wgd\.txt";
my $outfile = "$organism/Families\_Criteria-\[$criteria\]_$organism\_ProcessedFinal_Cliques$wgd\.txt";
# ----------------------------------------------------------------------------------

print "  Pair file   : $pairFile\n";
print "  Input       : $familyFile\n";
print "  Outfile     : $outfile\n";
print "  -> Warning: Double check that the correct nodes are in place (including or exclusing Sarcopterygii and Neopterygii).\n\n";

# Make an outfile 
open OUT, ">$outfile" or die $!;

# Read the ohno pair file and make a hash in both side to check connections
open PAIRS, $pairFile or die $!;
my %Pairs;
foreach (<PAIRS>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Pairs{$line[0]}{$line[1]} = '';
	$Pairs{$line[1]}{$line[0]} = '';
}

# Read the family file and go through them one by one
open FAM, $familyFile or die $!, $familyFile;
foreach (<FAM>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	# get family size
	my $famsize = shift @line;
	#print join ("\t", $famsize, @line), "\n";
	
	if ($famsize == 2){ # If he family size is 2, just print and go to next family
		
		print OUT join ("\t", $famsize, @line), "\n";
		next;
	}
	
	my $i = 1;
	my %Family; # This will have family id number starting at 1 as key, and ensembl id array as value. This will be used for clique identification.
	my %EnsemblId;
	
	#print "populating pair hash ...\n";
	foreach (@line){ # get each SSD members too
		
		my @members = split '\|', $_;
		foreach (@members){
			
			push @{$Family{$i}}, $_; # Populate a family hash with numbers starting from 1 for each member. This is because I found an easy code to identify clique and this is what it needs. 
			$EnsemblId{$i} = $_;     # To get the Ensembl id back from number. 

		}
		$i++;
	}
	
	# test print
	#foreach (sort {$a <=> $b} keys %Family){
	#	print "$_\t@{$Family{$_}}\n";
	#}	
	
	# *****************
	my @V = (1..$famsize); # This is the vertex array for clique identification
	my %E;                 # This is the edge array for clique identification
	my @E_ens;			   # Edge array with ensembl id
	
	#print "Getting all edges...\n";
	# Populate the %E with all ohnolog pairs 
	for (my $i = 1; $i < $famsize; $i++) {
		for (my $j = $i+1; $j <= $famsize; $j++) {

			#print "-> $i	$j\n";
			
			foreach my $id1 (@{$Family{$i}}){
				foreach my $id2 (@{$Family{$j}}){
					
					if (exists $Pairs{$id1}{$id2}){
						# print "$id1\t$id2*\n";
						$E{$i}{$j} = '';
						push @E_ens, [$id1, $id2];
						last;
					}
				}
			}
		}
	}
	
	# test print
 	#print "@V\n"; 	
	#foreach my $k1 (sort {$a <=> $b} keys %E){
	#	foreach my $k2 (sort  {$a <=> $b} keys %{$E{$k1}}){
	#		print "$k1 $k2 | ";
	#	}
	#}
	#print "\n";
 	#print "Total edges: ", scalar @E_ens, "\n";
 	#print "Total veretices: ", scalar @V, "\n";


	# ------------ okay, edges are ready now. Now I need to generate all combinations of cliques
	# e.g. For a family with 4 members: ABCD, I need to test if there are cliques of size 4 and 3
	# like ABCD, ABC, ABD, BCD. SO I use the package combinatorics to genegrate all combinations and then test each of them one by one.
	
	my @Cliques; # This will have all the cliques
	for my $k (reverse 3 .. $famsize){ 	# Minimum clique size will be 3, so in descending order
		
		# This is to generate all comnbinations using the package Combinatorics
		# http://blogs.perl.org/users/dana_jacobsen/2015/02/short-survey-of-modules-for-combinations-and-permutations.html
		my $citer = combinations(\@V, $k);
 		#print "Getting cliques of size $k\n";
		
		# check if this has enough edges to be a clique. The idea is that if the edges are short of total required for a clique for size $k, there is no point in checking it. This will save me time, for large families.
		my $required_edges = ($k * ($k - 1)) / 2; # this is how many edges are needed fir a clique of size $k
		if (scalar @E_ens < $required_edges){
			#print "Need $required_edges for a clique of $k, found only ", scalar @E_ens, "\n";
			next;
		}
	
		while (my $c = $citer->next) {
	
			#print "@$c\t";	
				
			# Check if the family already exists with a bigger clique e.g. if ABCD exists, skip ABC or ABD etc.
			my $flag = 0;
			foreach my $f (@Cliques){
				
				my @intersection =	grep { defined } 
				                    @{ { map { lc ,=> $_ } @{$f} } }
                                    	 { map { lc } @$c };
                 
                #print "* @{$f} ";
                 
             	if ((scalar @intersection == scalar @$c)){
             		#print "\t ~~~> A bigger family exists.\n";
					$flag = 1;
					last;
             	}
             	#print "\n";
			}
			if ($flag == 1){next;}					
			
			# ************ CHECK FOR CLIQUE **************
			my $clicRef = getCliques(\@$c, \%E, $k); # Subroutein to check if it's a clique or not. It returns members of the clique.
			my @cFam = @{$clicRef};
			
			if (!@cFam){ # If this is empty, it means it's not a clique
				#print " -> Not clique\n";
			}
			else { # Else this is a clique
				#print " ==> Clique!!! *@cFam*\n";
				push @Cliques, \@cFam; 
			}
		}

	}
	
	# The cliques are done. I will print them now. Here's what I need to do
	# 1. Print the cliques. ***** REMEMBER THAT THE CLIQUES MAY BE OVERLAPPING *****
	# 2. Then print the pairs that are not part of any clique. Some of the genes may also be in cliques.

	my %coveredFamIndex; # This will have index of the pairs that are in the cliques
	
	if (!@Cliques){ # If the array is empty, just print the pairs
		for (my $i = 0; $i <= $#E_ens; $i++){
			print OUT "2\t${$E_ens[$i]}[0]\t${$E_ens[$i]}[1]\n";
		}
	}
	else { # If there are cliques
		foreach my $f (@Cliques){

			# print the clique
 			#print join ("\t", @$f), "\n";
 			print OUT join ("\t", (scalar @$f),
				(map {join '|', @{$Family{$_}}} @$f)), "\n";

			# Then identify the pairs that are covered by these cliques, below until ~~~
			my %localfam;
			foreach my $clc (@$f){
				$localfam{$_} = 1 foreach @{$Family{$clc}};
			}
			#print join ("__", keys %localfam), "\n";
			
			for (my $i = 0; $i <= $#E_ens; $i++){
				if ((exists $localfam{${$E_ens[$i]}[0]}) && (exists $localfam{${$E_ens[$i]}[1]})){
					#print "*${$E_ens[$i]}[0]\t${$E_ens[$i]}[1]\n";
					$coveredFamIndex{$i} = 1;
				}
			}
			# ~~~
		}

		# Now print the pairs that are not in clique families
		for (my $i = 0; $i <= $#E_ens; $i++){
			if (not exists $coveredFamIndex{$i}){
				print OUT "2\t${$E_ens[$i]}[0]\t${$E_ens[$i]}[1]\n";
			}
		}	
	}
}


# Subroutein to get clique. It counts all the connections for a set of edges and checks if all of them are connected.
sub getCliques {
	
	my $Vref = shift;
	my $Eref = shift;
	my $k = shift;      # The clique size
	my @V = @{$Vref};   # The vertex array
	my %E = %{$Eref};   # The edge hash
	my @clq;
	
	# Compute how many total vertex in @V have a connection in %E
	my $edgecount = 0;
	for (my $i = 0; $i <= $#V; $i++) {
		for (my $j = $i+1; $j <= $#V; $j++) {
			if ((exists $E{$V[$i]}{$V[$j]}) || (exists $E{$V[$j]}{$V[$i]})){
				$edgecount++; 
			}
		}
	}
	
	# If this is a clique, the connections must be same as total required for a clique of size $k
	my $required_edges = ($k * ($k - 1)) / 2; # this is how many edges are needed fir a clique of size $k
	if ($edgecount == $required_edges){
		#print "*$edgecount* = $required_edges !!! Clique\n";
		@clq = @V;
	}
	elsif ($edgecount < $required_edges) {
		#print "*$edgecount* != $required_edges. NOT Clique\n";	
	}
	else {
		print "*$edgecount*, $required_edges *** Error *** \n"; # Print error if the number of edges are more than the clique, something's must be wrong.
	}	

	return (\@clq);

}
