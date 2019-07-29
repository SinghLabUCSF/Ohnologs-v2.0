# Reads ohno pair files and makes ohno families using depth first search algorithm.
#

use strict;
use warnings;

local $| = 1;

# --------------------------------------------------------------------------------------------------------------------
my @allorgs = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus', 'gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');
my @tetrapods = ('acarolinensis', 'fcatus','ggallus', 'ptroglodytes', 'btaurus', 'cfamiliaris', 'ggorilla', 'ecaballus', 'hsapiens', 'mmulatta', 'cjacchus', 'mmusculus', 'panubis', 'mdomestica', 'pabelii', 'sscrofa', 'ocuniculus', 'rnorvegicus', 'oaries', 'mgallopavo', 'csabaeus', 'tguttata', 'loculatus');
my @fish = ('gaculeatus', 'olatipes', 'tnigroviridis', 'drerio', 'nfurzeri');

# get command line arguments
my $subGenomeName = $ARGV[0];
my $criteria = $ARGV[1];
my $wgd = $ARGV[2]; # _3R or _2R for fish

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
if ($subGenomeName ~~ @tetrapods){$wgd = '';} # Remove the WGD type for tetrapods
print "> Processing $subGenomeName criteria $criteria $wgd\n";

# --------------------------------------------------------------------------------------------------------------------
my $subGenomeFile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$subGenomeName.'_Ens84.txt';
#my $subGenomeFile = '../1_All_Genes/3_Prepare_final_gene_files_all-scaffolds/AllGenes_'.$subGenomeName.'_NCBI.txt'; ## For N furzeri

my $polyploidChrGenes = "../5_Run_synteny/3_SyntenyOutputFiles/$subGenomeName\/self$wgd\/".$subGenomeName.'_Chromosomes-Genes.txt';
my $outFile = $subGenomeName.'/Families_Criteria-['.$criteria.']_'."$subGenomeName"."$wgd\.txt";
my $pairFile = "../7_FilterOhnologs/$subGenomeName/$subGenomeName\_Criteria-\[$criteria\]-Pairs$wgd\.txt";
my $polyploidWindow = 100;
my $shrink = 'n';
# --------------------------------------------------------------------------------------------------------------------

print "  PAIRS: $pairFile\n";
print "  GENOME: $subGenomeFile\n";
print "  CHROMOSOMES: $polyploidChrGenes\n";
print "  OUT: $outFile\n";
print "  Generating families...";
# Read file with polyploid chromosomes and genes on them
open PP, $polyploidChrGenes or die $!;
my @ppchrgn = <PP>;
close(PP);

my %ppChrGenes;
foreach (@ppchrgn){

	my @line = split "\t", $_;
	map {$_=~s/\n|\r//g} @line;
	
	$ppChrGenes{$line[0]} = $line[1];
}


# Open Ohnolog Pairs file and take them into an array ------------------------------
open OHNOPAIRS, $pairFile or die "$! \n$pairFile\n";
my @OHNOpairs = <OHNOPAIRS>;
shift @OHNOpairs;
close (OHNOPAIRS);


my %Explored; # Explored nodes
my %Relations; # Genes and all it's ohnologs
my %Locations;
my @family; # Ohnolog family identified by depth first search

# Read ohno pair file
foreach (@OHNOpairs){
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\r//g} @line;
	
	push @{$Relations{$line[0]}}, $line[1];  # Hash having each gene as id and all its ohnolog partners as value array
	push @{$Relations{$line[1]}}, $line[0];
}

# Read WGD genome file and get locations of genes on chromosomes
open WGDGENOME, $subGenomeFile or die $!;
my @ppGenome = <WGDGENOME>; 
close (WGDGENOME);
shift @ppGenome;

my $locRef = getLocationsOnChromosomes(\@ppGenome);
%Locations = %$locRef;

#foreach (keys %Locations){
#	print "$_\t$Locations{$_}\n";
#}


# Open combined family file
open FAMILY, ">$outFile" or die $!, " $outFile";

#print scalar keys %Relations;
#print "@{$Relations{'ENSGACG00000010605'}}\n";

#foreach my $kk (keys %OutgroupSupport){
#	
#	my %hash = %{$OutgroupSupport{$kk}};
#	
#	foreach (keys %hash){
#		print "$kk\t$_\t$hash{$_}\n";
#	}
#}


foreach (keys %Relations){ # For each ohnolog gene
	
	my @line = split "\t", $_;
	map {$_=~s/\n|\r//g} @line;
	
	if (not exists $Explored{$_}){                       # If it has not been explored already
		
		my @fam = getConnections(\%Relations, $_);   # get the family from this subroutein. Pass the current node and all relations
		
		if (scalar @fam > 1){
			#print "@fam\n";
			#if ($fam[0] eq 'ENSG00000100450') {clubSSDTogether(\@family, \%Locations);}           # In this family club together the SSD genes from this subroutein
			clubSSDTogether(\@family, \%Locations);
		}
		
		@family = (); # Empty the family array as it's a global one
	}

}

sub getConnections {
	
	my $graphRef = shift;
	my %Graph = %$graphRef; # Has all genes and ohnolog relations
	my $node = shift;       # The curent node
	#my @family; # Ohnolog family identified by depth first search
	
	#print "**$node** : ";
	
	if (not exists $Explored{$node}){push @family, $node;} # If the node has not been explored, start family using this node
	$Explored{$node} = 1;                                  # and mark it explored
	#print "@{$Graph{$_}}\n\n";
	foreach (@{$Graph{$_}}){    # For each of its ohnolog partners in graph. $_ here is $node. Check it properly later !!!
		
		my $nbor = $_;      # Take one partner or ohnolog neighbour at a time

		if ((not exists $Explored{$nbor})){ # if this one is not explored
			
			$Explored{$nbor} = 1; # mark it as explored
			
			push @family, $nbor;  # put it in family
	#		print "*\t$node\t$nbor\t$ogsup\n";
	#		print "Family: @family\n\n";
			getConnections(\%Relations, $nbor); # and recursively call the function to get all it's relations. It will keep on doing it as long as there are no nodes unexplored.
		}
	}
	#print "-------------\n< @family >\n------------\n";
	return @family; # Return this family
}


sub clubSSDTogether {
	
	my $ref = shift;
	my @Fam = @$ref;
	my $loc = shift;
	my %Position = %$loc; # This position hash has Ids as key and chr~pos as value.
	my %Ids = reverse %Position; # I use reverse to reverse its key and values
	my %OHNO;
	
	my %sorted; # hash to hold gene ids sorted based on chromosomes and gene count
	
	foreach (@Fam){
		
		my ($c, $p) = split '~', $Position{$_}; # get chromosome and position
		
		push @{$sorted{$c}}, $p; # Sort ids based on chromosomes
		#print "$_\t$c - $p\n";
	}
	
	foreach my $chr (keys %sorted){ # Foreach chromosome
		
		#print "$chr\t=>\t";
		my @pos = sort {$a <=> $b} @{$sorted{$chr}}; # Sort genes based on their position
		#print join '-', @pos;
		#print "\n";
		
		if (scalar @pos == 1){ # if there is just one gene on chromosome...
			#print FAMILY "$chr - $pos[0]\n";
			$OHNO{$Ids{"$chr~$pos[0]"}} = ''; # Push it in family
		}
		else { # if there are more than one genes
			my %ohno;
			my $count = 0; # This count is just like the 0,1,2 positions on chromosomes like earlier
			
			for (my $i = 0; $i < (scalar (@pos)-1); $i++){ # For all positions
				
				# Get window boundaries for both the consecutive genes
				my $leftWd_i; my $rightWd_i;
				my $leftWd_ip1;my $rightWd_ip1;
				if ($shrink eq 'n'){
					($leftWd_i, $rightWd_i) = getAsymmatricWindowPosition($pos[$i], $ppChrGenes{$chr}, ($polyploidWindow/2));
					($leftWd_ip1, $rightWd_ip1) = getAsymmatricWindowPosition($pos[$i+1], $ppChrGenes{$chr}, ($polyploidWindow/2));
				}
				else {
					($leftWd_i, $rightWd_i) = getSymmatricWindowPosition($pos[$i], $ppChrGenes{$chr}, ($polyploidWindow/2));
					($leftWd_ip1, $rightWd_ip1) = getSymmatricWindowPosition($pos[$i+1], $ppChrGenes{$chr}, ($polyploidWindow/2));
				}
				
				$ohno{$count}{$pos[$i]} = ''; # Push the current count and position in the nested hash
				#print "\t$rightWd_i  $leftWd_ip1\n";
				if ($rightWd_i < $leftWd_ip1){ # If the genes lie outside the window boundary
					$count++;              # Increase the count
					$ohno{$count}{$pos[$i+1]} = ''; # Push the count and position of next gene in a nested hash
					#print "$count\t$pos[$i]\t$pos[$i+1]\n";
				}
				else {
					$ohno{$count}{$pos[$i+1]} = ''; # else mark next gene also in the same count category
					#print "$count\t$pos[$i]\t$pos[$i+1]\n";
				}
			}
			
			foreach (keys %ohno){ #  Foreach count
				
				my $str = ''; # Initialize the variable str to hold all genes on current chromosome
				
				foreach (keys %{$ohno{$_}}){
					
					#print FAMILY "$chr - $_|";
					$str = $str.$Ids{"$chr~$_"}.'|'; # Foreach position mark all genes as SSD, join them by |
				}
				#print FAMILY "\n";
				$str =~s/\|$//g;
				$OHNO{$str} = ''; # Push structure in OHNO hash
			}
		}	
	} # end chr loop 
	
	#print FAMILY "\n"; # Separater for new family
	
	print FAMILY scalar keys %OHNO, "\t"; # join different families by \t and print
	print FAMILY join ("\t", keys %OHNO),"\n";
	#exit();
}
print ".......... DONE\n\n";


sub getAsymmatricWindowPosition {

	my $pos = shift; # position
	my $N = shift;   # Total number of genes
	my $wd = shift;  # This is already a half window
	
	#print "$pos\t$N\t$wd\n";
	
	my $leftWd; my $rightWd;
	my $lOffset; my $rOffset;
	
	# generate regular symmatric windows
	$leftWd = $pos - $wd;
	$rightWd = $pos + $wd;
	
	# If left window has a problem reset it and the right wd
	if ($leftWd <= 0){
		
		$lOffset = abs($leftWd) + 1;
		$rightWd = $rightWd + $lOffset;
		$leftWd = 1;
	}
	
	# Fix the right window is its greater than chromosome
	if ($rightWd > $N){
			
		$rOffset = $rightWd - $N;
		$leftWd = $leftWd - $rOffset;
		$rightWd = $N;
	}
	
	# If still they exceed the chromosome limit - due to chromosome smaller than windows - reset them to boundaries
	if ($leftWd <= 0){$leftWd = 1};
	if ($rightWd > $N){$rightWd = $N};
	
	return ($leftWd, $rightWd);
}

sub getSymmatricWindowPosition {

	my $pos = shift;
	my $N = shift;
	my $wd = shift;
	
	my $leftWd; my $rightWd;
	
	# Get the window boundaries simply by adding the half window to both sides of the gene
	# So the total number of the genes will be complete window +1 
	# eg if position is at 755 window of size 200  will be 755+100 and 755-100 = 655 to 855 (total 201 genes including 655th and 855th gene)
	# For the windows to not overlap the anchors sitting in the middle should be more than wd size apart for symmatric windows and middle segments
	
	$leftWd = $pos - $wd;
	$rightWd = $pos + $wd;
	
	# If the window boundaries exceed the segment length reset them to segment boundaries
	if ($leftWd <= 0){$leftWd = 1};
	if ($rightWd > $N){$rightWd = $N};
	
	return ($leftWd, $rightWd);
	
}


sub getLocationsOnChromosomes {

	my @wgdGenome = @{$_[0]};
	my %localPos;
	my $wgdGenomeCount = 1;
	
	for (my $i = 0; $i < scalar (@wgdGenome); $i++){

		# Read two lines of the file at a time
		my @currentLine = split "\t", $wgdGenome[$i];		
		my @previousLine = split "\t", $wgdGenome[$i-1];
		map {$_=~s/\n|^\s+|\s+$|\r//g} @currentLine;
		map {$_=~s/\n|^\s+|\s+$|\r//g} @previousLine;

		#print "$currentLine[2]\t$previousLine[2]\t$wgdGenomeCount\n";

		# If the chromosome number is different then reset the count to 1
		if ($currentLine[2] ne $previousLine[2]){
			$wgdGenomeCount = 1;
		}

		$localPos{$currentLine[0]} = $currentLine[2].'~'.$wgdGenomeCount;
		$wgdGenomeCount++;
	}
	
	return (\%localPos);
}
