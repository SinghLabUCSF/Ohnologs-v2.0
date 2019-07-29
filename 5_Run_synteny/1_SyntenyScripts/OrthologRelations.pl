# Make ortholog relations based on the gene position information in
# the polyploid and the outgrouop genome file
#
# Input : 1. Outgroup genome file with genes sorted based on chromosomes and then the order they appear on chromosomes
#         2. Plyploid genome file with genes sorted based on chromosomes and then the order they appear on chromosomes
#         3. Ortholog file having ids and positions of the ortholgs in 2 columns
#            
# Output: 1. Ortholog relation file
#
# IMPORTANT : 1. Order of genes and column do not matter in ortholog file
#             2. Genes should be sorted based on the chromosomes and positions from 1 to n for genome files
#             3. Genes from one chromosome should be together in genome files
#             4. Orthologs which do not exist in genome files will be ignored

use strict;
use warnings;
use diagnostics; 

my $startTime = time;

my $refGenomeName = $ARGV[0];
my $refGenomeFile = $ARGV[1];
my $subGenomeName = $ARGV[2];
my $subGenomeFile = $ARGV[3];
my $orthologyFile = $ARGV[4];
my $referenceWindow = $ARGV[5];
my $polyploidWindow = $ARGV[6];
my $minOrthologs  = $ARGV[7];
my $outFilesDir = $ARGV[8];
print "**$outFilesDir**\n";
my $orthologRelationFile = $refGenomeName.'-'.$subGenomeName.'_Ortholog_Relations.txt';
my $polyploidChrGenes = $subGenomeName.'_Chromosomes-Genes.txt';
my $outgroupChrGenes = $refGenomeName.'_Chromosomes-Genes.txt';
#----------------------------------------------------------------------------------------------------------

print "Input Parameters: 

Outgroup Genome Name             : $refGenomeName
Outgroup Genome File             : $refGenomeFile
Polyploid Genome Name            : $subGenomeName
Polyploid Genome File            : $subGenomeFile
Orthology Relationship File      : $orthologyFile
Window Size in Polyploid Genome  : $polyploidWindow       
Window Size in Reference Genome  : $referenceWindow       
Minimum Ortholog Pairs In Window : $minOrthologs
\n";

print "Creating ortholog relation file...";

#------------------- OPEN FILES ------------------------------#
# Outgroup genome: you want to keep at X-axis in synteny plots
open OUTGROUP, $refGenomeFile or die $!;
my @outgroup = <OUTGROUP>; close (OUTGROUP);
# A wgd Genome genome
open WGDGENOME, $subGenomeFile or die $!;
my @wgdGenome = <WGDGENOME>; close (WGDGENOME);
# Orthologs
open ORTHOLOGS, $orthologyFile or die $!;
my @orthologs = <ORTHOLOGS>; close (ORTHOLOGS);
# Output file
open OUTFILE, ">$outFilesDir\/$orthologRelationFile" or die $!;
open PPCHRGN, ">$outFilesDir\/$polyploidChrGenes" or die $!;
open OGCHRGN, ">$outFilesDir\/$outgroupChrGenes" or die $!;
#--------------------------------------------------------------#



# shift header lines
shift @wgdGenome;
shift @outgroup;
shift @orthologs;


my %Outgroup; my %WGDGenome; my %Orthologs; my %OutgroupPositions; my %WGDGenomePositions;

my $wgdGenomeCount = 1;
for (my $i = 0; $i < scalar (@wgdGenome); $i++){

	# Read two lines of the file at a time
	my @currentLine = split "\t", $wgdGenome[$i];		
	my @previousLine = split "\t", $wgdGenome[$i-1];
	map {$_=~s/\n|^\s+|\s+$//g} @currentLine;
	map {$_=~s/\n|^\s+|\s+$//g} @previousLine;

	#print "$currentLine[2]\t$previousLine[2]\t$wgdGenomeCount\n";

	# If the chromosome number is different then reset the count to 1
	if ($currentLine[2] ne $previousLine[2]){
		$wgdGenomeCount = 1;
	}

	$WGDGenome{$currentLine[0]} = \@currentLine; # Hash holding the gene ids and the corresponding lines
	$WGDGenomePositions{$currentLine[2]}{$currentLine[0]} = $wgdGenomeCount; # Nested hash holding chromosome, Id and its position

	$wgdGenomeCount++;
}

my $refCount = 1;
for (my $i = 0; $i< scalar (@outgroup); $i++){

	# Read two lines of the file at a time
	my @currentLine = split "\t", $outgroup[$i];
	my @previousLine = split "\t", $outgroup[$i-1];
	map {$_=~s/\n|^\s+|\s+$//g} @currentLine;
	map {$_=~s/\n|^\s+|\s+$//g} @previousLine;

	# If the chromosome number is different then reset the count to 1
	if ($currentLine[2] ne $previousLine[2]){
		$refCount = 1;
	}

	$Outgroup{$currentLine[0]} = \@currentLine;  # Hash holding the gene ids and the corresponding lines
	$OutgroupPositions{$currentLine[2]}{$currentLine[0]} = $refCount; # Nested hash holding chromosome, Id and its position

	$refCount++;
}

#print scalar keys %Outgroup,"\n";
#print join " ", keys %WGDGenome,"\n";
#print join " ", keys %Outgroup,"\n";
#print join "\n",keys %{$WGDGenomePositions{1}};
#print $WGDGenomePositions{'Y'}{'ENSG00000172283'};

# Print WGD chromosomes and gene counts in a file
foreach (keys %WGDGenomePositions){
	
	print PPCHRGN "$_\t";
	print PPCHRGN scalar keys %{$WGDGenomePositions{$_}},"\n";
}
close (PPCHRGN);

# Print outgroup chromosomes and gene counts in a file
foreach (keys %OutgroupPositions){
	
	print OGCHRGN "$_\t";
	print OGCHRGN scalar keys %{$OutgroupPositions{$_}},"\n";
}
close (OGCHRGN);

my %Duplicates; # This will mark and remove duplicate homolog pairs if any from the ortholog/paralog list
foreach (@orthologs){

	my @line = split "\t", $_;
	map {$_=~s/\n|^\s+|\s+$//g} @line;
	
	# Because user can place outgroup or genome as 1st or 2nd column, anywhere I am checking both ways.
	# Also, here the genes which do not exist in outgroup genome but are there in ortholog file(if any), will be ignored. This step is also problematic if a gene is ortholog with 2 genes and one of them is not there. 
	# Remember - the duplicates will be removed here
	
	if (not exists $Duplicates{$line[0]}{$line[1]}){
	
		if (exists $Outgroup{$line[0]}){ 
			push @{$Orthologs{$line[0]}}, $line[1];
		}
		if (exists $Outgroup{$line[1]}){
			push @{$Orthologs{$line[1]}}, $line[0];
		}
		# Mark duplicate pairs in both the directions
		$Duplicates{$line[0]}{$line[1]} = '';
		$Duplicates{$line[1]}{$line[0]} = '';
	}
}

print OUTFILE "$refGenomeName Chromosome\t$refGenomeName Position\t$subGenomeName Chromosome\t$subGenomeName Position\t$refGenomeName Gene Id\t$subGenomeName Gene Id\t$refGenomeName Orientation\t$subGenomeName Orientation\n";

foreach my $RefGene (keys %Outgroup){

	if (exists $Orthologs{$RefGene}){

		my @orthologLines = @{$Orthologs{$RefGene}};
		foreach (@orthologLines){

			my $refChr = ${$Outgroup{$RefGene}}[2];
			my $refPos = $OutgroupPositions{$refChr}{$RefGene};
			my $refGeneId = ${$Outgroup{$RefGene}}[0];
			my $refGeneOr = ${$Outgroup{$RefGene}}[1];
			if (exists $WGDGenome{$_}){  # here the genes which do not exist in polyploid genome (if any) will be ignored
				my $wgdGenomeChr = ${$WGDGenome{$_}}[2];
				my $wgdGenomeGeneId = ${$WGDGenome{$_}}[0];
				my $wgdGenomeOr = ${$WGDGenome{$_}}[1];
				my $wgdGenomePos = $WGDGenomePositions{$wgdGenomeChr}{$_};

				print OUTFILE "$refChr\t$refPos\t$wgdGenomeChr\t$wgdGenomePos\t$refGeneId\t$wgdGenomeGeneId\t$refGeneOr\t$wgdGenomeOr\n";
			}
			else {print OUTFILE "$refChr\t$refPos\t\t\t$refGeneId\t\t$refGeneOr\t\n";}
		}
	}
	else {

		my $refChr = ${$Outgroup{$RefGene}}[2];
		my $refPos = $OutgroupPositions{$refChr}{$RefGene};
		my $refGeneId = ${$Outgroup{$RefGene}}[0];
		my $refGeneOr = ${$Outgroup{$RefGene}}[1];
		print OUTFILE "$refChr\t$refPos\t\t\t$refGeneId\t\t$refGeneOr\t\n";
	}
}


print "................[DONE] in ";
if (eval (time - $startTime) >= 3600){printf "%.2f", eval(time - $startTime)/3600; print " hours\n";}
elsif (eval (time - $startTime) >= 60){printf "%.2f", eval(time - $startTime)/60; print " minutes\n";}
else {print time - $startTime, " seconds\n";}



