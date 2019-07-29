# Master script for outgroup comparison.
#
use strict;
use warnings;
use Getopt::Long;
use diagnostics; 

#my $inputFileDir = '../2_ParameterFiles/hsapiens/outgroup'; # Path to input files folder
#my $tempOutFilesDir = '../3_SyntenyOutputFiles/hsapiens/outgroup';   # Path to temporary outfiles folder 

if (scalar @ARGV < 2){die "Usage: script.pl input_file_relative_path output_file_relative_path\n";}

my $inputFileDir = $ARGV[0]; # Path to input files folder
my $tempOutFilesDir = $ARGV[1];   # Path to temporary outfiles folder
print "Input dir: $ARGV[0]\nOutput dir: $ARGV[1]\n\n";

print "Ohnolog run started at: ";
print `date`;
print "\n";

my $outgroupGenomeFile = "$inputFileDir\/OutgroupGenomes.in";
my $polyploidGenomeFile = "$inputFileDir\/PolyploidGenomes.in";
my $orthologFile = "$inputFileDir\/OrthologFiles.in";
my $parameterFile = "$inputFileDir\/RunParameters.in";
my $shrink = 'n';
my $probability = 1; # Probabiity cutoff to be used to filetr blocks. Blocks with P <= this cutoff will be filtered
my $self = 'y';

GetOptions ("outgroup=s"    => \$outgroupGenomeFile,
            "polyploid=s"   => \$polyploidGenomeFile,
            "orthologs=s"   => \$orthologFile,
            "parameters=s"  => \$parameterFile,
            "shrink=s" => \$shrink,
            "probability=f" => \$probability
);


#print "$outgroupGenomeFile\n$polyploidGenomeFile\n$OrthologFile\n$parameterFile\n";

open OUTGP, "$outgroupGenomeFile" or die $!;
my %outgroupGenomes = %{readGenomeFiles (\*OUTGP)};
#print join "\n", keys %outgroupGenomes;

open POLYP, "$polyploidGenomeFile" or die $!;
my %polyploidGenomes = %{readGenomeFiles (\*POLYP)};
#print join "\n", keys %polyploidGenomes;

open ORTH, "$orthologFile" or die $!;
my %orthologs = %{readParameterFiles(\*ORTH)};
#print join "\n", keys %orthologs;

open PARAM, "$parameterFile" or die $!;
my %parameters = %{readParameterFiles(\*PARAM)};
#print join "\n", keys %parameters;

foreach my $polyploidName (sort keys %polyploidGenomes){
	
	foreach my $outgroupName (sort keys %outgroupGenomes){
		
#		print "$outgroupName\t$outgroupGenomes{$outgroupName}\n";
#		print "$polyploidName\t$polyploidGenomes{$polyploidName}\n";
#		print "${$orthologs{$polyploidName}{$outgroupName}}[2]\t";
#		print join ("\t", @{$parameters{$polyploidName}{$outgroupName}}[2..4])."\n\n";

		my $outgroupFile = $outgroupGenomes{$outgroupName};
		my $polyploidFile = $polyploidGenomes{$polyploidName};
		my $orthologFile = ${$orthologs{$polyploidName}{$outgroupName}}[2];
		my $outgroupWindow = ${$parameters{$polyploidName}{$outgroupName}}[2];
		my $polyploidWindow = ${$parameters{$polyploidName}{$outgroupName}}[3];
		my $minRequiredOrthologs = ${$parameters{$polyploidName}{$outgroupName}}[4];

		print `date`, "\n";
		print `perl OrthologRelations.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $tempOutFilesDir`;
		print `date`, "\n";
		print `perl Blocks_OutgroupComparison_20121229.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $shrink $tempOutFilesDir`;		
		print `date`, "\n";
		print `perl filterBlocksWithMultipleOrthologPairs.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $tempOutFilesDir`;
		print `date`, "\n";

		print `perl CalculateProbability_EntireGenome_WithoutLog.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $shrink $tempOutFilesDir`;
		print `date`, "\n";
		
		print `perl getAnchorsBelongingtoMultipleBlocks_OutgroupComparison.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $shrink $probability $tempOutFilesDir`;
		print `date`, "\n";
		print `perl getPositionsForPlots.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $shrink $tempOutFilesDir`;
		print `date`, "\n";
		print `perl getOhnoPairs_OutgroupComparison.pl $outgroupName $outgroupFile $polyploidName $polyploidFile $orthologFile $outgroupWindow $polyploidWindow $minRequiredOrthologs $tempOutFilesDir`;
		print `date`, "\n";

	}
}


sub readGenomeFiles {
	
	my $FH = shift;
	my %genome;
	
	foreach (<$FH>){
		
		if ($_ !~/^[#\n\s+\t]/g){
			
			my @line = split "\t", $_;
			map {$_=~s/\n|^\s+|\s+$|\t+//g} @line;
			$genome{$line[0]} = $line[1];
		}
	}
	return (\%genome);
}

sub readParameterFiles {
	
	my $FH = shift;
	my %params;
	
	foreach (<$FH>){
		
		if ($_ !~/^[#\n\s+\t]/g){
					
			my @line = split "\t", $_;
			map {$_=~s/\n|\s+|\t+//g} @line;
			#print "$line[0]\t$line[1]\n";
			$params{$line[0]}{$line[1]} = \@line;
			$params{$line[1]}{$line[0]} = \@line;
		}
	}
	#foreach (keys %params){print "$_\n";}
	return (\%params);
}


print "\nAll done at: ";
print `date`;
print "\n";
