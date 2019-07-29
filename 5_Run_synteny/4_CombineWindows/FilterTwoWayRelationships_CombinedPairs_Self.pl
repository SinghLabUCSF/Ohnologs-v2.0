use strict;
use warnings;

my $polyploidName = $ARGV[0];
my $wgd = $ARGV[1];

open PAIR, "$polyploidName\/$polyploidName-$polyploidName\_CombinedPairs_AllWindows$wgd\.txt" or die $_;
open OUT, ">$polyploidName\/$polyploidName\_Ohnologs_CombinedFromSelf_OneWay$wgd\.txt" or die $_;

my %Pairs;

foreach (<PAIR>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Pairs{$line[0].'~'.$line[1]} = \@line; # Make pairs in any orientation
}
print "Pairs without filtering: ", scalar keys %Pairs, "\n";

my %explored; # Hash to mark checked gene already

foreach (keys %Pairs){ # Foreach pair
	
	if (not exists $explored{$_}){ # If it has not been explored -  to avoid counting pairs twice
		
		my ($P1, $P2) = split '~', $_; # get the genes

		if (exists $Pairs{$P2.'~'.$P1}){ # if reverse is there
			
			print OUT join ("\t", @{$Pairs{$P1.'~'.$P2}}[0..1]), "\t"; # print the ids
			#print join ("\t", @{$Pairs{$P2.'~'.$P1}}), "\n--------------\n";
			
			for (my $i = 2; $i <=7; $i++){ # Foreach probability
			
				#print ${$Pairs{$P1.'~'.$P2}}[$i],"\t";
				#print ${$Pairs{$P2.'~'.$P1}}[$i],"\t";
				
				# Calculate and print the geometric mean for respective probability
				my $multi = abs (${$Pairs{$P1.'~'.$P2}}[$i]) * abs(${$Pairs{$P2.'~'.$P1}}[$i]);
				my $geoMean = $multi ** (1/2);
				print OUT "$geoMean\t";
			}
		}
		else { # If reverse pair is not there - die with error because it must be there
			die "Error at $P1 - $P2\n";
		}
		print OUT "\n";
		# And mark both these pairs explored so that we dnt count them again
		$explored{$P2.'~'.$P1} = '';
		$explored{$P1.'~'.$P2} = '';
	}
}


