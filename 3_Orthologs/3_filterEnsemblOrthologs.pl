# Filter the Ensmebl orthologs to remove blank lines and NAs
use strict;
use warnings;

# Sort them in different directories based on 2R or 3R
# The directories must be pre made or you'll get an error
my $wgd = 'Vertebrate'; # '2R', '3R', or 'Vertebrate'

foreach (<$wgd\_Orthologs\/*.txt>){
	
		print "$_\n";
		$_ =~/.+\/(.+)_(.+)_orthologs_.+.txt/g;
		my $org1 = $1; my $org2 = $2;
		#print "$org1\t$org2\n";
		
		open FH, "$_" or die $!;
		my @orth = <FH>;
		shift @orth;
		
		# open outfile
		open OUT, ">$wgd\_Orthologs_filtered\/$org1\_$org2\_orthologs_filtered_Ens84.txt" or die $!;
		print OUT "$org1\t$org2\n";
		
		foreach (@orth){
				
				my @line = split "\t", $_;
				map {$_=~s/\n//g} @line;
				
				if ($line[0] ne '' && $line[1] ne '' && $line[0] ne 'NA' && $line[1] ne 'NA'){
					
					print OUT "$line[0]\t$line[1]\n";	
				}
		}
		close (OUT);
}
