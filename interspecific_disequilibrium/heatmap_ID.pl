use strict ;
use warnings ;
use Sort::Naturally ;

## Usage: perl heatmap_ID.pl ID_download/Sv_ALLsample-ID_stats-heatmap.txt < ID_download/Sv_ALLsample-ID_stats.txt

my$heatmap = $ARGV[0] ;

my%heatmap_table ;
my%test_table ;
my@mito_columns ;

while (<STDIN>) {

    if ($_ =~ m/^#/) {
        next ;
    }

    chomp ;

    my@split = split("\t", $_) ;

    $heatmap_table{$split[0] . "." . $split[1]}{$split[2] . "." . $split[3]} = $split[14] ;

	my$mito_coord = $split[2] . "." . $split[3] ;
    if (grep(/$mito_coord/, @mito_columns)) {
        next ;
    }

    else {
        push @mito_columns, $split[2] . "." . $split[3] ;
    }
}

my@sorted_mito = nsort @mito_columns ;

print "done reading in ID stats\n" ;

print "print out ID for all coordinates\n" ;

open OUT, ">$heatmap" or die "couldn't open $heatmap\n" ;

print OUT "symbiont_locus\t" ;

foreach (@sorted_mito) {
	print OUT $_, "\t" ;
}

print OUT "\n" ;

foreach my$sym (nsort keys %heatmap_table) {
	print OUT $sym, "\t" ;

		foreach my$mito (@sorted_mito) {				
			if ($heatmap_table{$sym}{$mito}) {
				print OUT $heatmap_table{$sym}{$mito}, "\t" ;
			}
	
			else {
				print OUT "NA\t" ;
			}
		}
	print OUT "\n" ;
}


close OUT ;


