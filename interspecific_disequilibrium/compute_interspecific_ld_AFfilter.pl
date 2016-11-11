use strict ;
use warnings ;
use Statistics::Distributions qw< chisqrprob > ;

my$min_AN = $ARGV[0] ;
my %geno ;

## retrieve qualified alleles and sites
while (<STDIN>) {

    #### skip headers
    if ( $_ =~ m/^#/ ) {
        next ;
    }

    my @split = split ( /\t/, $_ ) ;

    ## exclude indels and tri-allelic sites
    if ( length($split[3]) > 1 || length($split[4]) > 1 ) {
        next ;
    }

    ## exclude sites with AF < 0.1
    my @siteinfo = split ( /;/, $split[7] ) ;

    if ( $siteinfo[1] =~ m/(\d+\.\d+)/ ) {
        if ($1 <= 0.100) {
            next;
        }
    }

    # exclude sites not present in > AN=$min_AN samples (~50% of samples)
    if ( $siteinfo[2] =~ m/(\d+)/ ) {
        if ($1 < $min_AN) {
            next ;
        }
    }

    foreach my $index ( 9..$#split ) {
        if ( $split[$index] =~ m/^(\d):/ ) {
            $geno{"$split[0]\t$split[1]"}{$index} = $1 ;
        }
    }
}

## if more than one sample call exists at a site, append that site to @mito_positions or @sym_positions depending on source
my@mito_positions ;
my@sym_positions ;
foreach my$pos ( keys %geno ) {
    if ( scalar keys %{ $geno{$pos} } > 1 ) {
        my@coords = split(/\t/, $pos) ;
        my$genome = $coords[0] ;

        if ($genome eq "Sv_mito_chromosome") {
            push @mito_positions, $pos ;
            }

        if ($genome =~ m/Sv_sym_scaffold/) {
            push @sym_positions, $pos ;
            }
        }
    }

# Print header for output
print "sym_position\tmito_position\tsample_total\tavg_00_freq\t00_sum\t01_sum\t10_sum\t11_sum\tfreqA\tfreqB\tD\tr\tr^2\tchi-square_value\tp-value\n" ;

# iterate through all pairwise combinations of positions
foreach my $position1 ( 0..$#sym_positions ) {
    foreach my $position2 ( 0..$#mito_positions ) {

        my %genotype ;
        my $total = 0 ;
        $genotype{"00"} = 0 ;
        $genotype{"01"} = 0 ;
        $genotype{"10"} = 0 ;
        $genotype{"11"} = 0 ;

        ## iterate through samples, if sample has a call at position 1 and 2, add count to genotype combo type and total count
        foreach my $ind ( keys %{ $geno{$sym_positions[$position1]} } ) {
            if ( exists( $geno{$sym_positions[$position1]}{$ind} ) && exists ( $geno{$mito_positions[$position2]}{$ind} ) ) {
                $genotype{"$geno{$sym_positions[$position1]}{$ind}"."$geno{$mito_positions[$position2]}{$ind}"} ++ ;
                $total ++ ;
            }
        }

        ## skip remaining steps for this pair of sites if no sample has calls at both positions
        if ( $total == 0 ) {
            next ;
        }

        ## calculate genotype combo frequency, and frequences for each allele independently
        my $x00 = $genotype{"00"} / $total ;
        my $freqA = ( $genotype{"00"} + $genotype{"01"} ) / $total ;
        my $freqB = ( $genotype{"00"} + $genotype{"10"} ) / $total ;

        ## skip if either allele at locus A or B is fixed
        if ( $freqA == 0 || $freqA == 1 || $freqB == 0 || $freqB == 1 ) {
            next ;
        }

        ## calculate d, the deviation of the observed haplotype from its expected quantity
        my $d = $x00 - $freqA * $freqB ;
        ## calculate r, d corrected for allele frequency (so independent of it)
        my $sq = sqrt($freqA*(1-$freqA)*$freqB*(1-$freqB)) ;
        my $r = "NA" ;
        if ( $sq > 0 ) {
            $r = $d/$sq ;
        }

        ## calculate chi squared p-value with df = 4 classes - 2 parameters estimated (freqA, freqB) - 1 = 1
        my$r2 = "NA" ;
        my$chisq = "NA" ;
        my$pval = "NA" ;
        my$df = 1 ;

        if ($r ne "NA") {
            $r2 = $r * $r ;
            $chisq = $r2 * $total ;
            $pval = chisqrprob($df, $chisq) ;
        }


        print $sym_positions[$position1], "\t", $mito_positions[$position2], "\t", $total, "\t", $x00, "\t", $genotype{"00"}, "\t", $genotype{"01"}, "\t", $genotype{"10"}, "\t", $genotype{"11"},  "\t", $freqA, "\t", $freqB, "\t", $d, "\t", $r, "\t", $r2, "\t", $chisq, "\t", $pval, "\n" ;
    }
}
