use strict ;
use warnings ;

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

## if more than one sample call exists at a site, append that site to @positions
my @positions ;
foreach ( keys %geno ) {
    if ( scalar keys %{ $geno{$_} } > 1 ) {
        push @positions, $_ ;
    }
}

# iterate through all pairwise combinations of positions
foreach my $position1 ( 0..$#positions-1 ) {
    foreach my $position2 ( $position1+1..$#positions ) {

        my %genotype ;
        my $total = 0 ;			# samples with calls at locus A and B
        $genotype{"00"} = 0 ;	# ref calls at locus A and B
        $genotype{"01"} = 0 ;	# refA, altB 
        $genotype{"10"} = 0 ;	# altA, refB
        $genotype{"11"} = 0 ;	# alt calls at locus A and B

        ## iterate through samples, if sample has a call at position 1 and 2, add count to genotype combo type and total count
        foreach my $ind ( keys %{ $geno{$positions[$position1]} } ) {
            if ( exists( $geno{$positions[$position1]}{$ind} ) && exists ( $geno{$positions[$position2]}{$ind} ) ) {
                $genotype{"$geno{$positions[$position1]}{$ind}"."$geno{$positions[$position2]}{$ind}"} ++ ;
                $total ++ ;
            }
        }

        ## skip remaining steps for this pair of sites if no sample has calls at both positions
        if ( $total == 0 ) {
            next ;
        }

        ## calculate genotype combo frequency, and frequences for each allele independently
        my $x00 = $genotype{"00"} / $total ;
        my $freqA = ( $genotype{"00"} + $genotype{"01"} ) / $total ;	# calc ref allele frequency for A
        my $freqB = ( $genotype{"00"} + $genotype{"10"} ) / $total ;	# calc ref allele frequency for B

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

        print $positions[$position1], "\t", $positions[$position2], "\t", $total, "\t", $x00, "\t", $genotype{"00"}, "\t", $genotype{"01"}, "\t", $genotype{"10"}, "\t", $genotype{"11"},  "\t", $freqA, "\t", $freqB, "\t", $d, "\t", $r, "\n" ;
    }
}