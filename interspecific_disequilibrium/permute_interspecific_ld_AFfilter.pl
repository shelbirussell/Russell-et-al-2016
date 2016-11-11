use strict ;
use warnings ;
use Statistics::Distributions qw< chisqrprob > ;

my$min_AN = $ARGV[0] ;
my$cutoff_r2 = $ARGV[1] ;
my$cutoff_sig = $ARGV[2] ;
my$prefix = $ARGV[3] ;
my %geno ;

print "parsing vcf file\n" ;

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

print "making list of positions with more than one call at a site\n" ;

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
#print "sym_position\tmito_position\tsample_total\tavg_00_freq\t00_sum\t01_sum\t10_sum\t11_sum\tfreqA\tfreqB\tD\tr\tr^2\tchi-square_value\tp-value\n" ;

print "starting permutations\n" ;

my$greater_than_r2 = 0 ;
my$greater_than_sig = 0 ;

open OUT, ">${prefix}_permutation.out" ;

foreach (0 .. 4999) {
	# Collect counts and averages
	my$permutation = $_ ;
	my$count = 0 ;
	my$sum_r2 = 0 ;
	my$sig_count = 0 ;
	my$insig_count = 0 ;

	# iterate through all pairwise combinations of positions
	foreach my $position1 ( 0..$#sym_positions ) {
    	foreach my $position2 ( 0..$#mito_positions ) {
	        my %genotype ;
     	   	my $total = 0 ;
			$genotype{"00"} = 0 ;
        	$genotype{"01"} = 0 ;
        	$genotype{"10"} = 0 ;
        	$genotype{"11"} = 0 ;

	        ## iterate through samples, randomize symbiont samples, if sample has a call at position 1 and 2, add count to genotype combo type and total count
    	    foreach my $ind ( keys %{ $geno{$sym_positions[$position1]} } ) {
				my@indices = keys %{ $geno{$sym_positions[$position1]} } ;
				my$random_mito = $indices[rand @indices] ;
            	if ( exists( $geno{$sym_positions[$position1]}{$ind} ) && exists ( $geno{$mito_positions[$position2]}{$random_mito} ) ) {
                	$genotype{"$geno{$sym_positions[$position1]}{$ind}"."$geno{$mito_positions[$position2]}{$random_mito}"} ++ ;
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

				$count ++ ;
				$sum_r2 += $r2 ;
		
				if ($pval < 0.05) {
					$sig_count ++ ;
				}
		
				else {
					$insig_count ++ ;
				}

        	}
#	        print $sym_positions[$position1], "\t", $mito_positions[$position2], "\t", $total, "\t", $x00, "\t", $genotype{"00"}, "\t", $genotype{"01"}, "\t", $genotype{"10"}, "\t", $genotype{"11"},  "\t", $freqA, "\t", $freqB, "\t", $d, "\t", $r, "\t", $r2, "\t", $chisq, "\t", $pval, "\n" ;
    	}
	}
	
	my$avg_r2 = $sum_r2 / $count ;

	if ($avg_r2 >= $cutoff_r2) {
		$greater_than_r2 ++ ;
	}
	
	my$total = $sig_count + $insig_count ;
	
	if ($sig_count/$total >= $cutoff_sig) {
		$greater_than_sig ++ ;
	} 
	print OUT $permutation, "\t", $count, "\t", $avg_r2, "\t", $sig_count, "\t", $insig_count, "\n" ;
}

close OUT ;

open OUT2, ">${prefix}_permutaions.results" ;

print OUT2 "Number of permutations with an average r2 greater than $cutoff_r2: ", $greater_than_r2, "\n" ;
print OUT2 "pvalue = ", ($greater_than_r2 + 1)/(5000 + 1), "\n" ;
print OUT2 "Number of permuations with the number of sites in significant LD greater than $cutoff_sig: ", $greater_than_sig, "\n" ;
print OUT2 "pvalue = ", ($greater_than_sig + 1)/(5000 + 1), "\n" ;

close OUT2 ;