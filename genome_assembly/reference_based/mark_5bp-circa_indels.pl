use strict ;
use warnings ;

# filter variant records 5 bp around indels
# Find consensus indels (in all samples) - print total
# Filter snps in vcfs 5 bp +/- consensus indels

## usage: perl mark_5bp-circa_indels.pl Sv_MS_stampy_UGploidy1_ASboth.vcf

my$variants = $ARGV[0] ;

(my$indel_count, my$filtered_count) = filter_near_indels_vcf($variants) ;
print "total number of indels: ", $indel_count, "\n", "snps filtered within 5 bp of indel: " . $filtered_count . "\n" ;

sub filter_near_indels_vcf {
    my$vcf = shift ;
    open VCF, "<$vcf" or die "Cannot open file" ;
    my@VCF = <VCF> ;
    close VCF ;

    my%indels ;
    my$indel_count = 0 ;

## make list of indels
    foreach (@VCF) {

        if ( $_ =~ m/^\#/) {
            next ;
        }

        my@split = split ( /\t/, $_ ) ;

        #Skip invariant sites
        if ($split[4] eq ".") {
            next ;
        }

        my@split2 = split (/,/, $split[4]) ;

        if (length($split[3]) ne length($split2[0])) {
            $indel_count ++ ;
            if (exists $indels{$split[0]}) {
                push @{$indels{$split[0]}}, $split[1] ;
            }
            else {
                @{$indels{$split[0]}} = $split[1] ;
            }
        }
    }

    foreach my$string (keys %indels) {
            print "indel positions on $string: @{$indels{$string}}\n" ;
    }

    my$out = $vcf ;
    $out =~ s/vcf/filtered.vcf/ ;
    open OUT, ">$out" ;

    my$count_sites = 0 ;

    foreach (@VCF) {
        chomp $_ ;

        if ( $_ =~ m/^\#/ ) {
            print OUT $_, "\n" ;
            next ;
        }

        my@split = split ( /\t/, $_ ) ;

        if ($split[4] eq ".") {
            print OUT $_, "\n" ;
            next ;
        }

        else {
            my$report = "" ;

            # skip indel site itself
            if ($split[1] ~~ @{$indels{$split[0]}}) {
                $report = "false" ;
            }

            else {
                my@five_bp_around = $split[1]-5 .. $split[1]+5 ;

                foreach my$circa (@five_bp_around) {
                    if ($circa ~~ @{$indels{$split[0]}}) {
                        $report = "true" ;
                    }
                }
            }

            if ($report eq "true") {
               $count_sites ++ ;
                $_ =~ s/PASS/IF/ ;
                print OUT $_, "\n" ;
            }

            else {
                print OUT $_, "\n" ;
            }
        }
    }

    close OUT ;
    return ($indel_count, $count_sites) ;

}