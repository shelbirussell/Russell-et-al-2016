use strict ;
use warnings ;

my %bin ;
my %total ;

## Find distance between sites and r^2 for each pair located on the same scaffold, summing over r^2 at same distance
while (<STDIN>) {
    chomp $_ ;
    my @split = split ( /\t/, $_ ) ;

    if ($split[0] eq "Sv_mito_chromosome" || $split[2] eq "Sv_mito_chromosome") {
        next ;
        }

    if ( $split[0] eq $split[2] ) {
        $bin{abs($split[3]-$split[1])} += $split[13] * $split[13] ;
        $total{abs($split[3]-$split[1])} ++ ;
    }
}

## Sort distances in increasing order and print each distance and its average r^2 value
foreach my $dist ( sort {$a<=>$b} keys %bin ) {
    print $dist, "\t", $bin{$dist}/$total{$dist}, "\n" ;
}