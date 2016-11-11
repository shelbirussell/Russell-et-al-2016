use strict ;
use warnings ;
use Statistics::RankCorrelation ;

my@dist ;
my@r2 ;

while (<STDIN>) {
    my@split = split (/\t/, $_) ;
    push @dist, $split[0] ;
    push @r2, $split[1] ;
    }

my$cor = Statistics::RankCorrelation->new( \@dist, \@r2 ) ;

my$rho = $cor->spearman ;
print "Spearman rank correlation coefficient: ", $rho, "\n" ;
