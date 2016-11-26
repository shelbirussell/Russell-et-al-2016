use strict ;
use warnings ;

my$sum = 0 ;
my$lines = 0 ;

while(<STDIN>) {
        chomp ;
        my@split = split(/\t/, $_) ;

        $sum += $split[2] ;
        $lines ++ ;
}


print "average coverage: ", $sum/$lines, "\n" ;