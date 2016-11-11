use strict ;
use warnings ;
use Sort::Naturally ;

## Usage: perl ID_stats.pl < Sv_ALLsample-ID_stats.txt > Sv_ALLsample-ID_stats-ordered.txt

my@r2 ;
my@pval ;

while (<STDIN>) {

    if ($_ =~ m/^#/) {
        next ;
    }

    chomp ;

    my@split = split("\t", $_) ;

    push @r2, $split[14] ;
    push @pval, $split[16] ;

}
my$sum = 0 ;
foreach my$rsq (@r2) {
    $sum += $rsq ;
}

my$avg = $sum / @r2 ;

my$diff ;

foreach my$rsq (@r2) {
    $diff += ($rsq - $avg)**2 ;
}

my$variance = $diff / (@r2-1) ;
my$stdev = sqrt($variance) ;
my$stderr = $stdev / sqrt(@r2) ;

my$nonsig = 0 ;
my$sig = 0 ;

foreach my$prob (@pval) {
    if ($prob > 0.05) {
        $nonsig ++ ;
    }

    else {
        $sig ++ ;
    }
}

my$frac_nonsig = $nonsig / @pval ;
my$frac_sig = $sig / @pval ;

print "AverageID:", $avg,  "\tStdev:", $stdev, "\tStderr:", $stderr, "\tsignificant_fraction:", $frac_sig, "\tnonsignificant_fraction:", $frac_nonsig, "\n" ;

