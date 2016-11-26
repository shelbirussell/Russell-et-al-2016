use strict ;
use warnings ;
use Sort::Naturally ;

my$coverage_dir = $ARGV[0] ;

opendir(COV, $coverage_dir) ;
my@coverage_txt = readdir(COV) ;

foreach my$file (@coverage_txt) {
    my%coverage ;
    my$sample ;

    if ($file =~ m/(.+)MitoSym_coverage_by_site.txt/) {
        $sample = $1 ;
    }

    else {
        next ;
    }

    open IN, "<$coverage_dir$file" ;

    my$sum = 0 ;
    my$count = 0 ;

    my$window = 1000 ;
    my%scaffold_counter ;

    while (<IN>) {
        chomp ;

        my@split = split(/\t/, $_) ;

        if ($split[0] eq "Sv_mito_chromosome") {
            next ;
        }

        my$scaffold = $split[0] ;
        my$position = $split[1] ;
        my$coverage = $split[2] ;

        $sum += $coverage ;
        $count ++ ;

        if (! $scaffold_counter{$scaffold}{"COUNTER"}) {
            $scaffold_counter{$scaffold}{"COUNTER"} = 1 ;
            $scaffold_counter{$scaffold}{"WINDOW_SUM"} = $coverage ;
        }

        if ($scaffold_counter{$scaffold}{"COUNTER"} < 1000) {
            $scaffold_counter{$scaffold}{"COUNTER"} ++ ;
            $scaffold_counter{$scaffold}{"WINDOW_SUM"} += $coverage ;
        }

        else {
            $scaffold_counter{$scaffold}{"COUNTER"} ++ ;
            $scaffold_counter{$scaffold}{"WINDOW_SUM"} += $coverage ;

            my$window_average = $scaffold_counter{$scaffold}{"WINDOW_SUM"} / $scaffold_counter{$scaffold}{"COUNTER"} ;

            $coverage{$scaffold}{$position}{"START"} = $position - 999 ;
            $coverage{$scaffold}{$position}{"AVGCOV"} = $window_average ;

            $scaffold_counter{$scaffold}{"COUNTER"} = 0 ;
            $scaffold_counter{$scaffold}{"WINDOW_SUM"} = 0 ;
        }
    }

    close IN ;

    my$average = $sum/$count ;

    open OUT, ">${sample}_avg1000bp_circos_coverage.txt" ;

    foreach my$scaff (nsort keys %coverage) {
        foreach my$pos (nsort keys %{$coverage{$scaff}}) {

            my$scaled_cov = $coverage{$scaff}{$pos}{"AVGCOV"} / $average * 100 ;

            print OUT $scaff, "\t", $coverage{$scaff}{$pos}{"START"}, "\t", $pos, "\t", $scaled_cov, "\n" ;
        }
    }

    close OUT ;
}
