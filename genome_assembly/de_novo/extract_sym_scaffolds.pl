use strict ;
use warnings ;
use File::Basename ;

### This script parses scaffolds in to files based upon their length and coverage.  The built-in lower length bound is 1kb.
### perl extract_scaffolds.pl fasta coverage.txt symGC=% [output_directory]

my$fasta = $ARGV[0] ;
my$cov_file = $ARGV[1] ;
my$sym_GC = $ARGV[2] ;

my$output_dir ;
if ($ARGV[4]) {$output_dir = $ARGV[3] ;}
else {$output_dir = "./" ;}

my@sample_info = split("_", basename($fasta)) ;
my$sample = $sample_info[0] ;

my$seq = read_fasta($fasta) ;
my($cov, $cov_mean, $cov_stdev) = read_cov($cov_file) ;

my%seq = %{$seq} ;
my%cov = %{$cov} ;


# 95% confidence interval bounds for coverage (Z=1.96); 99% confidence interval bounds (Z=2.576)
my$lowerbound = $cov_mean - 1 * $cov_stdev ;
my$upperbound = $cov_mean + 1 * $cov_stdev ;


### Filter scaffolds by coverage
my%symbiont ;
my%mitochondrian ;

my$filtered_GC = 0 ;
my$filtered_cov = 0 ;

foreach my$header (keys %seq) {
    my$scaff_seq = $seq{$header} ;
    my$length = length($scaff_seq) ;
    my$scaff_cov = $cov{$header} ;

    if ($length < 1000) {
        next ;
    }

    my$Gs = () = $scaff_seq =~ /G/g ;
    my$Cs = () = $scaff_seq =~ /C/g ;
    my$gc = (($Gs + $Cs) / $length)*100 ;

    if ($gc < $sym_GC - 10 || $sym_GC + 10 < $gc) {
        next;
        $filtered_GC ++ ;
    }

    if ($lowerbound < $scaff_cov && $scaff_cov < $upperbound) {
        $symbiont{$header}{"SEQ"} = $scaff_seq ;
        $symbiont{$header}{"LENGTH"} = $length ;
        $symbiont{$header}{"GC"} = $gc ;
    }

    else {
        $filtered_cov ++ ;
    }
}

my$sym_outfile = $output_dir . $sample . "_extracted-" . $cov_mean . "x_" . $sym_GC . "%GC_sym_scaff.fasta" ;
open OUT, ">$sym_outfile" or print "Can't open file: ", $sym_outfile, "\n" ;

foreach my$header (sort {$symbiont{$b}{"LENGTH"} <=> $symbiont{$a}{"LENGTH"}} keys %symbiont) {
    print OUT ">", $header, "\t", $sample, "_symbiont\tlength:", $symbiont{$header}{"LENGTH"}, "\t%GC:", $symbiont{$header}{"GC"}, "\n" ;
    print OUT "$_\n" foreach ($symbiont{$header}{"SEQ"} =~ /.{1,80}/g) ;
}

close OUT ;

print "Filtered scaffolds:\n${filtered_GC} out of GC range\n${filtered_cov} out of coverage range\n" ;

sub read_cov {
    open COVERAGE, "<$_[0]" ;

    my%coverage ;
    my@depths = () ;

    while (<COVERAGE>) {

        if ( $_ =~ m/^#/ ) {next ;}

        my@split = split(/\t/, $_) ;

        my$scaff = $split[0] ;
        my$site = $split[1] ;
        my$depth = $split[2] ;

        $coverage{$scaff}{$site} = $depth ;

        push @depths, $depth ;
    }

    close COVERAGE ;

    my%avg_list ;

    foreach my$scaff (keys %coverage) {
        my$sites = 0 ;
        my$total_depth = 0 ;

        foreach my$site (keys %{$coverage{$scaff}}) {
            $sites ++ ;
            $total_depth = $total_depth + $coverage{$scaff}{$site} ;
        }

        my$avg_cov = $total_depth / $sites ;

        $avg_list{$scaff} = $avg_cov ;
    }

    ### Calculate depth statistics to filter on
	my$sum ;
	$sum += $_ for @depths ;
	my$mean = $sum / @depths ;

	my$dev ;
	foreach (@depths) {$dev += ($_ - $mean)**2 ;}
	my$stdev = sqrt($dev / @depths) ;

    return \%avg_list, \$mean, \$stdev ;

}

sub read_fasta {
    open FASTA, "<$_[0]" ;

    my%seqs ;
    my$header ;
    my$seq ;

    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
                $seqs{$header} = $seq ;
                }

            $header = $_ ;
            $header =~ s/^>// ;
            $header =~ s/\s+$// ;

            $seq = "" ;
            }

        else {
            $_ =~ s/\s+//g ;
            $seq .= $_ ;
            }
        }

    close FASTA ;

    if ($seq) {
        $seqs{$header} = $seq ;
        }

    return \%seqs ;
    }