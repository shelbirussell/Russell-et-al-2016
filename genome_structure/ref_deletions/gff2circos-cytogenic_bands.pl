use strict ;
use warnings ;
use Sort::Naturally ;

my$annotation = $ARGV[0] ;
my$output = $ARGV[1] ;

my$recs = read_Svsym_gff($annotation) ;
my%recs = %{$recs} ;

my$scaff1 = 1213830 ;
my$scaff2 = 892554 ;
my$scaff3 = 537612 ;
my$scaff4 = 28016 ;

open OUT, ">$output" ;

my$prev_start = 0 ;

foreach my$scaff (nsort keys %recs) {
    foreach my$start (nsort keys %{$recs{$scaff}}) {
        if ($start-1 > $prev_start) {
            print OUT "band ", $scaff, " ", "chromosome", " ", "chromosome", " ", $prev_start, " ", $start-2, " black\n" ;
            $prev_start = $recs{$scaff}{$start}{"STOP"} ;
        }

        print OUT "band ", $scaff, " ", $recs{$scaff}{$start}{"ID"}, " ", $recs{$scaff}{$start}{"ID"}, " ", $start-1, " ", $recs{$scaff}{$start}{"STOP"}-1, " red\n" ;
    }

    if ($scaff eq "Sv_sym_scaffold1") {
        if ($prev_start < $scaff1) {
            print OUT "band ", $scaff, " ", "chromosome", " ", "chromosome", " ", $prev_start, " ", $scaff1, " black\n" ;
        }
    }

    elsif ($scaff eq "Sv_sym_scaffold2") {
        if ($prev_start < $scaff2) {
            print OUT "band ", $scaff, " ", "chromosome", " ", "chromosome", " ", $prev_start, " ", $scaff2, " black\n" ;
        }
    }

    elsif ($scaff eq "Sv_sym_scaffold3") {
        if ($prev_start < $scaff3) {
            print OUT "band ", $scaff, " ", "chromosome", " ", "chromosome", " ", $prev_start, " ", $scaff3, " black\n" ;
        }
    }

    else {
        print OUT "band ", $scaff, " ", "chromosome", " ", "chromosome", " ", $prev_start, " ", $scaff4, " black\n" ;
    }

    $prev_start = 0 ;
}

sub read_Svsym_gff {
    my$gff = $_[0] ;
	open REF_GFF, "<$gff" ;

    my %gff ;
	my%counts ;

    while (<REF_GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;
        my$scaffold = $split[0] ;

        if ( $scaffold eq "Sv_mito_chromosome" ) { next ; }

		if ($counts{$scaffold}) {$counts{$scaffold} ++ ;}
		
		else {$counts{$scaffold} = 1 ;}

        my $id = "" ;
        if ( $split[8] =~ m/locus_tag=(.+)(;product|;note)/ ) {
            $id = $1 ;
        }
        if ( $split[8] =~ m/locus_tag=(.+);note/ ) {
            $id = $1 ;
        }

        my $gene = "" ;
        if ( $split[8] =~ m/product=([\(\)a\/,a-zA-Z0-9_-]+)/ ) {
            $gene = $1 ;
        }

        $gff{$scaffold}{$split[3]}{"COUNT"} = $counts{$scaffold} ;
        $gff{$scaffold}{$split[3]}{"ID"} = $id ;
        $gff{$scaffold}{$split[3]}{"STOP"} = $split[4] ;
        $gff{$scaffold}{$split[3]}{"STRAND"} = $split[6] ;
        $gff{$scaffold}{$split[3]}{"PRODUCT"} = $gene ;

       }
    close REF_GFF ;
    return \%gff ;
    }