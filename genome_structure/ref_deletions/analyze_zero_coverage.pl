use strict ;
use warnings ;
use Sort::Naturally ;
use Data::Dumper ;

## Usage: perl analyze_zero_coverage.pl zerocov_window annotation.gff tes.gff ./
## Ex. perl ../analyze_zero_coverage.pl 150 ~/Reference/Solemya_velum_MitoSym_top4.gff /n/regal/cavanaugh_lab/shelbirussell/TE_analysis/Svsym_mobile_elements.gff ../cov_files/

my$window = $ARGV[0] ;
my$annotation = $ARGV[1] ;
my$te_annotation = $ARGV[2] ;
my$coverage_dir = $ARGV[3] ;

## Load annotation
my$anno = read_Sv_gff($annotation) ;
my$tes = read_gff($te_annotation) ;
my%anno = %{$anno} ;
my%tes = %{$tes} ;

my$zero_cov = parse_coverage($window, $coverage_dir) ;
my%zero_cov = %{$zero_cov} ;

## Inventory zero coverage regions in the Sv sym genome
my%zero_regions ;
foreach my$sample (keys %zero_cov) {
    foreach my$scaff (keys %{$zero_cov{$sample}}) {

        if (! $zero_regions{$scaff}) {
            @{$zero_regions{$scaff}} = () ;
        }

        foreach my$span (@{$zero_cov{$sample}{$scaff}}) {
            my@split = split("-", $span) ;
            my@range = $split[0] .. $split[1] ;

            push @{$zero_regions{$scaff}}, @range ;
        }
    }
}

## Annotate zero coverage regions
my%annotated_regions ;
foreach my$scaff (keys %zero_regions) {
    my@filtered = uniq(@{$zero_regions{$scaff}}) ;
    my@sorted_filtered = nsort @filtered ;

    my$start = "" ;

    foreach my$index (0 .. $#sorted_filtered) {
        if ($start eq "") {
            $start = $sorted_filtered[$index] ;
        }

        # skip to next site if sites are in sequential order
        if ($sorted_filtered[$index] == $sorted_filtered[$index-1]+1) {
            next ;
        }

        # skip to next site if $index-1 loops around to the last index (ie $index is the first index)
        elsif ($sorted_filtered[$index] < $sorted_filtered[$index-1]) {
            next ;
        }

        # record range if contiguous zero coverage region ended
        else {
            my$range = $start . "-" . $sorted_filtered[$index-1] ;
            my@range = $start .. $sorted_filtered[$index-1] ;

            my@genes  = () ;
            foreach my$id (nsort keys %{$anno{$scaff}}) {
                if ($anno{$scaff}{$id}{"START"} ~~ @range || $anno{$scaff}{$id}{"STOP"} ~~ @range) {
                    my$gene = $anno{$scaff}{$id}{"GENE"} . ":" . $anno{$scaff}{$id}{"START"} . "-" . $anno{$scaff}{$id}{"STOP"} ;
                    push @genes, $gene ;
                }
            }

            my@tes= () ;
            foreach my$id (nsort keys %{$tes{$scaff}}) {
                if ($tes{$scaff}{$id}{"START"} ~~ @range || $tes{$scaff}{$id}{"STOP"} ~~ @range) {
                    my$te = $tes{$scaff}{$id}{"GENE"} . ":" . $tes{$scaff}{$id}{"START"} . "-" . $tes{$scaff}{$id}{"STOP"} ;
                    push @tes, $te ;
                }
            }

            $annotated_regions{$scaff}{$range}{"GENES"} = \@genes ;
            $annotated_regions{$scaff}{$range}{"TES"} = \@tes ;

            $start = $sorted_filtered[$index] ;
        }
    }
}

# Characterize each sample for each annotated zero-coverage region
my%sample_regions ;
foreach my$sample (keys %zero_cov) {
    foreach my$scaff (keys %{$zero_cov{$sample}}) {
        foreach my$span (@{$zero_cov{$sample}{$scaff}}) {
            my@sample_split = split("-", $span) ;

            foreach my$range (keys %{$annotated_regions{$scaff}}) {
                if (! $sample_regions{$scaff}{$sample}{$range}) {
                    @{$sample_regions{$scaff}{$sample}{$range}} = () ;
                }

                my@inventory_split = split("-", $range) ;

                if ($sample_split[0] >= $inventory_split[0] && $sample_split[1] <= $inventory_split[1]) {
                    push @{$sample_regions{$scaff}{$sample}{$range}}, $span ;
                }
            }
        }
    }
}

foreach my$scaff (nsort keys %annotated_regions) {
    print "##", $scaff, "\n" ;

    foreach my$region (nsort keys %{$annotated_regions{$scaff}}) {
        print "#region:", $region, "\tGenes:", join(",", @{$annotated_regions{$scaff}{$region}{"GENES"}}), "\n" ;
        print "#region:", $region, "\tTEs:", join(",", @{$annotated_regions{$scaff}{$region}{"TES"}}), "\n" ;
    }

    my@regions = nsort keys %{$annotated_regions{$scaff}} ;
    print "sample\t", join("\t", @regions), "\n" ;

    foreach my$sample (nsort keys %{$sample_regions{$scaff}}) {
        print $sample ;

        foreach my$region (nsort keys %{$annotated_regions{$scaff}}) {
            if (@{$sample_regions{$scaff}{$sample}{$region}}) {
                if (scalar @{$sample_regions{$scaff}{$sample}{$region}} > 1) {
                    print "\t-(", join(",", @{$sample_regions{$scaff}{$sample}{$region}}), ")" ;
                }

                else {
                    print "\t-" ;
                }
            }

            else {
                print "\t+" ;
            }
        }

        print "\n" ;
    }
}


sub read_Sv_gff {
    open GFF, "<$_[0]" ;

    my %data ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;

        if ( $split[0] eq "Sv_mito_chromosome" ) { next ; }

        if ( $split[2] ne "CDS" ) { next ; }

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

        $data{$split[0]}{$id}{"START"} = $split[3] ;
        $data{$split[0]}{$id}{"STOP"} = $split[4] ;
        $data{$split[0]}{$id}{"STRAND"} = $split[6] ;
        $data{$split[0]}{$id}{"GENE"} = $gene ;

       }
    close GFF ;
    return \%data ;
    }

sub read_gff {
    open GFF, "<$_[0]" ;

    my %data ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;

        if ( $split[0] eq "Sv_mito_chromosome" ) { next ; }

#        if ( $split[2] ne "CDS" ) { next ; }

        my $id = "" ;
        if ( $split[8] =~ m/ID=(.+);Name=/ ) {
            $id = $1 ;
        }

        my $gene = "" ;
        if ( $split[8] =~ m/;Name=(.+);Dbxref/ ) {
            $gene = $1 ;
        }

        $data{$split[0]}{$id}{"START"} = $split[3] ;
        $data{$split[0]}{$id}{"STOP"} = $split[4] ;
        $data{$split[0]}{$id}{"STRAND"} = $split[6] ;
        $data{$split[0]}{$id}{"GENE"} = $gene ;

       }
    close GFF ;
    return \%data ;
    }

sub parse_coverage {
    ## Find zero coverage regions for each sample
    my$span = $_[0] ;
    opendir(COV, $_[1]) ;
    my@coverage_txt = readdir(COV) ;

    my%zero_cov ;
    foreach my$file (@coverage_txt) {
        my$sample ;

        if ($file =~ m/(.+)MitoSym_coverage_by_site.txt/) {
            $sample = $1 ;
        }

        else {
            next ;
        }

        open IN, "<$coverage_dir$file" ;

        my%coverage ;
        while (<IN>) {
            chomp ;

            my@split = split(/\t/, $_) ;

            if ($split[0] eq "Sv_mito_chromosome") {
                next ;
            }

            my$scaffold = $split[0] ;
            my$position = $split[1] ;
            my$coverage = $split[2] ;

            if (! $coverage{$scaffold}) {
                @{$coverage{$scaffold}} = () ;
            }

            push @{$coverage{$scaffold}}, $coverage ;
        }

        ## Find contiguous regions of zero coverage (> $span bp long)
        foreach my$scaff (keys %coverage) {
            if (! $zero_cov{$sample}{$scaff}) {
                @{$zero_cov{$sample}{$scaff}} = () ;
            }

            my$count = 0 ;
            foreach my$index (0..$#{$coverage{$scaff}}) {
                if ($coverage{$scaff}[$index] > 1) {

                    if ($count >= $span) {
                        my$end = $index ;
                        my$start = $index - $count + 1 ;
                        push @{$zero_cov{$sample}{$scaff}}, $start . "-" . $end ;

#                        print $start, ":", $coverage{$scaff}[$start-1], "\t", $end, ":", $coverage{$scaff}[$end-1], "\n" ;
                    }

                    $count = 0 ;
                    next ;
                }

                else {
                    $count ++ ;
                }
            }
        }

        close IN ;
    }

    return\%zero_cov ;

}

sub uniq {
    my%seen ;
    grep !$seen{$_} ++, @_ ;
}