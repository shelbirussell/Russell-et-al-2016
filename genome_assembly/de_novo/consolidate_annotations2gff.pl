use strict ;
use warnings ;
use File::Basename ;

## This script takes multiple blast results and consolidates them into one gff file, using another gff file as a backbone for coordinates

my$cds = $ARGV[0] ;
my$nr = $ARGV[1] ;
my$uniref90 = $ARGV[2] ;
my$trembl = $ARGV[3] ;

my$nr_hits = parse_blast($nr) ;
my$uniref_hits = parse_blast($uniref90) ;
my$trembl_hits = parse_blast($trembl) ;

my%nr_hits = %{$nr_hits} ;
my%uniref_hits = %{$uniref_hits} ;
my%trembl_hits = %{$trembl_hits} ;

my$sample = basename($cds) ;
$sample =~ s/_.+$// ;

open OUT, ">${sample}_sym_cds.gff" ;

open CDS, "<$cds" ;

while (<CDS>) {
    if ($_ =~ m/^#/) {
        print OUT $_ ;
        next ;
    }

    my@split = split(/\t/, $_) ;
    my$scaffold = $split[0] ;
    my$source = $split[1] ;
    my$type = $split[2] ;
    my$start = $split[3] ;
    my$end = $split[4] ;
    my$score = $split[5] ;
    my$strand = $split[6] ;
    my$phase = $split[7] ;

    my@info = split(";", $split[8]) ;
    my$id = $info[0] ;

    $id =~ s/ID=(scaffold_\d+_\d+)/$1/ ;

    my$gene = $id ;

    if ($nr_hits{$gene} || $uniref_hits{$gene} || $trembl_hits{$gene}) {
        my$evalue ;
        my$db ;
        my$hit ;
        my$product ;
        my$pident ;

        my$nr_ident = 0 ;
        my$uniref_ident = 0 ;
        my$trembl_ident = 0 ;
        if ($nr_hits{$gene}{"PIDENT"}) {
            $nr_ident = $nr_hits{$gene}{"PIDENT"} ;
        }
        if ($uniref_hits{$gene}{"PIDENT"}) {
            $uniref_ident = $uniref_hits{$gene}{"PIDENT"} ;
        }
        if ($trembl_hits{$gene}{"PIDENT"}) {
            $trembl_ident = $trembl_hits{$gene}{"PIDENT"} ;
        }

        my@idents = ($nr_ident, $uniref_ident, $trembl_ident) ;
        my$max_ident = 0 ;
        $idents[$max_ident] > $idents[$_] or $max_ident = $_ for 1 .. $#idents ;

        if ($max_ident == 0) {
            $evalue = $nr_hits{$gene}{"Evalue"} ;
            $db = "NCBI_nr" ;
            $hit = $nr_hits{$gene}{"HIT"} ;
            $product = $nr_hits{$gene}{"PRODUCT"} ;
            $pident = $nr_hits{$gene}{"PIDENT"} ;
        }

        elsif ($max_ident == 1) {
            $evalue = $uniref_hits{$gene}{"Evalue"} ;
            $db = "UniRef90" ;
            $hit = $uniref_hits{$gene}{"HIT"} ;
            $product = $uniref_hits{$gene}{"PRODUCT"} ;
            $pident = $uniref_hits{$gene}{"PIDENT"} ;
        }

        elsif ($max_ident == 2) {
            $evalue = $trembl_hits{$gene}{"Evalue"} ;
            $db = "TrEMBL" ;
            $hit = $trembl_hits{$gene}{"HIT"} ;
            $product = $trembl_hits{$gene}{"PRODUCT"} ;
            $pident = $trembl_hits{$gene}{"PIDENT"} ;
        }

        else {
            print "script not working\n" ;
            last ;
        }

        print OUT $scaffold, "\t", $source, "\t", $type, "\t", $start, "\t", $end, "\t", $evalue, "\t", $strand, "\t", $phase, "\t", "ID=", $gene, ";Name=", $product, ";Dbxref=", $db, ":", $hit, ";Note=percent_identity:", $pident, "\n" ;
    }

    else {
        print OUT $scaffold, "\t", $source, "\t", $type, "\t", $start, "\t", $end, "\t", $score, "\t", $strand, "\t", $phase, "\t", "ID=", $gene, ";Name=hypothetical_protein;Dbxref=na;Note=na\n" ;
    }
}

close OUT ;

sub parse_blast {
    my%hits ;
    open REPORT, "<$_[0]" ;

    while (<REPORT>) {
        if ($_ =~ m/^#/) {
             next ;
        }

        my@info = split(/\t/, $_) ;
        my$cds = $info[0] ;
        my$hit = $info[1] ;
        my$product = $info[2] ;
        $product =~ s/ /_/g ;
        my$pident = $info[3] ;
        my$evalue = $info[12] ;

        $hits{$cds}{"HIT"} = $hit ;
        $hits{$cds}{"PRODUCT"} = $product ;
        $hits{$cds}{"PIDENT"} = $pident ;
        $hits{$cds}{"Evalue"} = $evalue ;
    }

    close REPORT;

    return \%hits ;
}