use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

my$cds = $ARGV[0] ;
my$rrna = $ARGV[1] ;
my$trna = $ARGV[2] ;

my$sample = basename($cds) ;
$sample =~ s/_.+$// ;

my%recs ;

open rRNA, "<$rrna" or die "can't open rRNA file\n" ;
my$rrna_count = 0 ;

while (<rRNA>) {
    if ($_ =~ /^#/) {
        next ;
    }

    chomp $_ ;
    my@split = split(/\t/, $_) ;
    my$scaffold = $split[0] ;
    my$source = $split[1] ;
    my$type = "rRNA" ;
    my$start = $split[3] ;
    my$end = $split[4] ;
    my$score = $split[5] ;
    my$strand = $split[6] ;
    my$frame = $split[7] ;
    my$attribute = $split[8] ;
    my$note = "na" ;

    my$out_string = "$scaffold\t$source\t$type\t$start\t$end\t$score\t$strand\t$frame\tID=$scaffold.$start.$attribute;Name=$attribute;Dbxref=na;Note=$note" ;
    $recs{$scaffold}{$start} = $out_string ;
    $rrna_count ++ ;
}

close rRNA ;

open tRNA, "<$trna" or die "can't open tRNA file\n" ;
my$trna_count = 0 ;

while (<tRNA>) {
    if ($_ =~ m/^scaffold.+/) {
        chomp $_ ;
        my@split = split(/\t/, $_) ;
        my$scaffold = $split[0] ;
        $scaffold =~ s/ // ;
        my$start = $split[2] ;
        $start =~ s/ // ;
        my$end = $split[3] ;
        $end =~ s/ // ;
        my$attribute = $split[4] ;
        my$anticodon = $split[5] ;
        my$score = $split[10] ;
        my$type = "tRNA" ;
        my$source = "tRNAscan" ;
        my$frame = 0 ;
        my$note = "anticodon:" . $anticodon ;
        my$strand ;

        if ($end-$start > 0) {
            $strand = "+" ;
        }

        if($end-$start < 0) {
            $strand = "-" ;
        }

    my$out_string = "$scaffold\t$source\t$type\t$start\t$end\t$score\t$strand\t$frame\tID=$scaffold.$start.$attribute;Name=$attribute-tRNA;Dbxref=na;Note=$note" ;
    $recs{$scaffold}{$start} = $out_string ;
    $trna_count ++ ;
    }
}
close tRNA ;

open CDS, "<$cds" ;
open OUT, ">${sample}_sym_genome.gff" ;

my@scaffolds ;
my$cds_count = 0 ;

while (<CDS>) {

    chomp $_ ;

    if ($_ =~ m/^#/) {
        if ($_ =~ m/^##gff-version/) {
            print OUT $_, "\n" ;
        }

        if ($_ =~ m/^#(scaffold_\d+)/) {
            my$scaffold = $1 ;
            my$start = 0 ;
            $recs{$scaffold}{$start} = $_ ;
        }

        if ($_ =~ m/^# Sequence Data: seqnum=\d+;seqlen=\d+;seqhdr="(scaffold_\d+)/) {
            my$scaffold = $1 ;
            my$start = 0 ;
            $recs{$scaffold}{$start} = $_ ;
        }

        next ;
    }

    my@split = split(/\t/, $_) ;
    my$scaffold = $split[0] ;
    my$start = $split[3] ;

    if (! grep (/$scaffold/, @scaffolds)) {
        push @scaffolds, $scaffold ;
    }

    $recs{$scaffold}{$start} = $_ ;
    $cds_count ++ ;
}

close CDS ;

my$rec_count = 0 ;

foreach my$scaffold (nsort keys %recs) {
    if (! exists $recs{$scaffold}{0}) {
        next ;
    }

    foreach my$start (nsort keys %{$recs{$scaffold}}) {
        print OUT $recs{$scaffold}{$start}, "\n" ;
        if ($start > 0) {
            $rec_count ++ ;
        }
    }
}

close OUT ;

print "rRNAs: ", $rrna_count, "\n" ;
print "tRNAs: ", $trna_count, "\n" ;
print "CDs: ", $cds_count, "\n" ;
print "total input records: ", $rrna_count + $trna_count + $cds_count, "\n" ;
print "records printed to genome.gff: ", $rec_count, "\n" ;