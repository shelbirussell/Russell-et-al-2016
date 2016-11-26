use strict ;
use warnings ;

my$blast = $ARGV[0] ;

my$sample = $blast ;
$sample =~ s/(.+)hits.txt/$1/ ;
my$output = $sample . "filtered_hits.txt" ;

open BLAST, "<$blast" ;

my%hits ;

while (<BLAST>) {
    if ($_ =~ m/^#/) {
        next ;
        }

    chomp ;

    my@hit = split(/\t/, $_) ;

    my$qseqid = $hit[0] ;
    my$sallacc = $hit[1] ;
    my$salltitles = $hit[2] ;
    my$query_coverage = $hit[3] ;
    my$pident = $hit[4] ;
    my$length = $hit[5] ;
    my$mismatch = $hit[6] ;
    my$gapopen = $hit[7] ;
    my$qstart = $hit[8] ;
    my$qend = $hit[9] ;
    my$sstart = $hit[10] ;
    my$send = $hit[11] ;
    my$evalue = $hit[12] ;
    my$bitscore = $hit[13] ;

    # Note for some blast outputs: This matches up with the blast output, but not the requested output format - ssequid was specified and misspelled in some blasts, and results in no output (and no error)
    if (! $hits{$qseqid}) {
        if ($query_coverage > 50 && $pident > 30 && $evalue < 1e-6) {
            $hits{$qseqid}{"SUBJECT"} = $sallacc ;
            $hits{$qseqid}{"TITLES"} = $salltitles ;
            $hits{$qseqid}{"QCOV"} = $query_coverage ;
            $hits{$qseqid}{"PIDENTITY"} = $pident ;
            $hits{$qseqid}{"LENGTH"} = $length ;
            $hits{$qseqid}{"MISMATCH"} = $mismatch ;
            $hits{$qseqid}{"GAPS"} = $gapopen ;
            $hits{$qseqid}{"QSTART"} = $qstart ;
            $hits{$qseqid}{"QEND"} = $qend ;
            $hits{$qseqid}{"SSTART"} = $sstart ;
            $hits{$qseqid}{"SEND"} = $send ;
            $hits{$qseqid}{"Evalue"} = $evalue ;
            $hits{$qseqid}{"BITSCORE"} = $bitscore ;
        }
    }

    else {
        # keep best hit only - new evalue must be lower and alignment length less than 3/4 of the previous alignment to be kept
        if ($evalue < $hits{$qseqid}{"Evalue"} && $length >= $hits{$qseqid}{"LENGTH"}*0.75) {
            $hits{$qseqid}{"SUBJECT"} = $sallacc ;
            $hits{$qseqid}{"TITLES"} = $salltitles ;
            $hits{$qseqid}{"QCOV"} = $query_coverage ;
            $hits{$qseqid}{"PIDENTITY"} = $pident ;
            $hits{$qseqid}{"LENGTH"} = $length ;
            $hits{$qseqid}{"MISMATCH"} = $mismatch ;
            $hits{$qseqid}{"GAPS"} = $gapopen ;
            $hits{$qseqid}{"QSTART"} = $qstart ;
            $hits{$qseqid}{"QEND"} = $qend ;
            $hits{$qseqid}{"SSTART"} = $sstart ;
            $hits{$qseqid}{"SEND"} = $send ;
            $hits{$qseqid}{"Evalue"} = $evalue ;
            $hits{$qseqid}{"BITSCORE"} = $bitscore ;
        }

        else {
            next ;
        }
    }
}

close BLAST ;

open OUT, ">$output" ;

foreach my$query (sort {$a cmp $b} keys %hits) {
    print OUT $query, "\t", $hits{$query}{"SUBJECT"}, "\t", $hits{$query}{"TITLES"}, "\t", $hits{$query}{"QCOV"}, "\t", $hits{$query}{"PIDENTITY"}, "\t", $hits{$query}{"LENGTH"}, "\t", $hits{$query}{"MISMATCH"}, "\t", $hits{$query}{"GAPS"}, "\t", $hits{$query}{"QSTART"}, "\t", $hits{$query}{"QEND"}, "\t", $hits{$query}{"SSTART"}, "\t", $hits{$query}{"SEND"}, "\t", $hits{$query}{"Evalue"}, "\t", $hits{$query}{"BITSCORE"}, "\n" ;
}