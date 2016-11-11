use strict ;
use warnings ;

# This script takes a phylip file produced from a mauve-generated xmfa and concatenates alignments for which at least @samples + 1 outgroup are represented
# Usage: perl parse_mauve_multiphylip.pl $num_outgroups $phylip

my$num_outgroups = $ARGV[0] ;
my$phylip = $ARGV[1] ;
$phylip =~ m/(.+)\.phy/ ;
my$out = $1 . "_" . $num_outgroups . "_concat.phy" ;

my($alignments, $samples) = read_phy($phylip) ;

my%alignments = %{$alignments} ;
my@samples = @{$samples} ;

my%concat_align ;
$concat_align{"SAMPLES"} = @samples ;
$concat_align{"LENGTH"} = 0 ;

# Concatenate alignments containing >= @samples + 1 outgroup
foreach my$number (sort {$a<=>$b} keys %alignments) {
    if ($alignments{$number}{"SAMPLES"} < @samples - $num_outgroups + 1 ) {
        print "Excluded: ", $number, "\t", $alignments{$number}{"SAMPLES"}, "\n" ;
        next ;
    }

    else {
        $concat_align{"LENGTH"} += $alignments{$number}{"LENGTH"} ;
        my@aligned_samples ;

        foreach my$sample (keys %{$alignments{$number}{"ALIGNMENT"}}) {
            my$sample_name = $alignments{$number}{"ALIGNMENT"}{$sample}{"SAMPLE"} ;

            push @aligned_samples, $sample_name ;

            $concat_align{"ALIGNMENT"}{$sample_name} .= $alignments{$number}{"ALIGNMENT"}{$sample}{"SEQ"} ;
        }

        foreach my$sample (@samples) {
            if ($sample ~~ @aligned_samples) {
                next ;
            }

            else {
                my$insert = "-" x $alignments{$number}{"LENGTH"} ;
                $concat_align{"ALIGNMENT"}{$sample} .= $insert ;
            }
        }
    }
}

open OUT, ">$out" ;
print OUT $concat_align{"SAMPLES"}, " ", $concat_align{"LENGTH"}, "\n" ;
foreach my$sample (@samples) {
    print OUT $sample, "  " ;
    print OUT " $_" foreach ($concat_align{"ALIGNMENT"}{$sample} =~ /.{1,10}/g), "\n" ;
}
close OUT ;

sub read_phy {
    open PHY, "<$_[0]" ;

    my%alignments ; # hash of all alignments in phy
    my@samples ;
    my$header ;
    my$samples ;
    my$length ;
    my%seqs ;   # hash of sequences in each alignment, added to %alignments and cleared when a new alignment is reached
    my$sample_count = 0 ;   # count number of samples for appending sequences
    my$count = 0 ;  # count number of alignments - for unique keys (b/c header might not be specific)

    while (<PHY>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        chomp $_ ;

        # Reach a new header - add last alignment (%seqs) to %alignments and clear %seqs
        if ( $_ =~ m/^ [0-9]/ ) {
#            print "actual_length: ", length($alignments{$count}{"ALIGNMENT"}{"1"}{"SEQ"}), "\n" ;

            $count ++ ;
            $header = $_ ;
            $header =~ m/^ (\d+) (\d+)/ ;
            $samples = $1 ;
            $length = $2 ;

            $alignments{$count}{"HEADER"} = $header ;
            $alignments{$count}{"SAMPLES"} = $samples ;
            $alignments{$count}{"LENGTH"} = $length ;
#            print "Reported_length: ", $length, "\n" ;

            $sample_count = 0 ;

            }

        # This is where the sequences for the alignment start - add sequence names and sequences to %seqs
        elsif ($_ =~ m/^(\w+.+)\s{2,}([ATCGNatcgn-]+)/) {
            my$sample = $1 ;
            my$seq = $_ ;
            $seq =~ s/^\w+.+\s{2,}// ;
            $seq =~ s/\s+//g ;
            $sample_count ++ ;

            $alignments{$count}{"ALIGNMENT"}{$sample_count}{"SAMPLE"} = $sample ;
            $alignments{$count}{"ALIGNMENT"}{$sample_count}{"SEQ"} = $seq ;

            if ($sample ~~ @samples) {
               next ;
            }

            else {
                push @samples, $sample ;
            }

        }

        # Reached the end of one segment of the alignment - clear sample count so next alignment segment can be added in the proper order
        elsif ($_ =~ m/^$/) {
            $sample_count = 0 ;
        }

        # This is within a continued interleaved alignment - append sequences to proper sample in %seqs
        elsif ($_ =~ m/^\s+[ATCGNatcgn-]/) {
            $sample_count ++ ;

            my$seq = $_ ;
            $seq =~ s/\s+//g ;

            $alignments{$count}{"ALIGNMENT"}{$sample_count}{"SEQ"} .= $seq ;
        }

        else {
            print "not_working", "\t", $_, "\n" ;
        }
    }

    close PHY ;
    return \%alignments, \@samples ;
    }