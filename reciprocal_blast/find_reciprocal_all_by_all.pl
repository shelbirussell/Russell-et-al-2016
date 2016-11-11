use strict ;
use warnings ;
use File::Basename ;

#Usage: perl find_reciprocal.pl base_spp(E.g.Sv) Svsym.gff Sp1_vs_Sp2.txt Sp2_vs_Sp1.txt Sp1_vs_Sp3.txt Sp3_vs_Sp1.txt Sp2_vs_Sp3.txt Sp3_vs_Sp2.txt ...
#(the species names must appear X_vs_Y to be recognized)
## E.g. perl find_reciprocal_all_by_all.pl Sv ~/Reference/Solemya_velum_MitoSym_top4.gff Se_vs_Sp_nt_filtered_hits.txt Sp_vs_Se_nt_filtered_hits.txt Sveles_vs_Se_nt_filtered_hits.txt Sv_vs_Se_nt_filtered_hits.txt Se_vs_Sv_nt_filtered_hits.txt Sp_vs_Sv_nt_filtered_hits.txt Sveles_vs_Sp_nt_filtered_hits.txt Sv_vs_Sp_nt_filtered_hits.txt Se_vs_Sveles_nt_filtered_hits.txt Sp_vs_Sveles_nt_filtered_hits.txt Sveles_vs_Sv_nt_filtered_hits.txt Sv_vs_Sveles_nt_filtered_hits.txt

my$base_seq = $ARGV[0] ;
my$base_gene_set = $ARGV[1] ;
my%blasts ;

my$base_genes = read_gff($base_gene_set) ;
my%base_genes = %{$base_genes} ;

foreach my$index (2..$#ARGV) {
    my$file = $ARGV[$index] ;
    my$reciprocal ;

    if ($file =~ m/^(\w+_vs_\w+)_nt_/) {
        $reciprocal = $1 ;
    }

    $blasts{$reciprocal} = parse_blast($file) ;
}

## Find reciprocal best matches for each reciprocal blast
## record in base_matches if $base_seq is one of pair, otherwise, save to $other_matches

my%base_matches ;
my%other_matches ;

foreach my$reciprocal (keys %blasts) {
    $reciprocal =~ m/(.+)_vs_(.+)/ ;
    my$sp1 = $1 ;
    my$sp2 = $2 ;
    my$opposite = "${sp2}_vs_${sp1}" ;

#	print $sp1, "\t", $sp2, "\n" ;
#    print $reciprocal, "\t", $opposite, "\n" ;

    my%sp1vssp2 = %{$blasts{$reciprocal}} ;   ## sp1vsp2
    my%sp2vssp1 = %{$blasts{$opposite}} ;   ## sp2vsp1

    @{$other_matches{$reciprocal}} = () ;

    foreach my$query (keys %{$blasts{$reciprocal}}) {

        my$subject1 = $sp1vssp2{$query} ; ## subject 1 = query for sp 2
        my$subject2 = $sp2vssp1{$subject1} ; ## subject 2 = query for sp 1

#		print $subject1, "\t", $subject2, "\n" ;
		
        ## Skip genes that didn't have a BLAST hit that passed the filters
        if (! $subject2) {
            next ;
        }

        if ($query eq $subject2) {
#            print $sp1, "\t", $subject2, "\t", $sp2, "\t", $subject1, "\t", $query, "\n" ;
            if ($sp1 eq $base_seq) {
                $base_matches{$subject2}{$sp2} = $subject1 ;
            }

            elsif ($sp2 eq $base_seq) {
                $base_matches{$subject1}{$sp1} = $subject2 ;
            }

            else {
                if ("$subject2\t$subject1" ~~ @{$other_matches{$reciprocal}}) {
                    next ;
                }

                push @{$other_matches{$reciprocal}}, "$subject2\t$subject1" ; ## this saves genes in same spp order as in $reciprocal name
            }
        }
    }
}

## Fill out base_matches hash for all samples
## This compares the reciprocal matches for non-base_name seqs inferred from reciprocal matches to the base_name seq
## to those from the actual reciprocal blasts for those seqs.  If they are not equal, the gene is deleted from the hash

my@spp = () ; # make list of species, in sorted order for output header

my%recips ;
foreach my$gene (keys %base_matches) {
    ##Make hash of inferred reciprocal matches
    foreach my$sp1 (sort keys %{$base_matches{$gene}}) {

        if (! grep(/$sp1/, @spp)) {
            push @spp, $sp1 ;
        }

        foreach my$sp2 (sort keys %{$base_matches{$gene}}) {
            if ($sp1 eq $sp2) {
                next ;
            }

            if ($sp1 eq $base_seq || $sp2 eq $base_seq) {
                next ;
            }

            $recips{"${sp1}_vs_${sp2}"} = "$base_matches{$gene}{$sp1}\t$base_matches{$gene}{$sp2}" ;
            $recips{"${sp2}_vs_${sp1}"} = "$base_matches{$gene}{$sp2}\t$base_matches{$gene}{$sp1}" ;
        }
    }
    ##See if blast reciprocal matches predicted reciprocal, if not, delete gene from list
    foreach my$pair (keys %recips) {
        if (exists $other_matches{$pair}) {
            if ($recips{$pair} ~~ $other_matches{$pair}) {
                next ;
            }

            else {
                my@split = split(/\t/, $recips{$pair}) ;
                my@blast_match1 = grep {/$split[0] /} @{$other_matches{$pair}} ;
                my@blast_match2 = grep {/$split[1] /} @{$other_matches{$pair}} ;

                print "Removed gene: ", "Sv_${gene}", "\tpair: ", $recips{$pair}, " doesn't equal reciprocal blast results: ", join(", ", @blast_match1), " or ", join(", ", @blast_match2), "\n" ;
                delete $base_matches{$gene} ;
                last ;
            }
        }
    }
}

my$output = join("_", @spp) . "_multireciprocal.out" ;

open OUT, ">$output" or die "can't open $output" ;

print OUT "gene\t", join("\t", @spp), "\n" ;

foreach my$gene (sort keys %base_matches) {
    print OUT $gene, "\t" ;

    foreach my$sp (@spp) {
        if ($base_matches{$gene}{$sp}) {
            print OUT $base_matches{$gene}{$sp}, "\t", ;
        }

        else {
            print OUT ".\t" ;
        }
    }

    print OUT "\n" ;
}


close OUT ;


sub read_gff {
    open GFF, "<$_[0]" ;

    my %data ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;

        my $id = "" ;
        if ( $split[8] =~ m/ID=(.+);/ ) {
            $id = $1 ;
        }

        my $gene = "" ;
        if ( $split[8] =~ m/Name=(.+);/ ) {
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

sub parse_blast {
    my%hits ;
    open REPORT, "<$_[0]" ;

    while (<REPORT>) {
        if ($_ =~ m/^#/) {
             next ;
        }

        my@info = split(/\t/, $_) ;
        my$query = $info[0] ;
        my$subject ;
		
		if ($info[1] =~ m/#/) {
			my@split = split("#", $info[1]) ;
			$subject = $split[0] ;
		}
		
		else {
			$subject = $info[1] ;
		}
		
        my$product = $info[2] ;
        $product =~ s/ /_/g ;
        my$pident = $info[4] ;
        my$evalue = $info[12] ;

        $hits{$query} = $subject ;
    }

    close REPORT;

    return \%hits ;
}