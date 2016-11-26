use strict ;
use warnings ;
use Data::Dumper ;

my$fasta = $ARGV[0] ;
my$vcf = $ARGV[1] ;
my$out_dir ;

if ($ARGV[2]) {
    $out_dir = $ARGV[2] ;
    }

else {
    $out_dir = "./" ;
    }

my$seq = read_fasta($fasta) ;

my($sample_names, $snps, $indels) = parse_vcf($vcf) ;

my$snp_seqs = insert_snps($seq, $sample_names, $snps) ;

my$variant_seqs = insert_indels($snp_seqs, $indels) ;

my%variant_seqs = %{$variant_seqs} ;

foreach my$sample (keys %variant_seqs) {
    my@sample_name = split("_", $sample) ;
    my$outfile = $out_dir . $sample_name[0] . "_MitoSym" . ".fasta" ;

    open OUT, ">$outfile" or print "Can't open file: ", $sample_name[0], "\n" ;

    foreach my$scaff (keys %{$variant_seqs{$sample}}) {

        my$seq_len = length($variant_seqs{$sample}{$scaff}) ;

        print OUT ">", $scaff, "\t", $sample_name[0], "\t", $seq_len, "\n" ;
        print OUT "$_\n" foreach ($variant_seqs{$sample}{$scaff} =~ /.{1,80}/g) ;

        }

    close OUT ;
    }

sub insert_indels {
    my%snp_seqs = %{$_[0]} ;
    my%indels = %{$_[1]} ;

    my%new_seqs ;

    if (!%indels) {
        print "no indels \n" ;
        %new_seqs = %snp_seqs ;
        }

    else {
        foreach my$sample (keys %indels) {
            foreach my$scaff (keys %{$indels{$sample}}) {

                my$scaff_length = length($snp_seqs{$sample}{$scaff}) ;
                my@seq = split("", $snp_seqs{$sample}{$scaff}) ;
                my$new_seq ;
                my$current_pos = 0 ;
                my$length_change = 0 ;

                foreach my$pos (sort {$a<=>$b} keys %{$indels{$sample}{$scaff}}) {

                    # Append reference sequence up until the indel site
                    if ($pos > $current_pos+1) {
                        my$length = $pos - ($current_pos+1) ;
                        my$seq_append = substr($snp_seqs{$sample}{$scaff}, $current_pos, $length) ;
                        $new_seq .= $seq_append ;
                        $current_pos = $current_pos + $length ;

                        my@alts ;

                        if (index($indels{$sample}{$scaff}{$pos}{"ALT"}, ",") != -1) {
                            @alts = split(",", $indels{$sample}{$scaff}{$pos}{"ALT"}) ;
                            }

                        else {
                            @alts = $indels{$sample}{$scaff}{$pos}{"ALT"} ;
                            }


                        if ($indels{$sample}{$scaff}{$pos}{"CALL"} == 1) {
                            $new_seq .= $alts[0] ;
                            $current_pos += length($indels{$sample}{$scaff}{$pos}{"REF"}) ;

                            my$indel_length = length($alts[0]) - length($indels{$sample}{$scaff}{$pos}{"REF"}) ;

                            $length_change += $indel_length ;
                            }

                        elsif ($indels{$sample}{$scaff}{$pos}{"CALL"} == 2) {
                            $new_seq .= $alts[1] ;
                            $current_pos += length($indels{$sample}{$scaff}{$pos}{"REF"}) ;

                            my$indel_length = length($alts[1]) - length($indels{$sample}{$scaff}{$pos}{"REF"}) ;

                            $length_change += $indel_length ;
                            }

                        else {
                            print $pos, " has more than 3 indel alleles!\n" ;
                        }
                    }
                }

                if ($current_pos < $scaff_length) {
                    $new_seq .= substr($snp_seqs{$sample}{$scaff}, $current_pos) ;
                    }

                if (length($new_seq) != $length_change+$scaff_length) {
                    print $sample, " new seq length: ", length($new_seq), " expected length: ", $scaff_length, " + ", $length_change, " = ", $length_change+$scaff_length, " new sequence is the wrong length!\n" ;
                    }

                $new_seqs{$sample}{$scaff} = $new_seq ;

                }
            }
        }
    return \%new_seqs ;
}

sub insert_snps {
    my%ref_seq = %{$_[0]} ;
    my%samples = %{$_[1]} ;
    my%snps = %{$_[2]} ;

    my%new_seqs ;

    if (!%snps) {
        print "no snps \n";
        foreach my$index (keys %samples) {
            foreach my$header (keys %ref_seq) {
                $new_seqs{$samples{$index}}{$header} = $ref_seq{$header} ;
                }
            }
        }

    else {
        foreach my$sample (keys %snps) {
            foreach my$scaff (keys %{$snps{$sample}}) {
                my@seq = split("", $ref_seq{$scaff}) ;

                my$new_seq ;
                my@var_pos ;

                foreach my$pos (keys %{$snps{$sample}{$scaff}}) {
                    push @var_pos, $pos ;
                    }

                foreach my$site (0..$#seq) {
                    my$adjusted_site = $site + 1 ;

                    if ($adjusted_site ~~ @var_pos) {

                        if (index($snps{$sample}{$scaff}{$adjusted_site}{"ALT"}, ",") != -1) {
                            my@alts = split(",", $snps{$sample}{$scaff}{$adjusted_site}{"ALT"}) ;

                            if ($snps{$sample}{$scaff}{$adjusted_site}{"CALL"} == 1) {
                                $new_seq .= $alts[0] ;
                                }

                            elsif ($snps{$sample}{$scaff}{$adjusted_site}{"CALL"} == 2) {
                                $new_seq .= $alts[1] ;
                                }

                            else {
                                print $sample, " ", $scaff, " ", $adjusted_site, " ", $snps{$sample}{$scaff}{$adjusted_site}{"CALL"}, " ", "call isn't 1 or 2\n" ;
                                }
                            }

                        else {

                            if ($snps{$sample}{$scaff}{$adjusted_site}{"CALL"} == 1) {
                                $new_seq .= $snps{$sample}{$scaff}{$adjusted_site}{"ALT"} ;
                                }

                            else {
                                print $sample, " ", $scaff, " ", $adjusted_site, " ", $snps{$sample}{$scaff}{$adjusted_site}{"CALL"}, " ", "call isn't 1 \n" ;
                                }
                            }
                        }

                    else {
                        $new_seq .= $seq[$site] ;
                        }
                    }

            $new_seqs{$sample}{$scaff} = $new_seq ;

            }
        }
    }

    return \%new_seqs ;

}

sub parse_vcf {
    open VCF, "<$_[0]" or die "cannot open $_[0]\n" ;
    my%sample_names ;
    my%snps ;
    my%indels ;

    while (<VCF>) {
        if ($_ =~ m/^##/) { next ; }

        if($_ =~ m/^#CHROM/) {
            chomp $_ ;
            my@split = split(/\t/, $_) ;

            foreach my$index (9..$#split) {
                $sample_names{$index} = $split[$index] ;
                }
            }

        my @split = split(/\t/, $_) ;

        my$scaff = $split[0] ;
        my$pos = $split[1] ;
        my$ref = $split[3] ;
        my$alt = $split[4] ;
        my%calls ;

        foreach my $index (9..$#split) {
            my$call ;
            my@sample_info = split(":", $split[$index]) ;
            my$GQ = $sample_info[3] ;

            if ( $split[$index] =~ m/^(\d):/ ) {
                $call = $1 ;
                my$sample = $sample_names{$index} ;
                $calls{$sample} = $call ;
            }
        }

        if (length($alt) != length($ref)) {
            # 2 scenarios: indels or multiallele snps

            if ($alt =~ m/([ATCGatcg]{2,})/ || $ref =~ m/([ATCGatcg]{2,})/) {
#                print $pos, " ", $ref, " ", $alt, " indel detected\n" ;
                foreach my$sample (keys %calls) {
                    if ($calls{$sample} == 0) {
                        next ;
                        }
                    else {
                        $indels{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                        $indels{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                        $indels{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
#                        print $sample, " ", $pos, " ", $calls{$sample}, " ", $ref, " ", $alt, "\n" ;
                        }
                    }
                }

            elsif ($alt =~ m/[ATCGatcg]{1},[ATCGatcg]{1}/) {
#                print $pos, " ", $ref, " ", $alt, " multiallele snp detected\n" ;
                foreach my$sample (keys %calls) {
                    if ($calls{$sample} == 0) {
                        next ;
                        }
                    else {
                        $snps{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                        $snps{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                        $snps{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
                        }
                    }
                }

            else {
                print $pos, " ", $ref, " ", $alt, " detected, but uncharacterized\n" ;
                }
            }

        elsif (length($alt) == length($ref)) {
#            print $pos, " ", $ref, " ", $alt, " snp detected\n" ;
            foreach my$sample (keys %calls) {
                if ($calls{$sample} == 0) {
                    next ;
                    }
                else {
                    $snps{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                    $snps{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                    $snps{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
                    }
                }
            }

        else {
            print $pos, " ", $ref, " ", $alt, " not counted! \n" ;
            }
        }
    close VCF ;

    return (\%sample_names, \%snps, \%indels) ;
}

sub read_fasta {
    open FASTA, "<$_[0]" or die "cannot open $_[0]\n" ;

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