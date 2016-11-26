use strict ;
use warnings ;
use Sort::Naturally ;

## usage: perl update_gene_reference.pl 100000 genome.fasta variants.vcf [/output_dir/]

my$window = $ARGV[0] ;
my$fasta = $ARGV[1] ;
my$vcf = $ARGV[2] ;

my$out_dir ;
if ($ARGV[3]) {$out_dir = $ARGV[3] ;}
else { $out_dir = "./" ; }

my$seqs = read_fasta($fasta) ;
my%seqs = %{$seqs} ;

my($sample_names, $snps, $indels) = parse_vcf($vcf) ;
my%sample_names = %{$sample_names} ;

foreach my$scaffold (nsort keys %seqs) {
	my$window_progress = 1 ;
	my$scaff_length = length($seqs{$scaffold}) ;
	
	while ($window_progress < $scaff_length) {

		my$window_end = $window_progress + $window-1 ;

		if ($window_end > $scaff_length) {$window_end = $scaff_length ;}
		
#		print $scaffold,"\t", $window_progress, ":", $window_end, "\t" ;
		
		my$seq = substr($seqs{$scaffold}, $window_progress, $window) ;
#		print length($seq), "\n" ;
		
    	my($snp_seqs) = insert_snps($seq, $sample_names, $snps, $scaffold, $window_progress, $window_end, $window_end-$window_progress+1) ;

	    my($variant_seqs) = insert_indels($snp_seqs, $indels, $scaffold, $window_progress, $window_end, $window_end-$window_progress+1) ;
		my%variant_seqs = %{$variant_seqs} ;
	
    	foreach my$sample (keys %variant_seqs) {
			open FH, ">${sample}_${scaffold}_${window_progress}_${window_end}.fasta" or die "can't open >${scaffold}_${window_progress}_${window_end}.fasta\n" ;

	       	my@sample_name = split("_", $sample) ;

        	my$window_seq = $variant_seqs{$sample} ;
        	my$seq_len = length($window_seq) ;

		    print FH ">", $sample_name[0], "\t", $scaffold, ":", $window_progress, "-", $window_end, "\tlength:", $seq_len, "\n" ;
    	    print FH "$_\n" foreach ($window_seq =~ /.{1,80}/g) ;
	
	    	close FH ;
    	}	
    	
		$window_progress += $window ;
	}
}

sub read_fasta {
    open FASTA, "<$_[0]" or die "can't open $_[0]\n" ;

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

sub parse_vcf {
    open VCF, "<$_[0]" or die "can't open $_[0]\n" ;
    my%sample_names ;
    my%snps ;
    my%indels ;

    while (<VCF>) {
        if ($_ =~ m/^##/) { next ; }

        if($_ =~ m/^#CHROM/) {
            chomp $_ ;
            my@split = split(/\t/, $_) ;

            foreach my$index (9..$#split) {
            	$split[$index] =~ m/(\w+\d+gill)/ ;
            	my$sample = $1 ;
                $sample_names{$index} = $sample ;
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
                    if ($calls{$sample} == 0) {next ;}

                    else {
                        $indels{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                        $indels{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                        $indels{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
                        }
                    }
                }

            elsif ($alt =~ m/[ATCGatcg]{1},[ATCGatcg]{1}/) {
#                print $pos, " ", $ref, " ", $alt, " multiallele snp detected\n" ;
                foreach my$sample (keys %calls) {
                    if ($calls{$sample} == 0) {next ;}

                    else {
                        $snps{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                        $snps{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                        $snps{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
                        }
                    }
                }

            else {print $pos, " ", $ref, " ", $alt, " detected, but uncharacterized\n" ;}
            }

        elsif (length($alt) == length($ref)) {
#            print $pos, " ", $ref, " ", $alt, " snp detected\n" ;
            foreach my$sample (keys %calls) {
                if ($calls{$sample} == 0) {next ;}

                else {
                    $snps{$sample}{$scaff}{$pos}{"REF"} = $ref ;
                    $snps{$sample}{$scaff}{$pos}{"ALT"} = $alt ;
                    $snps{$sample}{$scaff}{$pos}{"CALL"} = $calls{$sample} ;
                    }
                }
            }

        else {print $pos, " ", $ref, " ", $alt, " not counted! \n" ;}
    }
    close VCF ;

    my@samples_2_delete = ("KB30gill_mitochondrial_and_symbiont_genomes","PJRI36testis_mitochondrial_and_symbiont_genomes","PJRI38testis_mitochondrial_and_symbiont_genomes","SRI1gill_mitochondrial_and_symbiont_genomes") ;

    foreach my$delete (@samples_2_delete) {
        delete $snps{$delete} ;
        delete $indels{$delete} ;
    }

    return (\%sample_names, \%snps, \%indels) ;
}

sub insert_snps {
    my$sequence = $_[0] ;
    my%samples = %{$_[1]} ;
    my%snps = %{$_[2]} ;
    my$scaff = $_[3] ;
    my$start = $_[4] ;
    my$end = $_[5] ;
    my$length = $_[6] ;

    my%new_seqs ;

    my@range = $start..$end ;
    my@seq = split("", $sequence) ;

    foreach my$sample (keys %snps) {
        my$new_seq ;
        
        foreach my$site (0..$#seq) {
            my$adjusted_site = $site + $start ;
				
			my@alts ;
				
            if ($snps{$sample}{$scaff}{$adjusted_site}) {

                if (index($snps{$sample}{$scaff}{$adjusted_site}{"ALT"}, ",") != -1) {		#searches for comma in alt record, indicating multiple alternate alleles
					@alts = split(",", $snps{$sample}{$scaff}{$adjusted_site}{"ALT"}) ;
				}

				else {
					@alts = $snps{$sample}{$scaff}{$adjusted_site}{"ALT"} ;
				}	

                if ($snps{$sample}{$scaff}{$adjusted_site}{"CALL"} == 1) {
                    $new_seq .= $alts[0] ;
                }

                elsif ($snps{$sample}{$scaff}{$adjusted_site}{"CALL"} == 2) {
                    $new_seq .= $alts[1] ;
                }

            	else {
                    print $sample, " ", $scaff, " ", $adjusted_site, " ", $snps{$sample}{$scaff}{$adjusted_site}{"CALL"}, " ", "call isn't 1 \n" ;
                }
            }

            else {
                $new_seq .= $seq[$site] ;
            }
        }

        $new_seqs{$sample} = $new_seq ;

    }

    return (\%new_seqs) ;

}

sub insert_indels {
    my%snp_seqs = %{$_[0]} ;
    my%indels = %{$_[1]} ;
    my$scaff = $_[2] ;
    my$start = $_[3] ;
    my$end = $_[4] ;
    my$length = $_[5] ;

    my%new_seqs ;
    my@range = $start..$end ;

    foreach my$sample (keys %indels) {
        my@seq = split("", $snp_seqs{$sample}) ;
        my$new_seq ;
        my$current_pos = $start ;
        my$length_change = 0 ;

        my@var_pos ;

        foreach my$pos (@range) {if ($indels{$sample}{$scaff}{$pos}) {push @var_pos, $pos ;}}

        if (!@var_pos) {$new_seqs{$sample} = $snp_seqs{$sample} ;}

        else {
            foreach my$pos (@var_pos) {

                # Append reference sequence up until the indel site
                if ($pos >= $current_pos) {
                    my$length = $pos - $current_pos ;
                    my$adjusted_pos = $current_pos - $start ;
                    my$seq_append = substr($snp_seqs{$sample}, $adjusted_pos, $length) ;
                    $new_seq .= $seq_append ;
                    $current_pos = $current_pos + $length ;
                    
                    my@alts ;

					# Make list of alternate alleles
                    if (index($indels{$sample}{$scaff}{$pos}{"ALT"}, ",") != -1) {@alts = split(",", $indels{$sample}{$scaff}{$pos}{"ALT"}) ;}

                    else {@alts = $indels{$sample}{$scaff}{$pos}{"ALT"} ;}

					# Assign allele based on call
                    if ($indels{$sample}{$scaff}{$pos}{"CALL"} == 1) {
                        $new_seq .= $alts[0] ;
                        $current_pos += length($indels{$sample}{$scaff}{$pos}{"REF"}) ;		# Advance position past reference length
                        my$indel_length = length($alts[0]) - length($indels{$sample}{$scaff}{$pos}{"REF"}) ;
                        $length_change += $indel_length ;
                    }

                    elsif ($indels{$sample}{$scaff}{$pos}{"CALL"} == 2) {
                        $new_seq .= $alts[1] ;
                        $current_pos += length($indels{$sample}{$scaff}{$pos}{"REF"}) ;		# Advance position past reference length
                        my$indel_length = length($alts[1]) - length($indels{$sample}{$scaff}{$pos}{"REF"}) ;
                        $length_change += $indel_length ;
                    }

                    else {print $pos, " has more than 3 indel alleles!\n" ;}
                }
            }

			# Append remaining sequence after all variant sites are accounted for
            if ($current_pos <= $range[-1]) {
                my$adjusted_pos = $current_pos - $start ;
                $new_seq .= substr($snp_seqs{$sample}, $adjusted_pos) ;
            }

            if (length($new_seq) != $length_change+$length) {
                print $sample, " ", $scaff, ":", $start, "-", $end, "; new seq length: ", length($new_seq), "; expected length: ", $length, " + ", $length_change, " = ", $length_change+$length, "; new sequence is the wrong length!\n" ;
            }

#            else {print $sample, " ", $gene, " completed successfully!\n" ;}

            $new_seqs{$sample} = $new_seq ;
            }
        }
    return \%new_seqs ;
}