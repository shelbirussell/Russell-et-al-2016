use strict ;
use warnings ;
use Sort::Naturally ;

## usage: perl update_gene_reference.pl variants.vcf genome.fasta genome.gff [/output_dir/]

my$vcf = $ARGV[0] ;
my$fasta = $ARGV[1] ;
my$gff = $ARGV[2] ;

my$out_dir ;
if ($ARGV[3]) {
    $out_dir = $ARGV[3] ;
    }

else { $out_dir = "./" ; }

my$seqs = read_fasta($fasta) ;
my%seqs = %{$seqs} ;

my$genes = read_gff($gff) ;
my%genes = %{$genes} ;

my($sample_names, $variant_sites, $snps, $indels) = parse_vcf($vcf) ;
my%sample_names = %{$sample_names} ;

my%filehandles ;

foreach my$index (keys %sample_names) {
	open my$fh, ">${sample_names{$index}}_nt.fasta" or die "can't open >${sample_names{$index}}_nt.fasta\n" ;
	$filehandles{$sample_names{$index}} = $fh ; 
}

foreach my$scaffold (nsort keys %genes) {
    foreach my$gene (nsort keys %{$genes{$scaffold}}) {

        my$length = abs($genes{$scaffold}{$gene}{"START"} - $genes{$scaffold}{$gene}{"STOP"}) + 1 ;
        my$cut = substr($seqs{$scaffold}, $genes{$scaffold}{$gene}{"START"}-1, $length) ;
	
    	my($snp_report, $snp_seqs) = insert_snps($cut, $sample_names, $snps, $scaffold, $genes{$scaffold}{$gene}{"START"}, $genes{$scaffold}{$gene}{"STOP"}, $length, $gene) ;

	    my($indel_report, $variant_seqs) = insert_indels($snp_seqs, $indels, $scaffold, $genes{$scaffold}{$gene}{"START"}, $genes{$scaffold}{$gene}{"STOP"}, $length, $gene) ;

    	my%snp_report = %{$snp_report} ;
    	my%indel_report = %{$indel_report} ;
    	my%variant_seqs = %{$variant_seqs} ;

	    if (exists $snp_report{"yes"} || exists $indel_report{"yes"}) {
    	    foreach my$sample (keys %variant_seqs) {
        	    my@sample_name = split("_", $sample) ;

            	my$gene_seq = $variant_seqs{$sample} ;
            	my$seq_len = length($gene_seq) ;

	            if ($genes{$scaffold}{$gene}{"STRAND"} eq "-" ) {
    	            $gene_seq = reverse_complement($gene_seq) ;
        	    }
				my$fh = $filehandles{$sample_name[0]} ;
								
		        print $fh ">", $gene, "\t", $scaffold, ":", $genes{$scaffold}{$gene}{"START"}, "-", $genes{$scaffold}{$gene}{"STOP"}, "\tlength:", $length, "\tstrand:", $genes{$scaffold}{$gene}{"STRAND"}, "\tproduct:", $genes{$scaffold}{$gene}{"GENE"}, "\n" ;
    	        print $fh "$_\n" foreach ($gene_seq =~ /.{1,80}/g) ;

    	    }
    	}

	    else {
    	    print "No variants in ", $gene, ", encoding ", $genes{$scaffold}{$gene}{"GENE"}, "\n" ;
        	next ;
    	}
	}
}

foreach my$file (keys %filehandles) {
	close $filehandles{$file} ;
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

sub read_gff {
    open GFF, "<$_[0]" or die "can't open $_[0]\n" ;

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

sub parse_vcf {
    open VCF, "<$_[0]" or die "can't open $_[0]\n" ;
    my%sample_names ;
    my%variant_sites ;
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

                # check for duplicate calls at same position
                $variant_sites{$sample}{$scaff}{$pos}{$call} = $GQ ;

                if ( exists $variant_sites{$index}{$scaff}{$pos}{$call} ) {
                    print "duplicate variant sites in ", $sample_names{$index}, " at pos: ", $scaff, "-", $pos, " with call: ", $call, "\n" ;
                    }
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

    my@samples_2_delete = ("KB30gill_mitochondrial_and_symbiont_genomes","PJRI36testis_mitochondrial_and_symbiont_genomes","PJRI38testis_mitochondrial_and_symbiont_genomes","SRI1gill_mitochondrial_and_symbiont_genomes") ;

    foreach my$delete (@samples_2_delete) {
        delete $snps{$delete} ;
        delete $indels{$delete} ;
        }

    return (\%sample_names, \%variant_sites, \%snps, \%indels) ;
}

sub insert_snps {
    my$gene_seq = $_[0] ;
    my%samples = %{$_[1]} ;
    my%snps = %{$_[2]} ;
    my$scaff = $_[3] ;
    my$gene_start = $_[4] ;
    my$gene_end = $_[5] ;
    my$gene_length = $_[6] ;
    my$gene = $_[7] ;

    my%new_seqs ;
    my%report ;

    my@gene_range = $gene_start..$gene_end ;
    my@seq = split("", $gene_seq) ;

    foreach my$sample (keys %snps) {
        my@var_pos ;
        my$new_seq ;

        foreach my$pos (@gene_range) {
            if ($snps{$sample}{$scaff}{$pos}) {
                push @var_pos, $pos ;
                }
            }

        if (!@var_pos) {
            push @{$report{"no"}}, $sample ;
            $new_seqs{$sample} = $gene_seq ;
            }

        else {
            push @{$report{"yes"}}, $sample ;
            foreach my$site (0..$#seq) {
                my$adjusted_site = $site + $gene_start ;
				
				my@alts ;

                if (grep /$adjusted_site/, @var_pos) {
                	if (! $snps{$sample}{$scaff}{$adjusted_site}{"ALT"}) {
                		print $sample, " ", $scaff, " ", $adjusted_site, " ", $snps{$sample}{$scaff}{$adjusted_site}{"ALT"}, "\n" ;
                	}                	

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
        }

    return (\%report, \%new_seqs) ;

}

sub insert_indels {
    my%snp_seqs = %{$_[0]} ;
    my%indels = %{$_[1]} ;
    my$scaff = $_[2] ;
    my$gene_start = $_[3] ;
    my$gene_end = $_[4] ;
    my$gene_length = $_[5] ;
    my$gene = $_[6] ;

    my%new_seqs ;
    my%report ;
    my@gene_range = $gene_start..$gene_end ;

    foreach my$sample (keys %indels) {
        my@seq = split("", $snp_seqs{$sample}) ;
        my$new_seq ;
        my$current_pos = $gene_start ;
        my$length_change = 0 ;

        my@var_pos ;

        foreach my$pos (@gene_range) {
            if ($indels{$sample}{$scaff}{$pos}) {
                push @var_pos, $pos ;
                }
            }

        if (!@var_pos) {
            push @{$report{"no"}}, $sample ;
            $new_seqs{$sample} = $snp_seqs{$sample} ;
            }

        else {
            push @{$report{"yes"}}, $sample ;
            foreach my$pos (@var_pos) {
                # Append reference sequence up until the indel site

                if ($pos >= $current_pos) {
                    my$length = $pos - $current_pos ;
                    my$adjusted_pos = $current_pos - $gene_start ;
                    my$seq_append = substr($snp_seqs{$sample}, $adjusted_pos, $length) ;
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

                else {
                    if ($pos eq $var_pos[-1]) {
                        $length_change += $current_pos - $gene_range[-1] - 1 ;
                        $current_pos = $gene_range[-1] + 1 ;
                    }

                    else {
                        next ;
                    }
                }
            }

            if ($current_pos <= $gene_range[-1]) {
                my$adjusted_pos = $current_pos - $gene_start ;
                $new_seq .= substr($snp_seqs{$sample}, $adjusted_pos) ;
                }

            if ($current_pos > $gene_range[-1]) {
                $length_change += $current_pos - $gene_range[-1] - 1 ;
                }

            if (length($new_seq) != $length_change+$gene_length) {
                print $sample, " ", $gene, "; new seq length: ", length($new_seq), "; expected length: ", $gene_length, " + ", $length_change, " = ", $length_change+$gene_length, "; new sequence is the wrong length!\n" ;
                }

            else {
#                print $sample, " ", $gene, " completed successfully!\n" ;
                }

            $new_seqs{$sample} = $new_seq ;
            }
        }
    return (\%report, \%new_seqs) ;
}

sub reverse_complement {
    my$dna = shift ;

    my$revcomp = reverse($dna) ;

    $revcomp =~ tr/ACGTacgt/TGCAtgca/ ;

    return $revcomp ;
}
