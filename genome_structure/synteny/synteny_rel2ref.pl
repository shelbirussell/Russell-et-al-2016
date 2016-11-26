use strict ;
use warnings ;
use Sort::Naturally ;
use File::Basename ;

### Usage: perl synteny.pl ~/Reference/Solemya_velum_MitoSym_top4.gff ./gff/DML4gill_sym_filtered_draft_genome.gff DML4_multireciprocal.out 

my$ref_gff = $ARGV[0] ;
my$sample_gff = $ARGV[1] ;
my$reciprocal_hits = $ARGV[2] ;

my$out_all = basename($reciprocal_hits) ;
my$out_subset = basename($reciprocal_hits) ;
$out_all =~ s/.out/_all_orthologs.out/ ;
$out_subset =~ s/.out/_rearranged_orthologs.out/ ;

my$ref_annotation = read_Svsym_gff($ref_gff) ;
my$sample_annotation = read_gff($sample_gff) ;

my%ref_annotation = %{$ref_annotation} ;
my%sample_annotation = %{$sample_annotation} ;

my%reciprocals ;
open IN, "<$reciprocal_hits" or die "cannot open $reciprocal_hits \n" ;
while (<IN>) {
	chomp ;
	my@split = split(/\t/, $_) ;
	$reciprocals{$split[1]} = $split[0] ;
}

close IN ;

open OUT_ALL, ">$out_all" or die "cannot open $out_all\n" ;
open OUT_SUBSET, ">$out_subset" or die "cannot open $out_subset\n" ;

my%map ;

##read along sample gff and make hash of orthologous ref coordinate for every site
foreach my$scaffold (keys %sample_annotation) {
	foreach my$count (keys %{$sample_annotation{$scaffold}}) {
				
		if ($sample_annotation{$scaffold}{$count}{"ID"} && $reciprocals{$sample_annotation{$scaffold}{$count}{"ID"}}) {
			my$ref_ortholog = $reciprocals{$sample_annotation{$scaffold}{$count}{"ID"}} ; 
			my$ref_count ;
			my$ref_scaffold ;
			foreach my$scaffold (keys %ref_annotation) {	
				if ($ref_annotation{$scaffold}{$ref_ortholog}) {
					$ref_count = $ref_annotation{$scaffold}{$ref_ortholog}{"COUNT"} ;
					$ref_scaffold = $scaffold ;
				}
			}
		
			$map{$scaffold}{$count} =  $ref_scaffold . "\t" . $ref_count ;
		}
	}
}

my%reordered ;

my$prev_ref_scaff = "" ;
my$prev_ref_count = "" ;
my$prev_sample_scaff = "" ;
my$prev_sample_count = "" ;

foreach my$scaffold (nsort keys %map) {
	foreach my$count (nsort keys %{$map{$scaffold}}) {
		
		my@split = split(/\t/, $map{$scaffold}{$count}) ;
		my$ref_scaff = $split[0] ;
		my$ref_count = $split[1] ;

		print OUT_ALL $scaffold, "\t", $count, "\t", $map{$scaffold}{$count}, "\n" ;
		
		if ($prev_ref_scaff eq $ref_scaff) {
			if ($prev_sample_scaff eq $scaffold) {
				## Skip sites in continuous order
				if ($count == $prev_sample_count + 1 && $ref_count == $prev_ref_count + 1 || $count == $prev_sample_count + 1 && $ref_count == $prev_ref_count - 1) {next ;}
				
				## Skip sites that jump the same number of positions in ref and sample
				elsif (abs($count-$prev_sample_count) == abs($ref_count-$prev_ref_count)) {next ;}

				## Write sites that are in the wrong order to OUT_SUBSET
				else {
					my$offset = abs($count-$prev_sample_count)-abs($ref_count-$prev_ref_count) ;
					$reordered{"mismatched_order"}{$offset}{"$scaffold\t$count"}{"GENE"} = $sample_annotation{$scaffold}{$count}{"ID"} ;
					$reordered{"mismatched_order"}{$offset}{"$scaffold\t$count"}{"REF"} = $map{$scaffold}{$count} ;					
				}
			}
		
			else {
				$reordered{"mismatched_sample_scaff"}{"$scaffold\t$count"}{"GENE"} = $sample_annotation{$scaffold}{$count}{"ID"} ;
				$reordered{"mismatched_sample_scaff"}{"$scaffold\t$count"}{"REF"} = $map{$scaffold}{$count} ;					
			}	
		}
		
		else {
			if ($prev_sample_scaff eq $scaffold) {
				$reordered{"mismatched_reference_scaff"}{"$scaffold\t$count"}{"GENE"} = $sample_annotation{$scaffold}{$count}{"ID"} ;
				$reordered{"mismatched_reference_scaff"}{"$scaffold\t$count"}{"REF"} = $map{$scaffold}{$count} ;					
			}
		}	

		$prev_ref_count = $ref_count ;
		$prev_ref_scaff = $ref_scaff ;
		$prev_sample_count = $count ;
		$prev_sample_scaff = $scaffold ;

	}
}

foreach my$category (nsort keys %reordered) {
	print OUT_SUBSET $category, "\n" ;
		
	if ($category =~ m/mismatched_\w+_scaff/) {
		foreach my$site (nsort keys %{$reordered{$category}}) {
			my@ref_coords = split(/\t/, $reordered{$category}{$site}{"REF"}) ;
			
			print OUT_SUBSET "\tSample_scaffold_count:", $site, "\tSample_gene_ID:", $reordered{$category}{$site}{"GENE"}, "\tReference_scaffold_count:", $reordered{$category}{$site}{"REF"}, "\tReference_gene_ID:", $reciprocals{$reordered{$category}{$site}{"GENE"}},"\tReference_gene_product:", $ref_annotation{$ref_coords[0]}{$reciprocals{$reordered{$category}{$site}{"GENE"}}}{"PRODUCT"}, "\n" ;
		}
	}	

	elsif ($category eq "mismatched_order") {
		foreach my$offset (nsort keys %{$reordered{$category}}) {
			print OUT_SUBSET "\t", $offset, "\n" ;
			
			foreach my$site (nsort keys %{$reordered{"mismatched_order"}{$offset}}) {
				my@ref_coords = split(/\t/, $reordered{$category}{$offset}{$site}{"REF"}) ;
				print OUT_SUBSET "\t\tSample_scaffold_count:", $site, "\tSample_gene_ID:", $reordered{$category}{$offset}{$site}{"GENE"}, "\tReference_scaffold_count:", $reordered{$category}{$offset}{$site}{"REF"}, "\tReference_gene_ID:", $reciprocals{$reordered{$category}{$offset}{$site}{"GENE"}},"\tReference_gene_product:", $ref_annotation{$ref_coords[0]}{$reciprocals{$reordered{$category}{$offset}{$site}{"GENE"}}}{"PRODUCT"}, "\n" ;
			}
		}
	}
	
	else {print "Unrecognized:", $category, "\n" ;}
}			

close OUT_ALL ;
close OUT_SUBSET ;

sub read_gff {
    my$cds = $_[0] ;
    open CDS, "<$cds" or die "can't open gff file: $_[0]" ;

    my%gff ;

	my%counts ;
		
    while (<CDS>) {
        if ($_ =~ m/^#/) {
            next ;
        }

        chomp ;

        my@split = split(/\t/, $_) ;
        my$scaffold = $split[0] ;
        my$type = $split[2] ;
        my$start = $split[3] ;
        my$stop = $split[4] ;
        my$strand = $split[6] ;
        my@info = split(";", $split[8]) ;
        my$id = $info[0] ;
        my$product = $info[1] ;
        my$hit_info = $info[1] ;
        my$ref = $info[2] ;

		if ($counts{$scaffold}) {$counts{$scaffold} ++ ;}
		
		else {$counts{$scaffold} = 1 ;}

        $product =~ s/Name=// ;
        $id =~ s/ID=// ;

        my$taxon ;
        if ($hit_info =~ m/Tax=([A-Za-z]+_)/) {
            $taxon = $1 ;
        }

        elsif ($hit_info =~ m/OS=([A-Za-z]+_)/) {
            $taxon = $1 ;
        }

        elsif ($hit_info =~ m/\[([A-Za-z]+_)/) {
            $taxon = $1 ;
        }

        else {
            $taxon = $hit_info ;
        }

        $gff{$scaffold}{$counts{$scaffold}}{"LINE"} = $_ ;
        $gff{$scaffold}{$counts{$scaffold}}{"ID"} = $id ;
        $gff{$scaffold}{$counts{$scaffold}}{"START"} = $start ;
        $gff{$scaffold}{$counts{$scaffold}}{"STOP"} = $stop ;
        $gff{$scaffold}{$counts{$scaffold}}{"TYPE"} = $type ;
        $gff{$scaffold}{$counts{$scaffold}}{"PRODUCT"} = $product ;
        $gff{$scaffold}{$counts{$scaffold}}{"TAXON"} = $taxon ;
		$gff{$scaffold}{$counts{$scaffold}}{"STRAND"} = $strand ;
		
        if ($ref eq "Dbxref=na" && $type eq "CDS") {
            $gff{$scaffold}{$start}{"SOURCE"} = "prodigal" ;
        }

        elsif ($ref eq "Dbxref=na" && $type eq "rRNA") {
            $gff{$scaffold}{$start}{"SOURCE"} = "RNAmmer" ;
        }

        elsif ($ref eq "Dbxref=na" && $type eq "tRNA") {
            $gff{$scaffold}{$start}{"SOURCE"} = "tRNAscan" ;
        }

        else {
            $ref =~ s/Dbxref=// ;
            $gff{$scaffold}{$start}{"SOURCE"} = $ref ;
        }

    }

    close CDS ;

    return \%gff ;

}

sub read_Svsym_gff {
    my$gff = $_[0] ;
	open GFF, "<$gff" ;

    my %gff ;
	my%counts ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;
        my$scaffold = $split[0] ;

        if ( $scaffold eq "Sv_mito_chromosome" ) { next ; }

		if ($counts{$scaffold}) {$counts{$scaffold} ++ ;}
		
		else {$counts{$scaffold} = 1 ;}

#        if ( $split[2] ne "CDS" ) { next ; }

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

        $gff{$scaffold}{$id}{"COUNT"} = $counts{$scaffold} ;
        $gff{$scaffold}{$id}{"START"} = $split[3] ;
        $gff{$scaffold}{$id}{"STOP"} = $split[4] ;
        $gff{$scaffold}{$id}{"STRAND"} = $split[6] ;
        $gff{$scaffold}{$id}{"PRODUCT"} = $gene ;

       }
    close GFF ;
    return \%gff ;
    }
