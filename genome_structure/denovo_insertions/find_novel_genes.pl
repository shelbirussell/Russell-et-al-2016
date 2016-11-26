use strict ;
use warnings ;
use Sort::Naturally ;
use File::Basename ;

# Usage: perl find_novel_genes.pl /Users/shelbirussell/Dropbox/lab/projects/svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/Svelum_syms_reordered/reordered_scaffolds/PJRI36gill_sym_genome_reordered.gff /Users/shelbirussell/Dropbox/lab/projects/svelum_projects/1-PopGen_transmission_project/synteny/parse_gene_order/results/blast_results/PJRI36_vs_Sv_nt_hits.txt

# Finds genes in de novo assemblies with no match in Sv symbiont reference
# 1. Read in sample gff and blast result against the Sv reference genome
# 2. Report genes in gff that had no hits in the Sv reference

my$sample_gff = $ARGV[0] ;
my$sample_blast = $ARGV[1] ;

my$sample_genes = read_gff($sample_gff) ;
my%sample_genes = %{$sample_genes} ;

my$blast = parse_blast($sample_blast) ;
my%blast = %{$blast} ;

my$output = basename($sample_blast) ;
$output =~ s/hits.txt/novel_genes.out/ ;
open OUT, ">$output" or die "cannot open $output\n" ; 

# Loop through annotation for genes that aren't in the blast results
foreach my$scaffold (nsort keys %sample_genes) {
	foreach my$id (nsort keys %{$sample_genes{$scaffold}}) {
		if (exists $blast{$id}) { next ;}
		
		else {
			# Exclude genes with no hits that are annotated w/ Sv sequences (evalue was below 1e-6 cutoff)
			if ($sample_genes{$scaffold}{$id}{"PRODUCT"} =~ m/Solemya_velum_gill_symbiont/) {next ;}
			
			print OUT $scaffold, "\t", $id, "\t", $sample_genes{$scaffold}{$id}{"PRODUCT"}, "\n" ;
		}
	}
}

sub read_gff {
    my$cds = $_[0] ;
    open CDS, "<$cds" or die "can't open gff file: $_[0]" ;

    my%gff ;

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

        $gff{$scaffold}{$id}{"LINE"} = $_ ;
        $gff{$scaffold}{$id}{"START"} = $start ;
        $gff{$scaffold}{$id}{"STOP"} = $stop ;
        $gff{$scaffold}{$id}{"TYPE"} = $type ;
        $gff{$scaffold}{$id}{"PRODUCT"} = $product ;
        $gff{$scaffold}{$id}{"TAXON"} = $taxon ;
		$gff{$scaffold}{$id}{"STRAND"} = $strand ;
		
        if ($ref eq "Dbxref=na" && $type eq "CDS") {
            $gff{$scaffold}{$id}{"SOURCE"} = "prodigal" ;
        }

        elsif ($ref eq "Dbxref=na" && $type eq "rRNA") {
            $gff{$scaffold}{$id}{"SOURCE"} = "RNAmmer" ;
        }

        elsif ($ref eq "Dbxref=na" && $type eq "tRNA") {
            $gff{$scaffold}{$id}{"SOURCE"} = "tRNAscan" ;
        }

        else {
            $ref =~ s/Dbxref=// ;
            $gff{$scaffold}{$id}{"SOURCE"} = $ref ;
        }

    }

    close CDS ;

    return \%gff ;

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