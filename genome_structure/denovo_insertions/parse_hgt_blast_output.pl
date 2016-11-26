use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# This script reads in blast hits and stores them by sample 
# Then, for each sample, the sample gff is read in, and filter out reports for scaffolds without > $min_Sv_matches Svelum symbiont sequence on it
# Filter out genes with hits < 60% similarity and < 60% coverage of the genome sequence <-- These are not strongly supportive of HGT
# Report (1) all results, (2) a summary of all results, (3) results with > 98% similarity, and (4) a summary of results with > 98% similarity

## Usage: perl parse_hgt_blast_output.pl results-HGT_genes/novel_vs_wgs_nt_filtered_hits.txt /Users/shelbirussell/Dropbox/Lab/projects/Svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/Svelum_syms results-HGT_genes/subject_ids.txt

my$blast_results_file = $ARGV[0] ;
my$gff_dir = $ARGV[1] ;
my$subject_ids = $ARGV[2] ;

my$min_Sv_matches = 4 ;

my$blast_results = parse_blast_results($blast_results_file) ;
my%blast_results = %{$blast_results} ;

opendir(DIR, $gff_dir) ;
my@files = readdir(DIR) ;
my@gff = () ;
foreach my$file (@files) {
	if ($file =~ m/.gff$/) {
		push @gff, $file ;
	}
}

my$ids = parse_tabs($subject_ids) ;
my%ids = %{$ids} ;

## Only keep novel gene reports for genes on scaffolds with strong support for being of symbiont origin
## --> This step checks for sequences matching the Sv symbiont
foreach my$file (@gff) {
	my$sample = $file ;
	$sample =~ s/^(\w+\d+).+\.gff/$1/ ;

	my($anno, $Sv_matches) = read_gff($gff_dir, $file) ;
	my%Sv_matches = %{$Sv_matches} ;
	my%anno = %{$anno} ;
	
	foreach my$scaffold (nsort keys %{$blast_results{$sample}}) {
		if (exists $Sv_matches{$scaffold} && $Sv_matches{$scaffold} > $min_Sv_matches) {
			foreach my$gene (nsort keys %{$blast_results{$sample}{$scaffold}}) {
				$blast_results{$sample}{$scaffold}{$gene}{"PRODUCT"} = $anno{$scaffold}{$gene}{"GENE"} ;
			}
		}
		else {delete $blast_results{$sample}{$scaffold} ;}
	}
}
		
my$out1 = basename($blast_results_file) ;
$out1 =~ s/filtered_hits.txt/output.txt/ ;
my$out2 = basename($blast_results_file) ;
$out2 =~ s/filtered_hits.txt/summary_output.txt/ ;
my$out3 = basename($blast_results_file) ;
$out3 =~ s/filtered_hits.txt/gt98id_output.txt/ ;
my$out4 = basename($blast_results_file) ;
$out4 =~ s/filtered_hits.txt/gt98id_summary_output.txt/ ;

open OUT1, ">$out1" or die "cannot open $out1\n" ;
open OUT2, ">$out2" or die "cannot open $out2\n" ;
open OUT3, ">$out3" or die "cannot open $out3\n" ;
open OUT4, ">$out4" or die "cannot open $out4\n" ;

## Output table of well-supported novel genes for each sample
foreach my$sample (nsort keys %blast_results) {
#	print $sample, "\n" ;
	
	my%report ;
	my%report_gt98 ;
	my$novel = 0 ;
	my$min_id = 100 ;
	my$max_id = 0 ;
	
	foreach my$scaffold (nsort keys %{$blast_results{$sample}}) {
		foreach my$gene (nsort keys %{$blast_results{$sample}{$scaffold}}) {
			my$product = lc($blast_results{$sample}{$scaffold}{$gene}{"PRODUCT"}) ;
			my$hit = $ids{$blast_results{$sample}{$scaffold}{$gene}{"MATCH"}}{"NAME"} ;
			my$rRNAsmsu = $ids{$blast_results{$sample}{$scaffold}{$gene}{"MATCH"}}{"16SID"} ;
			$rRNAsmsu =~ s/[~<=]// ;
			my$id = $blast_results{$sample}{$scaffold}{$gene}{"ID"} ;
			my$cov = $blast_results{$sample}{$scaffold}{$gene}{"COV"} ;

			# Remove low % similarity hits with < 60% of gene covered by hit 
			if ($id < 60 && $cov > 60) { next; }
			
			print OUT1 $sample, "\t", $scaffold, "\t", $gene, "\t", $hit, "\t", $id, "\t", $cov, "\t", $rRNAsmsu, "\t", $product, "\n" ;
			
			# Report recent HGT events
			if ($id > 98 && $rRNAsmsu < 90) {
				print OUT3 $sample, "\t", $scaffold, "\t", $gene, "\t", $hit, "\t", $id, "\t", $cov, "\t", $rRNAsmsu, "\t", $product, "\n" ;
				
				if ($report_gt98{"PRODUCT"}{$product}) {$report_gt98{"PRODUCT"}{$product} ++ ;}	
				else {$report_gt98{"PRODUCT"}{$product} = 1 ;}
				
				if ($report_gt98{"HITS"}{$hit}) {$report_gt98{"HITS"}{$hit} ++ ;}
				else {$report_gt98{"HITS"}{$hit} = 1 ;}
			}
					
			if ($report{"PRODUCT"}{$product}) {$report{"PRODUCT"}{$product} ++ ;}	
			else {$report{"PRODUCT"}{$product} = 1 ;}
			if ($report{"HITS"}{$hit}) {$report{"HITS"}{$hit} ++ ;}
			else {$report{"HITS"}{$hit} = 1 ;}

			if ($id > $max_id) {$max_id = $id ;}
			if ($id < $min_id) {$min_id = $id ;}
		
			$novel ++ ;
		}
	}

	my@functions = () ;
	foreach my$product (nsort keys %{$report{"PRODUCT"}}) {
		push @functions, $product . "(" . $report{"PRODUCT"}{$product} . ")" ;
	}
	
	my@hits = () ;
	foreach my$hit (nsort keys %{$report{"HITS"}}) {
		push @hits, $hit . "(" . $report{"HITS"}{$hit} . ")" ;
	}	
	
	print OUT2 $sample, "\t", $novel, "\t", $min_id, "-", $max_id, "\t", join(",", @functions), "\t", join(",", @hits), "\n" ;

	my@functions_gt98 = () ;
	foreach my$product (nsort keys %{$report_gt98{"PRODUCT"}}) {
		push @functions_gt98, $product . "(" . $report_gt98{"PRODUCT"}{$product} . ")" ;
	}
	
	my@hits_gt98 = () ;
	foreach my$hit (nsort keys %{$report_gt98{"HITS"}}) {
		push @hits_gt98, $hit . "(" . $report_gt98{"HITS"}{$hit} . ")" ;
	}	

	if (scalar(@hits_gt98) > 0) {
		print OUT4 $sample, "\t", scalar(@hits_gt98), "\t", join(",", @functions_gt98), "\t", join(",", @hits_gt98), "\n" ;
	}
}		

close OUT1 ;
close OUT2 ;
sub parse_blast_results {
	my$file = $_[0] ;
	
	my%results ;
	
	open IN, "<$file" or die "cannot open $file\n" ;
	while (<IN>) {
		chomp ;
		
		my@split = split(/\t/, $_) ;
		
		my$id = $split[0] ;
		my$match = $split[1] ;
		my$cov = $split[2] ;
		my$ident = $split[3] ;
		my$length = $split[4] ;
		my$mismatch = $split[5] ;
		my$gaps = $split[6] ;
		my$qstart = $split[7] ;
		my$qend = $split[8] ;
		my$sstart = $split[9] ;
		my$send = $split[10] ;
		my$evalue = $split[11] ;
		my$bitscore = $split[12] ;

		my$sample = $id ;
		$sample =~ s/(\w+\d+)(_sym_scf\d+)-\d+/$1/ ;
		my$scaff = $1 . $2 ;
		
		$results{$sample}{$scaff}{$id}{"MATCH"} = $match ;
		$results{$sample}{$scaff}{$id}{"COV"} = $cov ;
		$results{$sample}{$scaff}{$id}{"ID"} = $ident ;
		$results{$sample}{$scaff}{$id}{"LENGTH"} = $length ;
	}
	
	close IN ;
	
	return \%results ;
}

sub read_gff {
	my$file = $_[0] . "/" . $_[1] ;
    open GFF, "<$file" or die "can not open gff $file\n" ;

    my%annotation ;
    my%Sv_sym_matches ;

    while (<GFF>) {

        if ( $_ =~ m/^#/ ) {next ;}

        chomp $_ ;
        my @split = split (/\t/, $_) ;

        if ( $split[0] eq "Sv_mito_chromosome" ) { next ; }

#        if ( $split[2] ne "CDS" ) { next ; }

        my $id = "" ;
        if ( $split[8] =~ m/ID=(.+);Name=/ ) {$id = $1 ;}

        my $gene = "" ;
        if ( $split[8] =~ m/;Name=(.+);Dbxref/ ) {$gene = $1 ;}
        $gene =~ s/_OS=.+$// ;
        $gene =~ s/_Tax=.+$// ;
        $gene =~ s/_n=\d+// ;
        $gene =~ s/\[.+\]// ;
        $gene =~ s/_$// ;
		
        $annotation{$split[0]}{$id}{"LINE"} = $_ ;
        $annotation{$split[0]}{$id}{"START"} = $split[3] ;
        $annotation{$split[0]}{$id}{"STOP"} = $split[4] ;
        $annotation{$split[0]}{$id}{"STRAND"} = $split[6] ;
        $annotation{$split[0]}{$id}{"GENE"} = $gene ;

		if ($Sv_sym_matches{$split[0]}) {
			if ($split[8] =~ m/Solemya_velum_gill_symbiont/) {$Sv_sym_matches{$split[0]} ++ ;}
		}
		else {
			if ($split[8] =~ m/Solemya_velum_gill_symbiont/) {$Sv_sym_matches{$split[0]} = 1 ;}
			else {$Sv_sym_matches{$split[0]} = 0 ;}
		}

       }
    close GFF ;
    return (\%annotation, \%Sv_sym_matches) ;
}

sub parse_tabs {
	my$file = $_[0] ;
	
	my%output ;
	
	open IN, "<$file" or die "cannot open $file\n" ;
	
	while (<IN>) {
		chomp ;
		
		my@split = split(/\t/, $_) ;
		
		my$id = $split[0] ;
		my$taxonomy = $split[1] ;
		my$rRNAsmsu = $split[2] ;
		
		$output{$id}{"NAME"} = $taxonomy ;
		$output{$id}{"16SID"} = $rRNAsmsu ;
	}
	
	close IN ;
	
	return \%output ;

}
