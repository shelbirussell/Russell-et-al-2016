use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

### The purpose of this script is to use gene annotation information to filter out scaffolds that are likely host contamination
### usage: perl filter_scaffolds_by_annotation.pl Selarraichensis_extracted-15x_50%GC_sym_scaff.fasta Selarraichensis_sym_genome.gff

my$assembly = $ARGV[0] ;
my$gff = $ARGV[1] ;

my$sample = basename($assembly) ;
$sample =~ s/_.+$// ;

#Eukaryotic taxa to remove
my@euks = ("Crassostrea_","Branchiostoma_","Saccoglossus_","Strongylocentrotus_","Nematostella_","Daphnia_","Danio_","Hydra_","Tetraodon_","Anoplophora_","Rhizophagus_","Capitella_","Ostrea_","Oreochromis_","Carassius_","Stegodyphus_","Amphimedon_","Lottia_","Caenorhabditis_","Metaseiulus_","Solenopsis_","Chaetopterus_","Azumapecten_","Acyrthosiphon_","Poecilia_","Poecilia_","Anolis_","Rhipicephalus_","Bison_","Trichuris_","Aplysia_","Callorhynchus_","Ixodes_","Ancylostoma_","Conothele_","Drosophila_","Limulus_","Reticulitermes_","Larimichthys_","Latimeria_","Nasonia_","Oncorhynchus_","Xiphophorus_","Acropora_","Biomphalaria_","Ophiophagus_","Xenopus_","Lygus_","Lytechinus_","Camponotus_","Rhizopus_","Diaphorina_citri_", "Bombyx_", "Pelodiscus_", "Danaus_", "Helobdella_", "Arabis_", "Oesophagostomum_", "Ciona_", "Stegastes_", "Tribolium_", "Trichoplax_", "Pachygrapsus_", "Tupaia_", "Eimeria_", "Cricetulus_", "Onchocerca_", "Platynereis_", "Dendroctonus_", "Geospiza_") ;

#Get scaffold names,lengths and GC content, get annotations
my$seqs = read_fasta($assembly) ;
my$annotations = read_gff($gff) ;

my%seqs = %{$seqs} ;
my%annotations = %{$annotations} ;

my$filtered_output = $sample . "_sym_removed_scaffolds.out" ;
my$flagged_output = $sample . "_sym_flagged_scaffolds.out" ;
open REMOVED, ">$filtered_output" or die "can't open $filtered_output\n" ;
print REMOVED "scaffold\tlength\t%GC\treason_filtered\n" ;
open FLAGGED, ">$flagged_output" or die "can't open $flagged_output\n" ;
print FLAGGED "scaffold\tlength\t%GC\treason_flagged\n" ;

#List of scaffolds to keep in assembly
my@keep ;

#Add scaffolds containing genes with homology to list to keep (i.e. remove scaffolds with no genes and scaffolds with only prodigal-predicted cds (maybe with RNAs) that didn't blast to anything)
foreach my$scaff (nsort keys %seqs) {
    if ( ! $annotations{$scaff}) {
        print REMOVED $scaff, "\t", $seqs{$scaff}{"LENGTH"}, "\t", $seqs{$scaff}{"GC"}, "\tno_genes\n" ;
        next ;
    }

    else {
        my$hypothetical_flag ;
        my$euk_flag ;
        my$chimera_flag ;

        foreach my$start (keys %{$annotations{$scaff}}) {
            if ($annotations{$scaff}{$start}{"PRODUCT"} =~ m/hypothetical_protein/ && $annotations{$scaff}{$start}{"SOURCE"} =~ m/prodigal/) {
                next ;
            }

            elsif ($annotations{$scaff}{$start}{"SOURCE"} =~ m/RNA/) {
                next ;
            }

            else {
                $hypothetical_flag = "keep" ;
            }

            if ($annotations{$scaff}{$start}{"TAXON"} ~~ @euks) {
                $chimera_flag = "chimera" ;
                next ;
            }

            else {
                print $scaff, " ", $annotations{$scaff}{$start}{"TAXON"}, "\n" ;
                $euk_flag = "keep" ;
            }
        }

        if ($hypothetical_flag && $euk_flag) {
            push @keep, $scaff ;
        }

        if ($chimera_flag && $euk_flag) {
            print FLAGGED $scaff, "\t", $seqs{$scaff}{"LENGTH"}, "\t", $seqs{$scaff}{"GC"}, "\tbacterial_euk_chimera\n" ;
        }

        if (! $hypothetical_flag) {
            print REMOVED $scaff, "\t", $seqs{$scaff}{"LENGTH"}, "\t", $seqs{$scaff}{"GC"}, "\tonly_predictions\n" ;
        }

        if (! $euk_flag) {
            print REMOVED $scaff, "\t", $seqs{$scaff}{"LENGTH"}, "\t", $seqs{$scaff}{"GC"}, "\teukaryotic\n" ;
        }
    }
}

close REMOVED ;
close FLAGGED ;

my$output_gff = $sample . "_sym_filtered.gff" ;
my$output_assembly = $sample . "_sym_filtered.fasta" ;

open OUT_GFF, ">$output_gff" or die "can't open $output_gff\n" ;
open OUT_FASTA, ">$output_assembly" or die "can't open $output_assembly\n" ;

foreach my$scaff (nsort @keep) {
    my$count = scalar( keys %{$annotations{$scaff}}) ;
    my$gff_scaff_header = "#" . $scaff . "\tlength:" . $seqs{$scaff}{"LENGTH"} . "\t%GC:" . $seqs{$scaff}{"GC"} . "\t#genes:" . $count ;
    print OUT_GFF $gff_scaff_header . "\n" ;

    foreach my$start (sort {$a <=> $b} keys %{$annotations{$scaff}}) {
        print OUT_GFF $annotations{$scaff}{$start}{"LINE"}, "\n" ;
    }

    print OUT_FASTA ">" . $seqs{$scaff}{"HEADER"}, "\n" ;
    print OUT_FASTA "$_\n" foreach ($seqs{$scaff}{"SEQ"} =~ /.{1,80}/g) ;
}

close OUT_GFF ;
close OUT_FASTA ;

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
        my@info = split(";", $split[8]) ;
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

        $gff{$scaffold}{$start}{"LINE"} = $_ ;
        $gff{$scaffold}{$start}{"TYPE"} = $type ;
        $gff{$scaffold}{$start}{"PRODUCT"} = $product ;
        $gff{$scaffold}{$start}{"TAXON"} = $taxon ;

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

sub read_fasta {
    open FASTA, "<$_[0]" or die "can't open fasta file: $_[0]" ;

    my%seqs ;
    my$header ;
    my$scaffold ;
    my$seq ;

    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
                $seqs{$scaffold}{"HEADER"} = $header ;
                $seqs{$scaffold}{"SEQ"} = $seq ;
                }

            $header = $_ ;
            $header =~ s/^>// ;
            $header =~ s/\s+$// ;
            $header =~ m/^(scaffold_\d+)\s+.+/ ;
            $scaffold = $1 ;

            $seq = "" ;
            }

        else {
            $_ =~ s/\s+//g ;
            $seq .= $_ ;
            }
        }

    close FASTA ;

    if ($seq) {
        $seqs{$scaffold}{"HEADER"} = $header ;
        $seqs{$scaffold}{"SEQ"} = $seq ;
        }

    my%seq_info ;
    foreach my$scaff (keys %seqs) {

        my$length = length($seqs{$scaff}{"SEQ"}) ;

        my$C = $seqs{$scaff}{"SEQ"} =~ tr/C|c// ;
        my$G = $seqs{$scaff}{"SEQ"} =~ tr/G|g// ;
        my$A = $seqs{$scaff}{"SEQ"} =~ tr/A|a// ;
        my$T = $seqs{$scaff}{"SEQ"} =~ tr/T|t// ;

        my$GC = ($C + $G) / ($C + $G + $A + $T) ;

        $seqs{$scaff}{"LENGTH"} = $length ;
        $seqs{$scaff}{"GC"} = $GC ;
    }

    return \%seqs ;
    }