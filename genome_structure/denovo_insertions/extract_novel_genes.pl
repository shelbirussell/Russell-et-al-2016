use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

## Usage: ../extract_novel_genes.pl /Users/shelbirussell/Dropbox/lab/projects/svelum_projects/1-PopGen_transmission_project/synteny/novel_genes/results/PJRI36_vs_Sv_nt_novel_genes.out /Users/shelbirussell/Dropbox/lab/projects/svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/Svelum_syms/PJRI36gill_sym_filtered_draft_genome.gff /Users/shelbirussell/Dropbox/lab/projects/svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/Svelum_syms/PJRI36gill_sym_filtered_draft_genome.fasta

my$novel_genes_file = $ARGV[0] ;
my$gff = $ARGV[1] ;
my$fasta = $ARGV[2] ;

my$seq = read_fasta($fasta) ;
my$ann = read_gff($gff) ;

my%seq = %{$seq} ;
my%ann = %{$ann} ;

open IN, "<$novel_genes_file" or die "cannot open $novel_genes_file\n" ;

while (<IN>) {
	chomp ;
	
	my@split = split(/\t/, $_) ;
	
	my$scaff = $split[0] ;
	my$id = $split[1] ;
	
	my$start = $ann{$scaff}{$id}{"START"} ;
	my$stop = $ann{$scaff}{$id}{"STOP"} ;
	
	my$seq = substr($seq{$scaff}, $start-1, $stop-$start) ;
	$seq = uc($seq) ;
	
    if ($ann{$scaff}{$id}{"STRAND"} eq "-" ) {
        $seq = reverse_complement($seq) ;
    }

#	print length($seq{$scaff}), "\n" ;
    print ">", "$id", "\t", $scaff, ":", $start, "-", $stop, "\tlength:", $stop-$start, "\tstrand:", $ann{$scaff}{$id}{"STRAND"}, "\tproduct:", $ann{$scaff}{$id}{"PRODUCT"}, "\n" ;
	print "$_\n" foreach ($seq =~ /.{1,80}/g) ;
}

close IN ;



sub read_fasta {
    open FASTA, "<$_[0]" ;

    my%seqs ;
    my$header ;
    my$seq ;
	my@split ;
	
    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
                $seqs{$split[0]} = $seq ;
                }

            $header = $_ ;
            $header =~ s/^>// ;
            $header =~ s/\s+$// ;
			@split = split(/\t/, $header) ;
						
            $seq = "" ;
            }

        else {
            $_ =~ s/\s+//g ;
            $seq .= $_ ;
            }
        }

    close FASTA ;

    if ($seq) {
        $seqs{$split[0]} = $seq ;
    }

    return \%seqs ;
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

sub reverse_complement {
    my$dna = $_[0] ;

    my$revcomp = reverse($dna) ;

    $revcomp =~ tr/ACGTacgt/TGCAtgca/ ;

    return $revcomp ;
    }