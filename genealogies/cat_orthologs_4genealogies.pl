use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# usage: perl cat_orthologs.pl orthologs.txt Svsym_pop_seqs/ Selarraichensis_nt.fasta Spervernicosa_nt.fasta Svelesiana_nt.fasta ./output/

my$reciprocal = $ARGV[0] ;
my$pop_dir = $ARGV[1] ;
my$sp2_nt = $ARGV[2] ;
my$sp3_nt = $ARGV[3] ;
my$sp4_nt = $ARGV[4] ;
my$output_dir ;

if ($ARGV[5]) {$output_dir = $ARGV[5] ;}

else {$output_dir = "./output/"}

my($orthologs, $names) = get_reciprocal($reciprocal, $pop_dir, $sp2_nt, $sp3_nt, $sp4_nt) ;
my%orthologs = %{$orthologs} ;
my%names = %{$names} ;

my%sample_fastas ;
my%spp_fastas ;

foreach my$sp (keys %names) {
	if ($pop_dir =~ m/$sp/) {
		my@orthos = @{$names{$sp}} ;

		opendir(POP, $pop_dir) ;
		my@samples = readdir(POP) ;
		foreach my$sample (@samples) {
			if ($sample =~ m/.fasta/) {
				my$sample_path = $pop_dir . "/" . $sample ;
				my$fasta = read_fasta($sample_path, \@orthos) ;
				$sample =~ s/^(\w+\d+)gill_.+/$1/ ;
				$sample_fastas{$sample} = $fasta ;
			}
		}
	}	
		
	elsif ($sp2_nt =~ m/$sp/) {
		my@orthos = @{$names{$sp}} ;
		$spp_fastas{$sp} = read_fasta($sp2_nt, \@orthos) ;
	}

	elsif ($sp3_nt =~ m/$sp/) {
		my@orthos = @{$names{$sp}} ;
		$spp_fastas{$sp} = read_fasta($sp3_nt, \@orthos) ;
	}
	
	elsif ($sp4_nt =~ m/$sp/) {		
		my@orthos = @{$names{$sp}} ;
		$spp_fastas{$sp} = read_fasta($sp4_nt, \@orthos) ;
	}
	
	else {
		print $sp, " doesn't match\n" ;
	}
}

foreach my$gene (nsort keys %orthologs) {
	open OUT, ">${gene}_multi.fasta" or die "can't open ${gene}_multi.fasta\n" ;
	
	foreach my$sample (nsort keys %sample_fastas) {
		my%record = %{$sample_fastas{$sample}} ;
		my$seq = $record{$gene} ;
		if ($seq) { 
#			print $gene, " ", $sample, " ", $seq, "\n" ;
			print OUT ">$sample\t$gene\n" ;
			print OUT "$_\n" foreach ($seq =~ /.{1,80}/g) ;
		}
		else {
			print "No variants at: ", $gene, " ", $sample, "\n" ;
		}
	}
	
	foreach my$sp (nsort keys %spp_fastas) {
		my$sp_gene = $orthologs{$gene}{$sp} ;
		my%record = %{$spp_fastas{$sp}} ;
		my$seq = $record{$sp_gene} ;
#		print $gene, " ", $sp, " ", $seq, "\n" ;
		print OUT ">$sp\t$sp_gene\n" ;
		print OUT "$_\n" foreach ($seq =~ /.{1,80}/g) ;
	}
	close OUT ;
}
	

sub read_fasta {
    open FASTA, "<$_[0]" or die "can't open $_[0]\n" ;
	my@orthologs = @{$_[1]} ;
	
    my%seqs ;
    my$header ;
    my$seq ;

    while (<FASTA>) {

        chomp $_ ;

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
            	if (grep /$header/, @orthologs) {
	                $seqs{$header} = $seq ;
	            }
            }

            $header = $_ ;
            $header =~ s/^>([a-zA-Z0-9_-]+)\s+.+$/$1/ ;
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

sub get_reciprocal {
	
	my%names ;
	
	my$sp1 = $_[1] ;
	$sp1 =~ s/^.+\/(Sv_sym_genes)\/$/$1/ ;
	@{$names{$sp1}} = () ;
	
	my$sp2 = basename($_[2]) ;
	$sp2 =~ s/^(\w+)_.+/$1/ ;
	@{$names{$sp2}} = () ;
	
	my$sp3 = basename($_[3]) ;
	$sp3 =~ s/^(\w+)_.+/$1/ ;
	@{$names{$sp3}} = () ;
	
	my$sp4 = basename($_[4]) ;
	$sp4 =~ s/^(\w+)_.+/$1/ ;
	@{$names{$sp4}} = () ;

	open RECIP, "<$_[0]" or die "cannot open $_[0]\n" ;
	
	my%orthologs ;
	while (<RECIP>) {
    	chomp $_ ;
    	if ($_ =~ m/^#/) { next ; }
    	my@split = split(/\t/, $_) ;
    	
    	$orthologs{$split[0]}{$sp2} = $split[1] ;
    	$orthologs{$split[0]}{$sp3} = $split[2] ;
    	$orthologs{$split[0]}{$sp4} = $split[3] ;

		push @{$names{$sp1}}, $split[0] ;
		push @{$names{$sp2}}, $split[1] ;
		push @{$names{$sp3}}, $split[2] ;
		push @{$names{$sp4}}, $split[3] ;
		
	}
	close RECIP ;
	
	return \%orthologs, \%names ;
}
