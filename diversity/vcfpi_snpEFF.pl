use strict ;
use warnings ;

## Uses snpEFF CS.vcf file to include all genotyped sites

# usage: perl vcfpi_snpEFF.pl Sv_sym_scaffold SvDML_MSGATKp1_snp.ann

my$chromosome_prefix = $ARGV[0] ;
my$vcf = $ARGV[1] ;

open VCF, "<$vcf" or die "can't open $vcf\n" ;

my$samples = 0 ;
my$segregating = 0 ;
my$indels = 0 ;
my$snps = 0 ;
my%alleles ;
my$length = 0 ;

my$synonymous_pi = 0 ;
my$synonymous_count = 0 ;
my$nonsynonymous_pi = 0 ;
my$nonsynonymous_count = 0 ;
my$intergenic_count = 0 ;
my$intergenic_pi = 0 ;
my$stop_count = 0 ;
my$stop_pi = 0 ;
my$start_count = 0 ;
my$start_pi = 0 ;
my$indel_effects = 0 ;

while (<VCF>) {
	
	if ($_ =~ m/#CHROM/) {
		my@header = split(/\t/, $_) ;
		
		$samples = $#header - 8 ;
	}
		
	if ($_ =~ m/#/) {next ;}
		
	my@split = split(/\t/, $_) ;
	
# Exclude chromosomes from other genomes (e.g. mitochondrial)
	if ($split[0] !~ m/$chromosome_prefix/) {
		next ;
	}
	
	$length ++ ;

# Skip sample-wide invariant sites
	if ($split[4] =~ m/\./) {
		next ;
	}

# Filter out sites with low support
	my@info = split(/;/, $split[7]) ;
	my$AN = $info[2] ;	
	$AN =~ s/AN=// ;

	if ($AN < $samples / 2) {
		next ;
	}

# Filter out invariant sites (retained in extracted samples, when variant was in another pop)
	my$AC = $info[0];
	$AC =~ s/AC=// ;
	my@ACs = split(/,/, $AC) ;

	if (scalar @ACs < 2 && $ACs[0] == 0) {
		next ;
	}

	if ($ACs[0] == $AN ) {
		next ;
	}

# Calculate pi for each allele
	my$site_pi ;	
	
	# Calculate reference allele count from other counts
	my$ref_count = $AN;
	foreach (@ACs) {
		$ref_count -= $_ ;
	}

	push @ACs, $ref_count ;
	
	# Biallelic pi
	# over all sites sum( 2*j(n-j)/ n(n-1) ) / L, where j = alternate allele count, n = number of called samples at site, and L=number of callable sites
	if (scalar @ACs == 2) {
		my$j = $ACs[0] ;
		$site_pi = 2 * $j * ($AN - $j) / ($AN * ($AN - 1)) ;		
	}
		
	# Multiallelic pi
	# over all site sum( over all alleles sum_i( Ji (n-Ji) ) / (n(n-1)) ) / L
	my$sum_i = 0 ;
	foreach my$allele (@ACs) {
		$sum_i += ($allele * ($AN-$allele)) ;
	}
	$site_pi = $sum_i / ($AN*($AN-1)) ;

# Read in snpEFF reports to classify alleles
	my$snpeff_reports = $info[-1] ;

	my@snpeff = split(/\|/, $snpeff_reports) ;
	my$type = $snpeff[1] ;
		
	if ($type eq "synonymous_variant") {
		$synonymous_count ++ ;
		$synonymous_pi += $site_pi ;
	}
		
	elsif (grep(/stream/, $type)){
		$intergenic_count ++ ;
		$intergenic_pi += $site_pi ;
	}
		
	elsif ($type eq "missense_variant") {
		$nonsynonymous_count ++ ;
		$nonsynonymous_pi += $site_pi ;
	}
	
	elsif (grep(/stop/, $type)) {
		$stop_count ++ ;
		$stop_pi += $site_pi ;
	}
	
	elsif (grep(/initiator/, $type)) {
		$start_count ++ ;
		$start_pi += $site_pi ;
	}
	
	elsif (grep(/frameshift/, $type)) {
		$indel_effects ++ ;
	}
	
	else {
		print $type, "\n" ;
	}
}

close VCF ;

my$avg_synonymous_pi = $synonymous_pi / $length ;
my$avg_nonsynonymous_pi = $nonsynonymous_pi / $length ;
my$avg_intergenic_pi = $intergenic_pi / $length ;
my$avg_stop_pi = $stop_pi / $length ;
my$avg_start_pi = $start_pi / $length ;

print "number of samples: ", $samples, "\n" ;
print "confidently called sites: ", $length, "\n" ;
print "synonymous variants: ", $synonymous_count, "\n" ;
print "synonymous pi: ", $avg_synonymous_pi, "\n" ;
print "nonsynonymous variants: ", $nonsynonymous_count, "\n" ;
print "nonsynonymous pi: ", $avg_nonsynonymous_pi, "\n" ;
print "intergenic variants: ", $intergenic_count, "\n" ;
print "intergenic pi: ", $avg_intergenic_pi, "\n" ;
print "stop variant count: ", $stop_count, "\n" ;
print "stop pi: ", $avg_stop_pi, "\n" ;
print "start variant count: ", $start_count, "\n" ;
print "start pi: ", $avg_start_pi, "\n" ;
print "indel effects: ", $indel_effects, "\n" ;
