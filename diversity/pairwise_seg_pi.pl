use strict ;
use warnings ;
use Sort::Versions ;

## Uses VS.vcf file to include all variant sites

# Returns: pairwise genotypes  
# usage: perl pairwise_pi.pl Sv_sym_scaffold 10 SvDML_MS_stampy_UGploidy1_CSboth.filtered.vcf

my$chromosome_prefix = $ARGV[0] ;
my$window = $ARGV[1] ;
my$vcf = $ARGV[2] ;

open VCF, "<$vcf" or die "can't open $vcf\n" ;

my%names ;
my%gts ;

while (<VCF>) {
	
	if ($_ =~ m/#CHROM/) {
		my@header = split(/\t/, $_) ;
			
		foreach (9 .. $#header) {
			my$index = $_ - 9 ;
			$names{$index} = $header[$_] ;
		}
	}
		
	if ($_ =~ m/#/) {next ;}
		
	my@split = split(/\t/, $_) ;
	
# Exclude chromosomes from other genomes (e.g. mitochondrial)
	if ($split[0] !~ m/$chromosome_prefix/) {next ;}

# Skip sample-wide invariant sites
	if ($split[4] =~ m/\./) {next ;}

# Filter out sites with low support
	my@info = split(/;/, $split[7]) ;
	my$AN = $info[2] ;	
	$AN =~ s/AN=// ;

	my$samples = $#split - 8 ;
	
	if ($AN < $samples / 2) {next ;}

# Filter out invariant sites (retained in extracted samples, when variant was in another pop)
	my$AC = $info[0];
	$AC =~ s/AC=// ;
	my@ACs = split(/,/, $AC) ;

	if (scalar @ACs < 2 && $ACs[0] == 0) {next ;}

	if ($ACs[0] == $AN ) {next ;}
	
# Get scaffold position	
	my$scaff = $split[0] ;
	$scaff =~ s/Sv_sym_scaffold// ;
	my$coord = $scaff . "\t" . $split[1] ;

# Get genotypes
	foreach my$id (9 .. $#split) {
		my$index = $id - 9 ;

		my$gt = $split[$id] ;
		
		if ($gt =~ m/^(\d):/) {
			$gt = $1 ;
		}
		
		else {
			$gt = "NA" ;
		}
		
		$gts{$coord}{$index} = $gt ;
	}
}

close VCF ;

my@names = keys %names ;
my@coord = sort{ versioncmp( $a, $b ) } keys %gts ; 
#print join(",", @coord), "\n" ;


foreach my$sample1 (0 .. $#names-1) {
	foreach my$sample2 ($sample1+1 .. $#names) {
		$names{$sample1} =~ m/^(\w+\d+)/ ;
		my$name1 = $1 ;
		$names{$sample2} =~ m/^(\w+\d+)/ ;
		my$name2 = $1 ;
		print $name1, "\t", $name2, "\t" ;
		
		open OUT, ">${name1}_${name2}.out" ;
		my$counter = 0 ;
		my$total_count = 0 ;
		my$total_pi = 0 ;
		my$pi = 0 ;
		my$last ;
		
		foreach my$pos (@coord) {
			
			$last = $pos ;

			if ($gts{$pos}{$sample1} eq "NA" || $gts{$pos}{$sample2} eq "NA") {
				next ;
			}
			
			elsif ($gts{$pos}{$sample1} eq $gts{$pos}{$sample2}) {
				$counter ++ ;
			}
			
			else {
				$pi ++ ;
				$counter ++ ;
			}
			
			if ($counter == $window) {
				print OUT $pos, "\t", $pi/$window, "\n" ;
				$total_count += $counter ;
				$counter = 0 ;
				$total_pi += $pi ;
				$pi = 0 ;
			}	
		}
		
		if ($counter < $window && $counter > 0) {
			$total_count += $counter ;
			print OUT $last, "\t", $pi/$counter, "\n" ;
		}
		close OUT ;

		print "average pi: ", $total_pi/$total_count, "\n" ;

	}
}