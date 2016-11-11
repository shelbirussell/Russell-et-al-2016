use strict ;
use warnings ;

# usage: perl bin_bins.pl 10 Sv_sym_ALL-LD_bin.txt

my$bin_size = $ARGV[0] ;
my$input = $ARGV[1] ;
my$output = $input ;
$output =~ s/bin.txt/bin$bin_size.txt/ ;

my$counter = 0 ;
my$ld = 0 ;
my$last_bin ;

open IN, "<$input" ;
open OUT, ">$output" ;
my@split ;

while (<IN>) {
	
	chomp ;
	@split = split(/\t/, $_) ;

	$counter ++ ;
	$ld += $split[1] ;	
	
	if ($counter == $bin_size) {
		print OUT $split[0], "\t", $ld/$bin_size, "\n"; 
		
		$counter = 0 ;
		$ld = 0 ;
	}	
	
}	

print OUT $split[0], "\t", $ld/$counter, "\n" ;

close IN ;
close OUT ;