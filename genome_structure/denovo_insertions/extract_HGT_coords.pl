use strict ;
use warnings ;

my$genes = $ARGV[0] ;
my$gff_dir = $ARGV[1] ;

## Get list of genes to be extracted
my%genes ;
open GENES, "<$genes" or die "cannot open $genes\n" ;

while (<GENES>) {
	chomp ;
	
	my@split = split(/\t/, $_) ;
	my@list = split(", ", $split[1]) ;
	
	$genes{$split[0]} = \@list ;
}

close GENES ;	
	


## Extract gene coordinates from gffs
opendir(GFFS, $gff_dir) ;
my@gffs = readdir(GFFS) ;

foreach my$file (@gffs) {

	if ($file !~ m/.gff$/) {next ;}
	
	my$sample = $file ;
	$sample =~ s/(\w+\d+)gill.+/$1/ ;
	
	if (! $genes{$sample}) {next ;}
	
	my%annotation ;
	
	open IN, "<$gff_dir/$file" or die "cannot open $file\n" ;
	
	while (<IN>) {
		if ($_ =~ m/^#/) {next ;}
		chomp ;
		
		my@split = split(/\t/, $_) ;
		
		my$start = $split[3] ;
		my$end = $split[4] ;
		
		my$info = $split[8] ;
		my@info = split(";", $info) ;
		
		my$ID = $info[0] ;
		$ID =~ s/ID=// ;
		
		my$product = $info[1] ;
		$product =~ s/Name=// ;
		
		$annotation{$ID}{"START"} = $start ;
		$annotation{$ID}{"END"} = $end ;
		$annotation{$ID}{"PRODUCT"} = $product ;
	}
			
	close IN ;
	
	open OUT, ">${sample}_HGT_genes_coords.txt" or die "cannot open ${sample}_HGT_genes_coords.txt\n" ;
	
	foreach my$gene (@{$genes{$sample}}) {
		print OUT $gene, "\t", $annotation{$gene}{"START"}, "\t", $annotation{$gene}{"END"}, "\t", $annotation{$gene}{"PRODUCT"}, "\n" ;
	}
	
	close OUT ;
}	