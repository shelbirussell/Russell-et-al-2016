use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# Usage: reorder_mauve_backbone.pl .backbone 
# samples after backbone file have to be in order as they are in the file (e.g. sequence1 = DML1)

my$backbone_file = $ARGV[0] ;
my$output = basename($backbone_file) ;
$output =~ s/.backbone/_reordered.backbone/ ;

my@samples = ("DML1", "DML2", "DML4", "DML5", "DML8", "DML9", "DML10", "DML12", "DML13", "DML14", "DML15", "DML16", "DML17", "DML18", "DML22", "KB16", "KB18", "KB25", "KB27", "KB29", "KB31", "KB34", "KB36", "PJRI36", "PJRI37", "PJRI38", "PJRI39", "PJRI44", "PJRI47", "PJRI48", "PJRI50", "PJRI51", "PJRI53", "PJRI55", "PJRI56", "PJRI57", "PJRI58", "PJRI59", "PJRI60", "PJRI61", "PJRI62", "PJRI63", "SB4", "SB7", "SB8", "SB11", "SB13", "SB16", "SB24", "SB25", "SB27", "SB28", "SRI2", "SRI3", "SRI4", "SRI6", "SRI9", "SRI12", "SRI16", "Svsym_ref") ;
my@reorder = ("KB16", "KB18", "KB25", "KB27", "KB29", "KB31", "KB34", "KB36", "PJRI36", "PJRI37", "PJRI38", "PJRI39", "PJRI44", "PJRI47", "PJRI48", "PJRI50", "PJRI51", "PJRI53", "PJRI55", "PJRI56", "PJRI57", "PJRI58", "PJRI59", "PJRI60", "PJRI61", "PJRI62", "PJRI63", "SRI2", "SRI3", "SRI4", "SRI6", "SRI9", "SRI12", "SRI16", "SB4", "SB7", "SB8", "SB11", "SB13", "SB16", "SB24", "SB25", "SB27", "SB28", "DML1", "DML2", "DML4", "DML5", "DML8", "DML9", "DML10", "DML12", "DML13", "DML14", "DML15", "DML16", "DML17", "DML18", "DML22", "Svsym_ref") ;

my%bbone ;
my@header = () ;
foreach my$sample (@samples) {
	$bbone{"${sample}_left"} = () ;
	$bbone{"${sample}_right"} = () ;
	push @header, "${sample}_left" ;
	push @header, "${sample}_right" ;
}

my@new_header = () ;
foreach my$sample (@reorder) {
	push @new_header, "${sample}_left" ;
	push @new_header, "${sample}_right" ;
}

my$lines = 0 ;
	
open IN, "<$backbone_file" or die "cannot open $backbone_file\n" ;

while (<IN>) {
	chomp ;
	
	if ($_ =~ m/seq/) {next;}
	
	$lines ++ ;
	
	my@split = split(/\t/, $_) ;
	
	foreach my$i (0..$#split) {
		push @{$bbone{$header[$i]}}, $split[$i] ;}

}

close IN ;

open OUT, ">$output" or die "cannot open $output\n" ;
foreach my$sample (@new_header) {
	print OUT $sample, "\t" ;
}

print OUT "\n" ;

foreach my$i (0..$lines) {
	foreach my$sample (@new_header) {
		print OUT $bbone{$sample}[$i], "\t" ;
	}
	
	print OUT "\n" ;
	
}

close OUT ;