use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

# This script formats link info for circos
# Takes: ref.gff, sym.gff, reordered.fasta, multireciprocal.out	
# And outputs 
# 1. karyotype info: in reordered order, each scaffold on a line, followed by its count, start, end, and color
# 2. link info: 

## Usage: perl format_synteny4circos.pl \
#/Users/shelbirussell/Dropbox/lab/reference_seqs/Svelum/Sv_sym/files/Solemya_velum_Sym_top4.gff \
#/Users/shelbirussell/Dropbox/Lab/projects/Svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/reorder/reordered_denovo_assemblies/reordered_scaffolds/KB16gill_sym_filtered_draft_genome_reordered_gt10000.gff \
#/Users/shelbirussell/Dropbox/Lab/projects/Svelum_projects/Genome_processing/denovo_assembly/FINAL_illumina_assemblies/reorder/reordered_denovo_assemblies/reordered_scaffolds/KB16gill_sym_filtered_draft_genome_reordered_gt10000.fasta \
#/Users/shelbirussell/Dropbox/lab/projects/svelum_projects/Genome_processing/orthology/Svelum_ref_orthologs/reciprocal_best_hits/KB31_multireciprocal.out \
#78,0,222

my$ref_gff = $ARGV[0] ;
my$sample_gff = $ARGV[1] ;
my$reordered_fasta = $ARGV[2] ;
my$multireciprocal = $ARGV[3] ;
my$color = $ARGV[4] ;

my$ref_anno = read_Svsym_gff($ref_gff) ;
my$anno = read_gff($sample_gff) ;
my$reordered_headers = read_fasta_headers($reordered_fasta) ;
my$orthologs = read_orthologs($multireciprocal) ;

my%ref_anno = %{$ref_anno} ;
my%anno = %{$anno} ;
my@reordered_headers = @{$reordered_headers} ;
my%orthologs = %{$orthologs} ;

my$reference_karyo = get_ref_karyo() ;

my@short_scfs = () ;

my$karyotype = basename($sample_gff) ;
$karyotype =~ s/.gff/.karyotype.txt/ ;
open KARYO, ">$karyotype" or die "cannot open $karyotype\n" ;

print KARYO $reference_karyo ;

my$count = 0 ;

foreach my$i (reverse 0..$#reordered_headers) {
	my@split = split(/\t/, $reordered_headers[$i]) ;
	if ($split[1] <= 10000) {
		push @short_scfs, $split[0] ;
		next ;
	}
	$count ++ ;
	print KARYO "chr - ", $split[0], "\t", $count, "\t", "0", "\t", $split[1]-1, "\t", $color, "\n" ;
}
close KARYO ;

my$links = basename($sample_gff) ;
$links =~ s/.gff/.links.txt/ ;
open LINKS, ">$links" or die "cannot open $links\n" ;
foreach my$scaff (nsort keys %ref_anno) {
	foreach my$id (nsort keys %{$ref_anno{$scaff}}) {
		if ($orthologs{$id}) {
			my$ref_id = $id ;
			my$sample_id = $orthologs{$id} ;
			
			# If using >10kb scaffold gff and fasta, not all ids in orthologs file are in sample gff
			if ($anno{$sample_id}) {
				my$ref_start = $ref_anno{$scaff}{$id}{"START"} - 1 ;
				my$ref_stop = $ref_anno{$scaff}{$id}{"STOP"} - 1 ;
			
				my$sample_scaff = $anno{$sample_id}{"SCAFFOLD"} ;
				my$sample_start = $anno{$sample_id}{"START"} - 1;
				my$sample_stop = $anno{$sample_id}{"STOP"} - 1;
			
				foreach my$short (@short_scfs) {
					if ($sample_scaff eq $short) { 
						next ;
					}
				}
			
				print LINKS $scaff, "\t", $ref_start, "\t", $ref_stop, "\t", $sample_scaff, "\t", $sample_start, "\t", $sample_stop, "\n" ;
			}
		}
	}
}

close LINKS ;



sub read_orthologs {
	my$reciprocal_hits = $_[0] ;
	my%reciprocals ;

	open IN, "<$reciprocal_hits" or die "cannot open $reciprocal_hits \n" ;

	while (<IN>) {
		chomp ;
		my@split = split(/\t/, $_) ;
		$reciprocals{$split[0]} = $split[1] ;
	}

	close IN ;
	
	return \%reciprocals ;
}

sub read_fasta_headers {
    open FASTA, "<$_[0]" ;

	my@headers = () ;
    my$header ;
    my$seq ;
	my@split ;
	
    while (<FASTA>) {

        if ( $_ =~ m/^#/ ) {
            next ;
            }

        if ( $_ =~ m/>/ ) {
            if ($seq) {
		    	my$length = length($seq) ;
				push @headers, "$split[0]\t$length" ;
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
    	my$length = length($seq) ;
		push @headers, "$split[0]\t$length" ;
    }

    return \@headers ;
}

sub read_gff {
    my$cds = $_[0] ;
    open SAMPLE_CDS, "<$cds" or die "can't open gff file: $_[0]" ;

    my%gff ;
		
    while (<SAMPLE_CDS>) {
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
		
        $gff{$id}{"SCAFFOLD"} = $scaffold ;
        $gff{$id}{"LINE"} = $_ ;
        $gff{$id}{"START"} = $start ;
        $gff{$id}{"STOP"} = $stop ;
        $gff{$id}{"TYPE"} = $type ;
        $gff{$id}{"PRODUCT"} = $product ;
        $gff{$id}{"TAXON"} = $taxon ;
		$gff{$id}{"STRAND"} = $strand ;
		
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

    close SAMPLE_CDS ;

    return \%gff ;

}

sub read_Svsym_gff {
    my$gff = $_[0] ;
	open REF_GFF, "<$gff" ;

    my %gff ;
	my%counts ;

    while (<REF_GFF>) {

        if ( $_ =~ m/^#/ ) {
            next ;
        }

        chomp $_ ;
        my @split = split ( /\t/, $_ ) ;
        my$scaffold = $split[0] ;

        if ( $scaffold eq "Sv_mito_chromosome" ) { next ; }

		if ($counts{$scaffold}) {$counts{$scaffold} ++ ;}
		
		else {$counts{$scaffold} = 1 ;}

        my $id = "" ;
        if ( $split[8] =~ m/locus_tag=(.+)(;product|;note)/ ) {
            $id = $1 ;
        }
        if ( $split[8] =~ m/locus_tag=(.+);note/ ) {
            $id = $1 ;
        }

        my $gene = "" ;
        if ( $split[8] =~ m/product=([\(\)a\/,a-zA-Z0-9_-]+)/ ) {
            $gene = $1 ;
        }

        $gff{$scaffold}{$id}{"COUNT"} = $counts{$scaffold} ;
        $gff{$scaffold}{$id}{"START"} = $split[3] ;
        $gff{$scaffold}{$id}{"STOP"} = $split[4] ;
        $gff{$scaffold}{$id}{"STRAND"} = $split[6] ;
        $gff{$scaffold}{$id}{"PRODUCT"} = $gene ;

       }
    close REF_GFF ;
    return \%gff ;
    }

sub get_ref_karyo {
	my$karyo = "chr - Sv_sym_scaffold2 2 0 892554 lgrey\nchr - Sv_sym_scaffold4 4 0 28016 lgrey\nchr - Sv_sym_scaffold1 1 0 1213830 lgrey\nchr - Sv_sym_scaffold3 3 0 537612 lgrey\nband Sv_sym_scaffold1 chromosome chromosome 0 4191 lgrey\nband Sv_sym_scaffold1 JV_ISf_40 JV_ISf_40 4192 5232 black\nband Sv_sym_scaffold1 chromosome chromosome 5233 125751 lgrey\nband Sv_sym_scaffold1 JV46_18070-IS JV46_18070-IS 125752 126714 black\nband Sv_sym_scaffold1 chromosome chromosome 126715 139099 lgrey\nband Sv_sym_scaffold1 JV46_12410-IS JV46_12410-IS 139100 140062 black\nband Sv_sym_scaffold1 chromosome chromosome 140063 141546 lgrey\nband Sv_sym_scaffold1 JV46_18130-IS JV46_18130-IS 141547 142509 black\nband Sv_sym_scaffold1 chromosome chromosome 142510 152470 lgrey\nband Sv_sym_scaffold1 JV46_12510-IS JV46_12510-IS 152471 153277 black\nband Sv_sym_scaffold1 chromosome chromosome 153278 153318 lgrey\nband Sv_sym_scaffold1 JV46_12520-IS JV46_12520-IS 153319 153597 black\nband Sv_sym_scaffold1 chromosome chromosome 153598 193081 lgrey\nband Sv_sym_scaffold1 JV46_18520-IS JV46_18520-IS 193082 194251 black\nband Sv_sym_scaffold1 chromosome chromosome 194252 311399 lgrey\nband Sv_sym_scaffold1 JV46_19450-IS JV46_19450-IS 311400 312362 black\nband Sv_sym_scaffold1 chromosome chromosome 312363 321559 lgrey\nband Sv_sym_scaffold1 JV_ISf_3-IS JV_ISf_3-IS 321560 321661 black\nband Sv_sym_scaffold1 chromosome chromosome 321662 377919 lgrey\nband Sv_sym_scaffold1 JV46_13380-IS JV46_13380-IS 377920 378882 black\nband Sv_sym_scaffold1 chromosome chromosome 378883 425624 lgrey\nband Sv_sym_scaffold1 JV46_13580-IS JV46_13580-IS 425625 426587 black\nband Sv_sym_scaffold1 chromosome chromosome 426588 524403 lgrey\nband Sv_sym_scaffold1 JV46_13950-IS JV46_13950-IS 524404 525366 black\nband Sv_sym_scaffold1 chromosome chromosome 525367 535648 lgrey\nband Sv_sym_scaffold1 JV46_21090-IS JV46_21090-IS 535649 536611 black\nband Sv_sym_scaffold1 chromosome chromosome 536612 571678 lgrey\nband Sv_sym_scaffold1 JV46_14200-IS JV46_14200-IS 571679 572641 black\nband Sv_sym_scaffold1 chromosome chromosome 572642 783810 lgrey\nband Sv_sym_scaffold1 JV46_15980-IS JV46_15980-IS 783811 784773 black\nband Sv_sym_scaffold1 chromosome chromosome 784774 980434 lgrey\nband Sv_sym_scaffold1 JV46_22850-IS JV46_22850-IS 980435 981397 black\nband Sv_sym_scaffold1 chromosome chromosome 981398 1048754 lgrey\nband Sv_sym_scaffold1 JV46_17490-IS JV46_17490-IS 1048755 1049144 black\nband Sv_sym_scaffold1 chromosome chromosome 1049145 1213830 lgrey\nband Sv_sym_scaffold2 chromosome chromosome 0 38245 lgrey\nband Sv_sym_scaffold2 JV46_04480-IS JV46_04480-IS 38246 39052 black\nband Sv_sym_scaffold2 chromosome chromosome 39053 39093 lgrey\nband Sv_sym_scaffold2 JV46_04510-IS JV46_04510-IS 39094 39372 black\nband Sv_sym_scaffold2 chromosome chromosome 39373 42220 lgrey\nband Sv_sym_scaffold2 JV46_04640-IS JV46_04640-IS 42221 43027 black\nband Sv_sym_scaffold2 chromosome chromosome 43028 43062 lgrey\nband Sv_sym_scaffold2 JV46_04700-ISf_5 JV46_04700-ISf_5 43063 43347 black\nband Sv_sym_scaffold2 chromosome chromosome 43348 82737 lgrey\nband Sv_sym_scaffold2 JV46_08080-IS JV46_08080-IS 82738 83700 black\nband Sv_sym_scaffold2 chromosome chromosome 83701 87535 lgrey\nband Sv_sym_scaffold2 Scaff2_ICESpuPO1 Scaff2_ICESpuPO1 87536 93714 black\nband Sv_sym_scaffold2 JV46_07670-ICE38 JV46_07670-ICE38 93715 94904 black\nband Sv_sym_scaffold2 Scaff2_ICESpuPO1 Scaff2_ICESpuPO1 94905 97287 black\nband Sv_sym_scaffold2 JV46_07690-ICE37 JV46_07690-ICE37 97288 97750 black\nband Sv_sym_scaffold2 Scaff2_ICESpuPO1 Scaff2_ICESpuPO1 97751 98678 black\nband Sv_sym_scaffold2 JV46_07700-ICE38 JV46_07700-ICE38 98679 98993 black\nband Sv_sym_scaffold2 Scaff2_ICESpuPO1 Scaff2_ICESpuPO1 98994 99723 black\nband Sv_sym_scaffold2 chromosome chromosome 99724 102229 lgrey\nband Sv_sym_scaffold2 JV46_02330-IS JV46_02330-IS 102230 103036 black\nband Sv_sym_scaffold2 chromosome chromosome 103037 103077 lgrey\nband Sv_sym_scaffold2 JV46_02340-IS JV46_02340-IS 103078 103356 black\nband Sv_sym_scaffold2 chromosome chromosome 103357 109340 lgrey\nband Sv_sym_scaffold2 JV_ISf_40-IS JV_ISf_40-IS 109341 109567 black\nband Sv_sym_scaffold2 JV_ISf_40/26-IS JV_ISf_40/26-IS 109568 109625 black\nband Sv_sym_scaffold2 JV_ISf_26-IS JV_ISf_26-IS 109626 109660 black\nband Sv_sym_scaffold2 JV46_07730-ISf_26 JV46_07730-ISf_26 109661 110437 black\nband Sv_sym_scaffold2 JV46_07730-IS JV46_07730-IS 110438 110463 black\nband Sv_sym_scaffold2 JV46_07730/02410-IS JV46_07730/02410-IS 110464 110467 black\nband Sv_sym_scaffold2 JV46_02410-IS JV46_02410-IS 110468 110781 black\nband Sv_sym_scaffold2 chromosome chromosome 110782 110881 lgrey\nband Sv_sym_scaffold2 JV46_02420-ICE40 JV46_02420-ICE40 110882 111194 black\nband Sv_sym_scaffold2 chromosome chromosome 111195 111812 lgrey\nband Sv_sym_scaffold2 JV46_02430-ICE16 JV46_02430-ICE16 111813 112090 black\nband Sv_sym_scaffold2 chromosome chromosome 112091 112872 lgrey\nband Sv_sym_scaffold2 JV46_07740-IS JV46_07740-IS 112873 113823 black\nband Sv_sym_scaffold2 chromosome chromosome 113824 115038 lgrey\nband Sv_sym_scaffold2 JV46_02470-ICE21 JV46_02470-ICE21 115039 115468 black\nband Sv_sym_scaffold2 chromosome chromosome 115469 116330 lgrey\nband Sv_sym_scaffold2 JV46_02480-ICE16 JV46_02480-ICE16 116331 116589 black\nband Sv_sym_scaffold2 chromosome chromosome 116590 120136 lgrey\nband Sv_sym_scaffold2 JV46_02530-ICE39 JV46_02530-ICE39 120137 120434 black\nband Sv_sym_scaffold2 chromosome chromosome 120435 121040 lgrey\nband Sv_sym_scaffold2 JV46_02540-IS JV46_02540-IS 121041 121847 black\nband Sv_sym_scaffold2 chromosome chromosome 121848 121888 lgrey\nband Sv_sym_scaffold2 JV46_02550-IS JV46_02550-IS 121889 122167 black\nband Sv_sym_scaffold2 chromosome chromosome 122168 199798 lgrey\nband Sv_sym_scaffold2 JV46_03080-IS JV46_03080-IS 199799 200605 black\nband Sv_sym_scaffold2 chromosome chromosome 200606 200646 lgrey\nband Sv_sym_scaffold2 JV46_03090-IS JV46_03090-IS 200647 200925 black\nband Sv_sym_scaffold2 chromosome chromosome 200926 230605 lgrey\nband Sv_sym_scaffold2 JV46_03320-IS JV46_03320-IS 230606 231568 black\nband Sv_sym_scaffold2 chromosome chromosome 231569 232274 lgrey\nband Sv_sym_scaffold2 JV46_08160-IS JV46_08160-IS 232275 232553 black\nband Sv_sym_scaffold2 chromosome chromosome 232554 232594 lgrey\nband Sv_sym_scaffold2 JV46_08170-IS JV46_08170-IS 232595 233401 black\nband Sv_sym_scaffold2 chromosome chromosome 233402 415937 lgrey\nband Sv_sym_scaffold2 JV46_09020-IS JV46_09020-IS 415938 416216 black\nband Sv_sym_scaffold2 chromosome chromosome 416217 416257 lgrey\nband Sv_sym_scaffold2 JV_ISf_31-IS JV_ISf_31-IS 416258 417058 black\nband Sv_sym_scaffold2 chromosome chromosome 417059 418607 lgrey\nband Sv_sym_scaffold2 JV46_09060-IS JV46_09060-IS 418608 419159 black\nband Sv_sym_scaffold2 chromosome chromosome 419160 419162 lgrey\nband Sv_sym_scaffold2 JV46_04570-IS JV46_04570-IS 419163 419969 black\nband Sv_sym_scaffold2 chromosome chromosome 419970 420010 lgrey\nband Sv_sym_scaffold2 JV46_04580-IS JV46_04580-IS 420011 420289 black\nband Sv_sym_scaffold2 chromosome chromosome 420290 474957 lgrey\nband Sv_sym_scaffold2 JV46_09310-IS JV46_09310-IS 474958 475920 black\nband Sv_sym_scaffold2 chromosome chromosome 475921 478006 lgrey\nband Sv_sym_scaffold2 JV46_04930-IS JV46_04930-IS 478007 478969 black\nband Sv_sym_scaffold2 chromosome chromosome 478970 482518 lgrey\nband Sv_sym_scaffold2 JV46_04960-IS JV46_04960-IS 482519 483481 black\nband Sv_sym_scaffold2 chromosome chromosome 483482 487013 lgrey\nband Sv_sym_scaffold2 JV46_09350-IS JV46_09350-IS 487014 487976 black\nband Sv_sym_scaffold2 chromosome chromosome 487977 491295 lgrey\nband Sv_sym_scaffold2 JV46_09400-IS JV46_09400-IS 491296 491574 black\nband Sv_sym_scaffold2 chromosome chromosome 491575 491615 lgrey\nband Sv_sym_scaffold2 JV46_09410-IS JV46_09410-IS 491616 492422 black\nband Sv_sym_scaffold2 chromosome chromosome 492423 586162 lgrey\nband Sv_sym_scaffold2 JV46_09940-IS JV46_09940-IS 586163 587125 black\nband Sv_sym_scaffold2 chromosome chromosome 587126 597897 lgrey\nband Sv_sym_scaffold2 JV46_10000-IS JV46_10000-IS 597898 598176 black\nband Sv_sym_scaffold2 chromosome chromosome 598177 598217 lgrey\nband Sv_sym_scaffold2 JV46_10010-IS JV46_10010-IS 598218 599024 black\nband Sv_sym_scaffold2 chromosome chromosome 599025 599906 lgrey\nband Sv_sym_scaffold2 JV46_10030-IS JV46_10030-IS 599907 600869 black\nband Sv_sym_scaffold2 chromosome chromosome 600870 601260 lgrey\nband Sv_sym_scaffold2 JV46_05690-IS JV46_05690-IS 601261 602421 black\nband Sv_sym_scaffold2 chromosome chromosome 602422 637868 lgrey\nband Sv_sym_scaffold2 JV46_10250-IS JV46_10250-IS 637869 638147 black\nband Sv_sym_scaffold2 chromosome chromosome 638148 638188 lgrey\nband Sv_sym_scaffold2 JV_ISf_32-IS JV_ISf_32-IS 638189 638989 black\nband Sv_sym_scaffold2 chromosome chromosome 638990 679692 lgrey\nband Sv_sym_scaffold2 JV46_10530-IS JV46_10530-IS 679693 679971 black\nband Sv_sym_scaffold2 chromosome chromosome 679972 680012 lgrey\nband Sv_sym_scaffold2 JV46_10540-IS JV46_10540-IS 680013 680819 black\nband Sv_sym_scaffold2 chromosome chromosome 680820 778286 lgrey\nband Sv_sym_scaffold2 JV46_10880-IS JV46_10880-IS 778287 779249 black\nband Sv_sym_scaffold2 chromosome chromosome 779250 892554 lgrey\nband Sv_sym_scaffold3 chromosome chromosome 0 323 lgrey\nband Sv_sym_scaffold3 JV46_26090-ISf_39 JV46_26090-ISf_39 324 1292 black\nband Sv_sym_scaffold3 chromosome chromosome 1293 182880 lgrey\nband Sv_sym_scaffold3 JV46_24240-IS JV46_24240-IS 182881 183843 black\nband Sv_sym_scaffold3 chromosome chromosome 183844 230978 lgrey\nband Sv_sym_scaffold3 JV46_24430-IS JV46_24430-IS 230979 231941 black\nband Sv_sym_scaffold3 chromosome chromosome 231942 251012 lgrey\nband Sv_sym_scaffold3 JV46_24510-IS JV46_24510-IS 251013 251975 black\nband Sv_sym_scaffold3 chromosome chromosome 251976 355270 lgrey\nband Sv_sym_scaffold3 JV46_24950-IS JV46_24950-IS 355271 356077 black\nband Sv_sym_scaffold3 chromosome chromosome 356078 356118 lgrey\nband Sv_sym_scaffold3 JV46_24960-IS JV46_24960-IS 356119 356397 black\nband Sv_sym_scaffold3 chromosome chromosome 356398 360465 lgrey\nband Sv_sym_scaffold3 JV46_27830-ISf_37 JV46_27830-ISf_37 360466 361434 black\nband Sv_sym_scaffold3 chromosome chromosome 361435 367398 lgrey\nband Sv_sym_scaffold3 JV46_27880-IS JV46_27880-IS 367399 368544 black\nband Sv_sym_scaffold3 chromosome chromosome 368545 368680 lgrey\nband Sv_sym_scaffold3 JV46_27890-IS JV46_27890-IS 368681 369442 black\nband Sv_sym_scaffold3 chromosome chromosome 369443 369445 lgrey\nband Sv_sym_scaffold3 JV46_20730-IS JV46_20730-IS 369446 370252 black\nband Sv_sym_scaffold3 chromosome chromosome 370253 370293 lgrey\nband Sv_sym_scaffold3 JV46_20430-IS JV46_20430-IS 370294 370572 black\nband Sv_sym_scaffold3 chromosome chromosome 370573 377857 lgrey\nband Sv_sym_scaffold3 JV46_27960-IS JV46_27960-IS 377858 378820 black\nband Sv_sym_scaffold3 chromosome chromosome 378821 473344 lgrey\nband Sv_sym_scaffold3 JV46_28780-IS JV46_28780-IS 473345 473779 black\nband Sv_sym_scaffold3 chromosome chromosome 473780 474445 lgrey\nband Sv_sym_scaffold3 JV46_25380-IS JV46_25380-IS 474446 475408 black\nband Sv_sym_scaffold3 chromosome chromosome 475409 495328 lgrey\nband Sv_sym_scaffold3 JV46_25590-IS JV46_25590-IS 495329 496135 black\nband Sv_sym_scaffold3 chromosome chromosome 496136 496176 lgrey\nband Sv_sym_scaffold3 JV46_25600-IS JV46_25600-IS 496177 496455 black\nband Sv_sym_scaffold3 chromosome chromosome 496456 502617 lgrey\nband Sv_sym_scaffold3 JV46_28850-ISf_38 JV46_28850-ISf_38 502618 503586 black\nband Sv_sym_scaffold3 chromosome chromosome 503587 537612 lgrey\nband Sv_sym_scaffold4 chromosome chromosome 0 28016 lgrey\n"
}