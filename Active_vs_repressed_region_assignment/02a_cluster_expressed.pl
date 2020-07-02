#!/usr/bin/perl -w
use strict;

# Define paths and parameters
my $inputfile = "PATH_TO_EXPRESSIONFILE"; # insert path to SRR4787042Aligned.sortedByCoord.out.bam (obtained by mapping hRPE-1 gene expression as downloaded from GEO accession GSE89413, with STAR aligner)
my $annotation = "PATH_TO_ANNOTATION"; # insert path to Homo_sapiens.GRCh38.98.gtf
my $inputfile2 = "PATH_TO_READSPERGENE"; # insert path to SRR4787042ReadsPerGene.out.tab (obtained by mapping hRPE-1 gene expression as downloaded from GEO accession GSE89413, with STAR aligner)
my $min_fpkm = 1;
my $max_gene_dist = 10_000;

# Load chromosome lengths from BAM file header
my %chrlen = ();
open F,'samtools view -H '.$inputfile.' |';
while ( <F> ) {
    next unless m/^\@SQ.*SN\:(\S+).*LN\:(\d+)/; # sequence header in BAM starts with @SQ
    my ( $chr, $len ) = ( $1, $2 );
    $chrlen{$chr} = $len;
    warn $chr, ' -> ', $len, "\n";
}
close F;

warn "Loading gene positions\n";
my %gene_coord = ();
my $i = 0;
open F, $annotation;
while ( <F> ) {
    next if m/^\#/; # skip header
    chomp;
    my ( $chr, $source, $feature, $start, $end, $score, $strand, $phase, $info ) = split /\t/;
    if ( $feature eq 'gene' ) { # process gene information
        my ( $gene_id ) = $info =~ m/gene_id \"([^\"]*)\"/;
        $gene_coord{$gene_id}{'chr'} = $chr;
        $gene_coord{$gene_id}{'start'} = $start;
        $gene_coord{$gene_id}{'end'} = $end;
        warn $i, " genes done\n" unless ++$i % 1000;
    }
    elsif ( $feature eq 'exon' ) { # process exone information
        my ( $gene_id ) = $info =~ m/gene_id \"([^\"]*)\"/;
        foreach my $pos ( $start .. $end ) { # for each gene store each exonic coordinate in hash. The total number of unique elements would be total exon length
            $gene_coord{$gene_id}{$pos}++;
        }
    }
}
close F;

my %fpkm = (); # hash to store FPKMs
my $reads = 0; # counter for total reads
warn "Reading expression values\n";
open F, $inputfile2;
while ( <F> ) {
    my ( $gene, $unstranded, $first, $second ) = split /\t/; # STAR outputs all 3 possibilites of library prep method 
    next if $gene =~ m/^N\_/; # skip special features like "nofeature", "ambigous" and "unmapped" 
    # The total number of exonic bases per gene will be equal to total number of elements in $gene_coord{$gene_id} sub-hash
    # We need to subtract 3 elements though (chr, start, end) that we also store in that hash
    # and divide it by 1000 (that is multiply by 0.001 so that we have it per kilobase rather than per base pair
    my $exons_kb = 0.001 * ( scalar( keys %{$gene_coord{$gene}} ) - 3); 
    $fpkm{$gene} = $second / $exons_kb; #Now fpkm has read per kilobase, but we store total reads in reads and will normalize it "per million" later
    $reads += $second;
}
close F;

# Apply normalization "per milliom reads", save active genes
my %genes = (); # to store chr, starts and ends of active genes (where expression > min_fpkm)
my $expressed = 0; # stores total number of expressed genes
foreach my $gene ( keys %fpkm ) {
    $fpkm{$gene} *= 1_000_000 / $reads; # what is stored in $fpkm already (reads/kb) is multiplied by a million and divided by total reads, so we are finally getting FPKMs :)
    if ( $fpkm{$gene} >= $min_fpkm ) { 
        $expressed++;
        $genes{$gene_coord{$gene}{'chr'}}{$gene_coord{$gene}{'start'}} = $gene_coord{$gene}{'end'};
    }
}
warn "expressed: ", $expressed, " out of ", scalar keys %fpkm, "\n";

MAIN:    # there are nested loops, so I call the main loop "MAIN"
while ( 1 == 1 ) { # Exit criteria is absence of merges between regions, checked in the end of the loop, not here
    my $changes = 0; 
    foreach my $chr ( keys %genes ) {
        my @starts = sort {$a<=>$b} keys %{$genes{$chr}}; # Sort genes by start, 1 chromosome at a time
        warn $chr, ' : ', scalar(@starts), "\n";
        foreach my $ele ( 0 .. $#starts-1 ) { # for each element of array check if we should not merge it with the following
            if ( $starts[$ele+1] <= $genes{$chr}{$starts[$ele]} + $max_gene_dist ) { # if they are close:
                my @arr = sort {$a<=>$b} (@starts[$ele, $ele+1], $genes{$chr}{$starts[$ele]}, $genes{$chr}{$starts[$ele+1]} ); # sort start and end positions of both genes. First element will be start and last the end of the merger
                delete $genes{$chr}{$starts[$ele+1]}; # remove second gene
                $genes{$chr}{$starts[$ele]} = $arr[-1]; # set the end of the gene to the biggest coordinate of the four coordinates
                $changes++; # to indicate the change is made
                next MAIN; # now we return to the main loop as we've changed %genes hash and need to re-initiate @starts, etc
            }
        }
    }
    last unless $changes; # exit only if cannot make more changes
}

warn "Saving\n";
my $total_bp = 0; # Stores total number of bases in active regions
open F, '>', '20191126_RPE1_active_passive.bed';
foreach my $chr ( sort keys %genes ) {
    my @coords = ( 0 ); # now hash stores 1st base of chr, it will be start of passive domain
    foreach my $start ( sort {$a<=>$b} keys %{$genes{$chr}} ) { # sort gene starts numerically
        $total_bp += $genes{$chr}{$start} - $start + 1; # calculate length of active region, add it to total_bp
        push @coords, $start-1; # BED start is zero-based, substract 1
        push @coords, $genes{$chr}{$start}; # this value indicates end of gene and is 1-based
    }
    push @coords, $chrlen{$chr}; # last element is end base of chromosome, also passive domain
    warn "Outputting $chr\n";
    foreach my $ele ( 0 .. $#coords-1 ) {
        my $type = $ele % 2 ? 'ACTIVE' : 'PASSIVE'; # active and passive domains are alterating in our array (even are passive, odd are active)
        print F join( "\t", $chr, $coords[$ele], $coords[$ele+1], $type), "\n" if $coords[$ele+1] > $coords[$ele];
    }
}
close F;
warn "Total bp in active regions: ", $total_bp, "\n";

