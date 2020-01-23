#!/usr/bin/perl -w
use strict;

my $min_fpkm = 1;
my $max_gene_dist = 10_000;

my %chrlen = ();
open F, 'samtools view -H SRR4787042Aligned.sortedByCoord.out.bam |';
while ( <F> ) {
    next unless m/^\@SQ.*SN\:(\S+).*LN\:(\d+)/;
    my ( $chr, $len ) = ( $1, $2 );
    $chrlen{$chr} = $len;
    warn $chr, ' -> ', $len, "\n";
}
close F;

warn "Loading gene positions\n";
my %gene_coord = ();
my $i = 0;
open F, 'Homo_sapiens.GRCh38.98.gtf';
while ( <F> ) {
    next if m/^\#/;
    chomp;
    my ( $chr, $source, $feature, $start, $end, $score, $strand, $phase, $info ) = split /\t/;
    if ( $feature eq 'gene' ) {
        my ( $gene_id ) = $info =~ m/gene_id \"([^\"]*)\"/;
        $gene_coord{$gene_id}{'chr'} = $chr;
        $gene_coord{$gene_id}{'start'} = $start;
        $gene_coord{$gene_id}{'end'} = $end;
        warn $i, " genes done\n" unless ++$i % 1000;
    }
    elsif ( $feature eq 'exon' ) {
        my ( $gene_id ) = $info =~ m/gene_id \"([^\"]*)\"/;
        foreach my $pos ( $start .. $end ) {
            $gene_coord{$gene_id}{$pos}++;
        }
    }
}
close F;

my %fpkm = ();
my $reads = 0;
warn "Reading expression values\n";
open F, 'SRR4787042ReadsPerGene.out.tab';
while ( <F> ) {
    my ( $gene, $unstranded, $first, $second ) = split /\t/;
    next if $gene =~ m/^N\_/;
    my $exons_kb = 0.001 * ( scalar( keys %{$gene_coord{$gene}} ) - 3);
    $fpkm{$gene} = $second / $exons_kb;
    $reads += $second;
}
close F;

my %genes = ();
my $expressed = 0;
foreach my $gene ( keys %fpkm ) {
    $fpkm{$gene} *= 1_000_000 / $reads;
    if ( $fpkm{$gene} >= $min_fpkm ) {
        $expressed++;
        $genes{$gene_coord{$gene}{'chr'}}{$gene_coord{$gene}{'start'}} = $gene_coord{$gene}{'end'};  
    }
}
warn "expressed: ", $expressed, " out of ", scalar keys %fpkm, "\n";

MAIN:
while ( 1 == 1 ) {
    my $changes = 0;
    foreach my $chr ( keys %genes ) {
        my @starts = sort {$a<=>$b} keys %{$genes{$chr}};
        warn $chr, ' : ', scalar(@starts), "\n";
        foreach my $ele ( 0 .. $#starts-1 ) {
            if ( $starts[$ele+1] <= $genes{$chr}{$starts[$ele]} + $max_gene_dist ) {
                my @arr = sort {$a<=>$b} (@starts[$ele, $ele+1], $genes{$chr}{$starts[$ele]}, $genes{$chr}{$starts[$ele+1]} );
                delete $genes{$chr}{$starts[$ele+1]};
                $genes{$chr}{$starts[$ele]} = $arr[-1];
                $changes++;
                next MAIN;
            }
        }
    }
    last unless $changes; 
}

warn "Saving\n";
my $total_bp = 0;
open F, '>', '20191126_RPE1_active_passive.bed';
foreach my $chr ( sort keys %genes ) {
    my @coords = ( 0 );
    foreach my $start ( sort {$a<=>$b} keys %{$genes{$chr}} ) {
        $total_bp += $genes{$chr}{$start} - $start + 1;
        push @coords, $start-1;
        push @coords, $genes{$chr}{$start};
    }
    push @coords, $chrlen{$chr};
    warn "Outputting $chr\n";
    foreach my $ele ( 0 .. $#coords-1 ) {
        my $type = $ele % 2 ? 'ACTIVE' : 'PASSIVE';
        print F join( "\t", $chr, $coords[$ele], $coords[$ele+1], $type), "\n" if $coords[$ele+1] > $coords[$ele];
    }
}
close F;
warn "Total bp in active regions: ", $total_bp, "\n";


