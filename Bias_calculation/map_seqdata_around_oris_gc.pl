#!/usr/bin/perl -w
use strict;
use warnings;
use Parallel::ForkManager;
use POSIX ();
use Math::Round;


# Script uses double-click-seq data to map histone deposition around estimated origins that are grouped into quantiles based on GC-percentage OR RT-score.
# The input BAM files must be sorted and duplicates must be marked.
# Script can deal with BAM files with both single end and paired end data.
# Origin information can be supplied through estimated_oris.txt (outputfile from get_significant_oris.R) where a last column (column 11) containing either the the GC-percentage OR the RT-score of the region around the origin MUST be added!
# Script can also be used to quantify OK-seq reads around estimated origins. Just supply path to OK-seq BAM file(s) instead of double-click-seq BAM file(s)
# Input: 
# 1. double-click seq BAM file(s). File(s) must be sorted, indexed and duplicates must be marked.
# 2. estimated_oris.txt (outputfile from get_significant_oris.R) where a last column (column 11) containing the GC-percentage OR the RT-score of the region around the origin MUST be added.
# Output: Script generates two outputfiles:
# 1. *_results.txt contains numbers of fragments mapping to the forward and reverse strand for each bin (that contains any reads) within -$region and $region around each origin. This file is used for bootstrap analysis with bootstrap_bias_gcquantiles.pl
# 2. *_summedresults.txt contains the summed numbers of fragments mapping to the forward and reverse strand for each bin within -$region and $region across all origins (grouped into qunatiles based on their GC percentage OR RT-score). This file is used for visualizing the average bias across origins.
# Outputformat:
# 1. $chr(=chromosome)  $ori_start(=origin start location)  $ori_stop(=origin_stop_location)   $ori_info{$chr}{$ori_center}{"gc"}(= GC percentage of origin region OR RT-score )    $quantile(= GC quantile OR RT quantile)   $ori_info{$chrom}{$center}{"fdr"}(=fdr of origin)   $b(=bin around origin)  $n_f(=number of fragments from forward strand in bin)   $n_r(=number of fragments from reverse strand in bin)   $frvalue(=(($n_f-$n_r)/($n_f+$n_r)). This value is eventually not used)
# 2. $b(=bin)  $n_f(=summed number of fragments from forward strand)    $n_r(=summed number of fragments from reverse strand)    $frvalue(=(($n_f-$n_r)/($n_f+$n_r)). This value is eventually not used)    $q(=GC quantile OR RT quantile). Columns are tab-separated.
#
# Author: Thamar Jessurun Lobo

# Paths and parameters
my @ori_files = ("PATH(S)_TO_ESTIMATED_ORIS.TXT"); # Insert path(s) to estimated_oris.txt (output from get_significant_oris.R) with an additional 11th column containing GC information OR RT score
my @files= ("PATH(S)_TO_DOUBLECLICKSEQ_BAMFILES"); # Insert path(s) to double-click-seq BAM file(s) or OK-seq BAM file(s)
my $outpath = "./"; # Insert path to output directory
my $min_mapq = 20; # Mimimal required mapping quality per mapped read. Default = 20
my $region = 200000;  # Define regionsize (in basepairs) to include to each side of the origin center. Default is 200000. MUST be a multiple of $windowsize!
my $min_tlen = 145; # SET TO ZERO when data is not paired! Minimum required TLEN when data is paired-end. Default is 145. Should prevent use of readpairs that are too small to originate from a nucleosome.
my $windowsize= 4000; # Windowsize in basepairs. Default is 4000.
my $type = "GC"; # GC or RT

# For each origin file:
foreach my $ori_file (@ori_files) {

# Vars
my $n_oris=0; 
my @gcs=(); 
my %ori_info;
my $n_chroms;

# Run in parallel
#my $pm = Parallel::ForkManager->new(1); # Uncomment to run in parallel

# Extract $ori_sample and $scale from $ori_file
my $orisample = $ori_file;$orisample =~ s/^.*\///;
$orisample = (split(/_/, $orisample))[0];
my $scale = $ori_file;
$scale =~ s/^.*bam_(.*)ndeltarfs.*$/$1/;
print "orisample is ", $orisample, " and scale is ", $scale, "\n";

# Save origin locations and info
open F, $ori_file;
warn "Saving ori info..\n";
while ( <F> ){
    chomp;
    my ( $chrom, $start, $stop, $deltafr, $pvalue, $fdr, $gc) = (split /\t/)[0,1,2,4,7,9,10];
    next if $chrom !~ m/^[0-9]{1,2}$/i; # Only autosomes
    my $center = POSIX::floor(($stop + $start) / 2); # Choose lowest if .5
    $ori_info{$chrom}{$center}{"start"}=$start;
    $ori_info{$chrom}{$center}{"stop"}=$stop;
    $ori_info{$chrom}{$center}{"gc"}=$gc; 
    $ori_info{$chrom}{$center}{"fdr"}=$fdr;
    $ori_info{$chrom}{$center}{"pvalue"}=$pvalue;
    push @gcs, $gc; 
    $n_oris++;
}
close F;
$n_chroms = scalar keys %ori_info;
warn "Saved info on $n_oris oris on $n_chroms chromosomes ..\n";

# Determine quantiles GC percentage 
my @sorted_gcs = sort {$a <=> $b } @gcs; # Sort gcs 
my $Q1 = @sorted_gcs[@sorted_gcs * 0.25]; 
my $Q2 = @sorted_gcs[@sorted_gcs * 0.5]; 
my $Q3 = @sorted_gcs[@sorted_gcs * 0.75]; 

# For each file:
foreach my $file(@files) {
    my %diads = ();
    warn "Working on $file ..\n";
    
    # Get filename from file path
    my $filen = $file;
    $filen =~ s/^.*\///g;
 
    if (-e $file){ # Dummy check
    
    # For each autosome:
    foreach my $chr (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22){

        # Run in parallel
        #$pm->start and next;  # Uncomment to run in parallel
            
        warn "Saving histone info from $file : $chr..\n";
        
        # Access alignments with a minimum mapping quality of $min_mapq on chromosome $chr
        open F, 'samtools view -q '.$min_mapq.' '.$file.' '.$chr.'|'; 
    
    	# Save diad centre locations and strands
    	while ( <F> ) {
    		chomp;
			my ( $read, $flag, $chr, $pos, $cigar, $tlen ) = (split /\t/)[0,1,2,3,5,8];
			next if $flag & 1024;  # Ignore marked duplicates
			next if $flag & 256; # Ignore secondary alignments (should not be reported by bowtie2)
			next if $flag&1 and ( $flag&16 or not($flag&2) );  # Paired End (PE), do not consider minus strand reads or discordantly mapped
			next if abs($tlen) < $min_tlen; # Ignore aligments where TLEN is smaller than $min_tlen
			my $diad_pos;
			if ($flag&1){ # If Paired End (PE) data      
				$diad_pos=$pos + POSIX::floor(abs($tlen)/2); # Define $diad_pos as the middle basepair of aligned readpair             
				if ($flag&64){                        
					# Read is first read, so fragment comes from plus strand
					push @{$diads{$chr}{"pos"}}, $diad_pos;
					push @{$diads{$chr}{"strand"}}, "f";
				}else{
					# Read is second read, so fragment comes from minus strand
					push @{$diads{$chr}{"pos"}}, $diad_pos;
					push @{$diads{$chr}{"strand"}}, "r";
				}
			}else{ # If Single End data:
				if ($flag & 16) { # Reverse strand
					$diad_pos = $pos -1;
					while ( $cigar =~ m/(\d+)[MD]/g ) { $diad_pos +=$1 };
					$diad_pos = $diad_pos -75; # Define $diad_pos as the middle basepair of the fragment the read originates from. We assume a fragment size of 150 basepairs
					push @{$diads{$chr}{"pos"}}, $diad_pos;
					push @{$diads{$chr}{"strand"}}, "r";
				} else{ # Forward strand
					$diad_pos = $pos + 75; # Define $diad_pos as the middle basepair of the fragment the read originates from. We assume a fragment size of 150 basepairs
					push @{$diads{$chr}{"pos"}}, $diad_pos ;
					push @{$diads{$chr}{"strand"}}, "f";
				}
			}
		}
		close F;

        warn "Quantifying histones from $file around origins : $chr \n";
        
        # Open outputfile
        open F1, '>'.$outpath.$filen."_".$type."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_mintlen".$min_tlen."_chr".$chr."_results.txt";

        # Get diad centers and strandinfo sorted by ascending position
        my @ds = @{ $diads{$chr}{"pos"} };
        my @strs = @{ $diads{$chr}{"strand"}};
			
        # Get origin centers sorted by ascending left end position
        my @ori_centers = sort {$a <=> $b} keys %{$ori_info{$chr}};
        
        # For each origin center
        foreach my $ori_center ( @ori_centers ){
			
            # Determine quantile of its GC percentage  
            my $gc = $ori_info{$chr}{$ori_center}{"gc"}; 
            my $quantile; 
            if ($gc > $Q2) {
                if($gc > $Q3){ 
                    $quantile = "Q4"; 
                }else{ 
                    $quantile = "Q3"; 
                } 
            }else{ 
                if($gc > $Q1){ 
                    $quantile = "Q2"; 
                }else{ 
                    $quantile = "Q1"; 
                } 
            }  
				
            # Pick up further origin location info
            my $ori_stop = $ori_info{$chr}{$ori_center}{"stop"};
            my $ori_start = $ori_info{$chr}{$ori_center}{"start"};
                
            while (@ds and $ds[0] < ($ori_center - ($region + 300))){ # Speed up runtime when diad center < left region ori center (incl 300 bp buffer)
                shift @ds;
                shift @strs;
            }
                
            # Create hash for saving origin info
            my %ori_diad_counts=();
				
            # Define vars
            my $position;
            my $bin_center;

            # Count how many forward-strand and reverse-strand diad centers map to each bin (size = $windowsize) to each side of the origin (limited by $region).
            foreach my $ele (0 .. $#ds){
                last if $ds[$ele] > ($ori_center + $region + 300); # Exit loop if current diad_center is larger than right $region + 300 bp buffer
                if (($ds[$ele] >= ($ori_center - $region) && $ds[$ele] <= ($ori_center + $region))){   # Diad_center lies within -region and region around origin center
                    $position = $ds[$ele] - $ori_center;
                    $bin_center = round($position/$windowsize);
                }else{
                    next; # Diad_center is located within buffer region
                }
                if($strs[$ele] eq "f"){
			$ori_diad_counts{$bin_center}{"n_f"} +=1; 
                }elsif($strs[$ele] eq "r"){
                        $ori_diad_counts{$bin_center}{"n_r"} +=1; # Increment the correct counter
                }
            }
            # For each bin:
            foreach my $b( sort {$a <=> $b} keys %ori_diad_counts){
                my $n_f = 0;
                my $n_r = 0;
                $n_f += $ori_diad_counts{$b}{"n_f"} if exists $ori_diad_counts{$b}{"n_f"};
                $n_r += $ori_diad_counts{$b}{"n_r"} if exists $ori_diad_counts{$b}{"n_r"};
                        
                # Calculate f-r/f+r and set $frvalue to NA if both n_f and n_r are 0
                my $frvalue;
                if($n_f + $n_r == 0){
                    $frvalue = "NA";
                }else{                           
                    $frvalue = (($n_f-$n_r)/($n_f+$n_r));                        
                }
                # Write result to outputfile
                print F1 join("\t", $chr, $ori_start, $ori_stop, $ori_info{$chr}{$ori_center}{"gc"}, $quantile, $ori_info{$chr}{$ori_center}{"fdr"}, $b, $n_f, $n_r, $frvalue), "\n"; 
            }
        }
        close F1;
       #$pm-> finish; # Uncomment to run in parallel
    }
    #$pm->wait_all_children; # Uncomment to run in parallel

    # Concatenate separate chromosomes outputfiles
    my $outpattern = $outpath.$filen."_".$type."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_mintlen".$min_tlen."_chr*_results.txt";
    my $outfn =  $outpath.$filen."_".$type."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_mintlen".$min_tlen."_results.txt";
    system join(' ',
            'cat',
            $outpattern, 
            '>',
            $outfn,
    );
            
    # Remove separate chromosomes outputfiles
    system join(' ',
            'rm',
            $outpattern
    );
		
    # Get summed diad counts per bin across origins
    open F2, $outfn;
    my %summed_diad_counts;
    while (<F2>){
        chomp;
        my ($chr, $ori_start, $ori_stop,  $gc, $gc_quantile, $fdr, $b, $n_f, $n_r, $frvalue) = split(/\t/);
        $summed_diad_counts{$gc_quantile}{$b}{"n_f"} += $n_f;
        $summed_diad_counts{$gc_quantile}{$b}{"n_r"} += $n_r;
    }
    close F2;
		
    # Print fr values for summed dyads to F3 outputfile   
    warn "Writing summed diad counts to outputfile\n";
    open F3, ">".$outpath.$filen."_".$type."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_mintlen".$min_tlen."_summedresults.txt";
    foreach my $q (keys %summed_diad_counts){ #gc
        foreach my $b (sort {$a <=> $b} keys %{$summed_diad_counts{$q}}){ #$q = gc
            my $n_f = 0;
            my $n_r = 0;
            $n_f += $summed_diad_counts{$q}{$b}{"n_f"} if exists $summed_diad_counts{$q}{$b}{"n_f"};
            $n_r += $summed_diad_counts{$q}{$b}{"n_r"} if exists $summed_diad_counts{$q}{$b}{"n_r"};
            
            # Calculate f-r/f+r and make it NA if both n_f and n_r are 0
            my $frvalue;
            if($n_f + $n_r == 0){
                $frvalue = "NA";
            }else{
                $frvalue = ($n_f-$n_r)/($n_f+$n_r);
            }
            
            # Write summed results to "summedresults.txt" outputfile
            print F3 join("\t", $b, $n_f, $n_r, $frvalue, $q), "\n";
        }
    }
    close F3;    
    }else{	
        warn "file does not exist: $file \n";
        next;
    }
}
warn "Done with histones around origins from $ori_file..";
}
warn "Done!\n";


