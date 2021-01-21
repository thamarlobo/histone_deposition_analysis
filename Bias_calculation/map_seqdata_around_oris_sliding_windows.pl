#!/usr/bin/perl -w
use strict;
use warnings;
use Parallel::ForkManager;
use POSIX ();
use List::Util qw( min max );
use Math::Round;

# Script uses double-click-seq data to map and smooth the histone deposition around estimated origins.
# Script can deal with BAM files with both single end and paired end data.
# Origin information can be supplied through estimated_oris.txt (outputfile from get_significant_oris.R)
# Script can also be used to quantify OK-seq reads around estimated origins. Just supply path to OK-seq BAM file(s) instead of double-click-seq BAM file(s)
# Input: 
# 1. double-click seq BAM file(s). File(s) must be sorted, indexed and duplicates must be marked.
# 2. estimated_oris.txt (outputfile from get_significant_oris.R)
# Output: Script generates two outputfiles:
# 1. *_results.txt contains numbers of fragments mapping to the forward and reverse strand for each bin (that contains any reads) within -$region and $region around each origin. This file is used for heatmap visualization of bias for the separate origins. The information is first converted into a matrix with convert_to_matrix.pl, after which it can more easily be visualized with pheatmap (R).
# 2. *_summedresults.txt contains the summed numbers of fragments mapping to the forward and reverse strand for each bin within -$region and $region across all origins. File can be used for visualizing of the smoothed bias, but was not used for our paper.
# Outputformat:
# 1. $chr(=chromosome)  $ori_start(=origin start location)  $ori_stop(=origin_stop_location)   $ori_info{$chr}{$ori_center}{"deltafr"}(=delta_rf of origin)    $ori_info{$chr}{$ori_center}{"pvalue"}(=pvalue of origin) $ori_info{$chr}{$ori_center}{"fdr"}(=fdr of origin)   $b(=bin around origin)  $n_f(=number of fragments from forward strand in bin)   $n_r(=number of fragments from reverse strand in bin)   $frvalue(=(($n_f-$n_r)/($n_f+$n_r)). This value is eventually not used)
# 2. $b(=bin)  $n_f(=summed number of fragments from forward strand)    $n_r(=summed number of fragments from reverse strand)    $frvalue(=(($n_f-$n_r)/($n_f+$n_r)). This value is eventually not used). Columns are tab-separated.
#
# Author: Thamar Jessurun Lobo

# Paths and parameters
my @ori_files = ("PATH(S)_TO_ESTIMATED_ORIS.TXT"); # Insert path(s) to estimated_oris.txt (output from get_significant_oris.R)
my @files= ("PATH(S)_TO_DOUBLECLICKSEQ_BAMFILES"); # Insert path(s) to double-click-seq BAM file(s) or OK-seq BAM file(s)
my $outpath = "./"; # Insert path to output directory
my $min_mapq = 20; # Mimimal required mapping quality per mapped read. Default = 20
my $min_tlen = 145; # SET TO ZERO when data is not paired! Minimum required TLEN when data is paired-end. Default is 145. Should prevent use of readpairs that are too small to originate from a nucleosome.
my $region = 204000;  # Define regionsize (in basepairs) to include to each side of the origin center. Default is 204000, meaning we will look into a 200000 region to each side of the origin. The additional 4000 basepairs are needed to smooth the signal in the bins furthest away from the origin. $region can be adjusted to any region that is a multiple of $slide + ($windowsize - $slide)/2 
my $slide = 1000; # sliding step = binsize in basepairs. Default is 1000.
my $windowsize= 9000; # Sliding windowsize in basepairs. Default is 9000, in combination with $slide =1000 this means smoothing will happen with 4 bins to each side of every bin (binsize=$slide). Can be adjusted to uneven number * $slide. 

my $n_slides_to_each_side= (($windowsize/$slide)-1)/2; # Do not change! Needed for sliding window
my $n_bins = ($region / $slide) - $n_slides_to_each_side; # Do not change! Needed for sliding window.

# For each origin file:
foreach my $ori_file (@ori_files) {

# Vars
my $n_oris=0;
my %ori_info;
my $n_chroms;

# Run in parallel
#my $pm = Parallel::ForkManager->new(1); # Uncomment to run in parallel

# Extract $ori_sample and $scale from $ori_file
my $orisample = $ori_file;
$orisample =~ s/^.*\///;
$orisample = (split(/_/, $orisample))[0];
my $scale = $ori_file;
$scale =~ s/^.*bam_(.*)ndeltarfs.*$/$1/;
print "orisample is ", $orisample, " and scale is ", $scale, "\n";

# Save origin locations and info
open F, $ori_file;
warn "Saving ori info..\n";
while ( <F> ){
        chomp;
        my ( $chrom, $start, $stop, $deltafr, $pvalue, $fdr) = (split /\t/)[0,1,2,4,7,9]; 
        next if $chrom !~ m/^[0-9]{1,2}$/i; # Only autosomes
        my $center = POSIX::floor(($stop + $start) / 2); # Choose lowest if .5
        $ori_info{$chrom}{$center}{"start"}=$start;
        $ori_info{$chrom}{$center}{"stop"}=$stop;
        $ori_info{$chrom}{$center}{"deltafr"}=$deltafr; 
        $ori_info{$chrom}{$center}{"fdr"}=$fdr;
        $ori_info{$chrom}{$center}{"pvalue"}=$pvalue;
        $n_oris++;
}
close F;
$n_chroms = scalar keys %ori_info;
warn "Saved info on $n_oris oris on $n_chroms chromosomes ..\n";

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
                next if $flag & 1024; # Ignore marked duplicates
                next if $flag & 256; # Ignore secondary alignments (should not be reported by bowtie2)
                next if $flag&1 and ( $flag&16 or not($flag&2) ); # Paired End (PE), do not consider minus strand reads or discordantly mapped
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
                        $diad_pos = $diad_pos -75; # Define $diad_pos as the middle basepair of the fragment the read originates from. We assume a fragment size of 150 basepairs.
                        push @{$diads{$chr}{"pos"}}, $diad_pos;
                        push @{$diads{$chr}{"strand"}}, "r";
                    } else{ # Forward strand
                        $diad_pos = $pos + 75; # Define $diad_pos as the middle basepair of the fragment the read originates from. We assume a fragment size of 150 basepairs.
                        push @{$diads{$chr}{"pos"}}, $diad_pos ;
                        push @{$diads{$chr}{"strand"}}, "f";
                    }
                }
            }
            close F;  
            
            warn "Quantifying histones from $file around origins : $chr \n";
            
            # Open outputfile
            open F1, '>'.$outpath.$filen."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_slidesize".$slide."_mintlen".$min_tlen."_chr".$chr."_results.txt";
            
            # Get diad centers and strandinfo sorted by ascending position
            my @ds = @{ $diads{$chr}{"pos"} };
            my @strs = @{ $diads{$chr}{"strand"}};

            # Get ori centers sorted by ascending left end position
            my @ori_centers = sort {$a <=> $b} keys %{$ori_info{$chr}};
            
            # For each origin center
            foreach my $ori_center( @ori_centers ){
            
                # Pick up further origin location info
                my $ori_stop = $ori_info{$chr}{$ori_center}{"stop"};
                my $ori_start = $ori_info{$chr}{$ori_center}{"start"};
				
				while (@ds and $ds[0] < ($ori_center - ($region + 300))){  # Speed up runtime when diad center < left region ori center (incl 300 bp buffer)
					shift @ds;
					shift @strs;
                }
                
				# Create hash for saving origin info
                my %ori_diad_counts=();
                
				# Define vars
                my $position;
                my $bin_center;
                my $bin_start;
                my $bin_stop;
                  
                # Count how many forward-strand and reverse-strand diad centers map to each bin (size = $slide) to each side of the origin (limited by $region), and add these numbers of diad centers up to n=$n_slides_to_each_side bins to each side of the original bin as well.
                foreach my $ele (0 .. $#ds){
					last if $ds[$ele] > ($ori_center + $region + 300); # Exit loop if current diad_center is larger than right region + 300 bp buffer
                    if (($ds[$ele] >= ($ori_center - $region) && $ds[$ele] <= ($ori_center + $region))){  # Diad_center lies within -region and region around origin center
						$position = $ds[$ele] - $ori_center;
                        $bin_center = round($position/$slide);
                        $bin_start = $bin_center - $n_slides_to_each_side;
                        $bin_stop = $bin_center + $n_slides_to_each_side;
                        if (abs($bin_start) > $n_bins){ 
							$bin_start= -1 * $n_bins;
                        }
                        if (abs($bin_stop) > $n_bins){ 
							$bin_stop= $n_bins;
                        }
                   }else{
                       next; # Diad_center is located within buffer region
                   }
                   my @bins = ($bin_start, $bin_stop);
                   my $min_bin= min @bins;
                   my $max_bin = max @bins;
                   if($strs[$ele] eq "f"){
						foreach my $b ($min_bin .. $max_bin){
                            $ori_diad_counts{$b}{"n_f"} +=1;
                        }
                   }elsif($strs[$ele] eq "r"){
                        foreach my $b ($min_bin .. $max_bin){
                            $ori_diad_counts{$b}{"n_r"} +=1; #Increment the correct counter
                        }
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
                    print F1 join("\t", $chr, $ori_start, $ori_stop, $ori_info{$chr}{$ori_center}{"deltafr"},  $ori_info{$chr}{$ori_center}{"pvalue"} ,  $ori_info{$chr}{$ori_center}{"fdr"}, $b, $n_f, $n_r, $frvalue), "\n";
                }
            }
            close F1;
			#$pm-> finish; # Uncomment to run in parallel
		}
        #$pm->wait_all_children; # Uncomment to run in parallel

        # Concatenate separate chromosomes outputfiles
		my $outpattern = $outpath.$filen."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_slidesize".$slide."_mintlen".$min_tlen."_chr*_results.txt";
		my $outfn =  $outpath.$filen."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_slidesize".$slide."_mintlen".$min_tlen."_results.txt";
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
			my ($chr, $ori_start, $ori_stop,  $deltafr, $pvalue, $fdr, $b, $n_f, $n_r, $frvalue) = split(/\t/);
		    $summed_diad_counts{$b}{"n_f"} += $n_f;
		    $summed_diad_counts{$b}{"n_r"} += $n_r;
		}
		close F2;
		
        # Print fr values for summed dyads to F3 outputfile       
        warn "Writing summed diad counts to outputfile\n";
        open F3, ">".$outpath.$filen."_".$orisample."_scale".$scale."_region".$region."_windowsize".$windowsize."_slidesize".$slide."_mintlen".$min_tlen."_summedresults.txt";
        foreach my $b (sort {$a <=> $b} keys %summed_diad_counts){
            my $n_f = 0;
            my $n_r = 0;
            $n_f += $summed_diad_counts{$b}{"n_f"} if exists $summed_diad_counts{$b}{"n_f"};
            $n_r += $summed_diad_counts{$b}{"n_r"} if exists $summed_diad_counts{$b}{"n_r"};
            
            # Calculate f-r/f+r and set $frvalue to NA if both n_f and n_r are 0
            my $frvalue;
            if($n_f + $n_r == 0){
				$frvalue = "NA";
            }else{
                $frvalue = ($n_f-$n_r)/($n_f+$n_r);
            }
            
            # Write summed results to "summedresults.txt" outputfile
            print F3 join("\t", $b, $n_f, $n_r, $frvalue), "\n";
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

                
        
