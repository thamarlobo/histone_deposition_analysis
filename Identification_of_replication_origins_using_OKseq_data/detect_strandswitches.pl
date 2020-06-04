#!/usr/bin/perl

# Script implements deltaWC algorithm to calculate strand-switching regions using bam files coming from OK-seq.
# Optionally, script can randomly shuffle strand information and call background strand-switches (to give an indication of noise effect). Random shuffling is activated by setting $randomize to "r".
# Outputfile .*deltarfs.txt should be run through get_significant_oris.R to get a list of origins.
# Input: OK-seq BAM file(s). File(s) MUST be sorted on position and duplicate reads MUST be marked!
# Output: 1. *deltarfs.txt contains strandswitches. 2. *deltarf_stats.txt contains the number and percenatge of calls for each delta_fr_ratio 
# Outputformat 1: chromosome    strandswitch_start  strandswitch_stop   delta_fr_ratio(=delta_fr/$scale)   delta_fr(=($n_r_right + $n_f_left) - $scale)   strandswitch_left_f(=number of reads from forward strand to the left of the strandswitch) strandswitch_right_r(=number of reads from reverse strand to the right of the strandswitch). Columns are tab-separated.
# Outputformat 2:  delta_fr_ratio   $delta_stats{$delta_fr_ratio}(=exact occurence) 100*$delta_stats{$delta_fr_ratio}/$i)(=occurence in percentages). Columns are tab-separated.
#
# Author: Thamar Jessurun Lobo

use warnings;
use strict;
#use List::Util qw(min max);
use 5.010;
#use List::MoreUtils qw(first_index);
#use Data::Dumper;

# Paths and Parameters 
my @bams = ("INSERT_PATH(S)_TO_YOUR_BAM_FILES");
my $min_mapq = 20; # Mimimal required mapping quality per mapped read. Default = 20
my @scales = (100);# Number of reads you want to consider on each side of a strand switch. Default is 100.
my $randomize = "n"; # Randomize strand information? yes = "r" or no = "n"
my $fragment_size=150; # Expected size of Okazaki fragments in basepairs. Default = 150

# For each BAM file:
foreach my $bam ( @bams) {
    my %valid_aligned_fragments;
    my %valid_aligned_fragments_info;
    my %delta_stats;
    my %delta_regions;
    my %valid_read_info;
    my %deltas;
    my @left;
    my @right;
    my $i;
    my $current_valid_read;
    my $shifted_off;
    my $n_f_left;
    my $n_r_right;
    my $delta_fr;
    
    # Open BAM file: only show alignments with a minimum mapping quality of $min_mapq
    open F2, 'samtools view -q '.$min_mapq.' '.$bam.'|'; 
    my $valid_bam_read = 0;
    warn "Determining and saving Okazaki fragment location for reads in $bam \n";
    
    # For each alignment:
    while ( <F2> ){
        my ( $read, $flag, $chr, $pos, $cigar, $seq ) = (split /\t/)[0,1,2,3,5,9];
        next if $chr !~ m/^[0-9]{1,2}$|^MT$|^Y$|^X$/i;
        next if $flag & 1024; # Ignore marked duplicates
        next if $flag & 256; # Ignore secondary alignments (should not be reported by bowtie2)
        my $strand;
        $strand =  $flag & 16 ? "-" : "+";
        
        $valid_bam_read += 1; # Current valid bam read number
        
        # Get right-end location of read
        my $stop = $pos -1;
        while ( $cigar =~ m/(\d+)[MD]/g ) { $stop +=$1 }
        
        # Determine right-end location of Okazaki fragment that reads must originate from
        my $fragment_stop = $strand eq "-"? $stop : $pos + ($fragment_size -1);
        
        # Determine left-end location of Okazaki fragment that reads must originate from
        my $fragment_start =$strand eq "-"? $stop - ($fragment_size - 1): $pos;
        
        # Randomize strand information if $randomize= "r"
	if ($randomize eq "r"){
            my $random_n = rand();
            $strand = $random_n < 0.5 ? "-" : "+";
	}
	
        # Save info on Okazaki fragment
        $valid_aligned_fragments{$chr}{$valid_bam_read}=$fragment_start;
        $valid_aligned_fragments_info{$valid_bam_read}{"start"}= $fragment_start;
        $valid_aligned_fragments_info{$valid_bam_read}{"stop"}= $fragment_stop;
        $valid_aligned_fragments_info{$valid_bam_read}{"chr"}= $chr;
        $valid_aligned_fragments_info{$valid_bam_read}{"strand"}= $strand;
    }
    close F2;
    warn "Saved fragment information on $valid_bam_read valid reads\n";
    
    # For each $scale:
    foreach my $scale (@scales) {
		warn "Scale = $scale \n";
		my %delta_stats;
		my %delta_regions;
		my %valid_read_info;
		my %deltas;
		my @left;
		my @right;
		my $i;
		my $current_valid_read;
		my $shifted_off;
		my $n_f_left;
		my $n_r_right;
		my $delta_fr;
    
		open F, ">".$bam.'_'.$scale.$randomize.'deltarfs.txt'; # Open outputfile
		
		# For each chromosome:
		foreach my $chromosome (sort keys %valid_aligned_fragments){
			
			# Reset variables for new chromosome
			%delta_regions = ();
			%valid_read_info = ();
			%deltas = ();
			@left = ();
			@right = ();

			# Set variables
			my %fragment_to_bam_read = ();
			my $current_valid_fragment = 0;
			
			warn "Calculating deltaFRs for $bam : chromosome $chromosome \n";
			
			# For each valid bam read, sorted by fragment start position:
			foreach my $current_valid_bam_read (sort { $valid_aligned_fragments{$chromosome}{$a} <=> $valid_aligned_fragments{$chromosome}{$b} } keys %{$valid_aligned_fragments{$chromosome}}) {
			
				$current_valid_fragment += 1 ;
				$fragment_to_bam_read{$current_valid_fragment}=$current_valid_bam_read;
				
				# Retrieve strand information for read
				my $strand = $valid_aligned_fragments_info{$current_valid_bam_read}{"strand"}; 
				
				# Push strand information into right array
				$strand eq "-" ? push (@right, "r"): push (@right, "f") ; 
				next if scalar @right < ($scale + 1);
				
				# Move the first value from the right array to the end of the left array (shift and push) when right array is full
				$shifted_off = shift @right; 
				push (@left, $shifted_off);
				next if scalar @left < ($scale + 1);
				
				# Move first value off left array (shift) once the left array is full
				shift @left;
				
				# Calculate deltafr
				$n_f_left = scalar(grep { $_ eq "f" } @left);
				$n_r_right = scalar(grep { $_ eq "r" } @right); 
				$delta_fr = ($n_r_right + $n_f_left) - $scale; # Maximum value of deltafr will be equal to $scale. Minimum value will be equal to -($scale). Random (no switch) will be 0. 
				
				# Save deltafrin hash with key=current valid fragment
				$deltas{$current_valid_fragment}=$delta_fr;
				
				# Save the strand switch region for current_valid_fragment
				$delta_regions{$current_valid_fragment}{"start"}= $valid_aligned_fragments_info{$fragment_to_bam_read{$current_valid_fragment - $scale}}{"stop"} + 1;
				$delta_regions{$current_valid_fragment}{"end"}= $valid_aligned_fragments_info{$fragment_to_bam_read{$current_valid_fragment - ($scale-1)} }{"start"} -1;
				
				# Save the number of reads mapping to the forward strand from the left array and the number of reads mapping to the reverse strand from the right array
				$delta_regions{$current_valid_fragment}{"n_f_left"}=$n_f_left;
				$delta_regions{$current_valid_fragment}{"n_r_right"}=$n_r_right;
			}
			warn "Writing strand switches to outputfile for $bam : chromosome $chromosome \n";

			# Get sorted fragments based on their linked absolute deltafrvalue (high to low absolute deltaRF)
			my @order = sort{abs($deltas{$b}) <=> abs($deltas{$a})} keys %deltas; 
			
			# For each fragment:
			while (@order){
				my $top_delta_fr_fragment = shift(@order);
				next unless exists($deltas{$top_delta_fr_fragment}); # Values may have been removed as a neighbour
				my $top_delta_fr = $deltas{$top_delta_fr_fragment}; # Retrieve deltafrvalue
				
				# Retrieve strand switch location
				my $strandswitch_start = $delta_regions{$top_delta_fr_fragment}{"start"};
				my $strandswitch_stop = $delta_regions{$top_delta_fr_fragment}{"end"};
				if ($strandswitch_start > $strandswitch_stop){ # If reads overlapped, assign the middle basepair of the overlap to the strandswitch
					my $mean = int(($strandswitch_start + $strandswitch_stop)/2);
					$strandswitch_start = $mean;
					$strandswitch_stop = $mean;
				}
				# Get the number of lagging strand reads to each side
				my $strandswitch_left_f = $delta_regions{$top_delta_fr_fragment}{"n_f_left"};
				my $strandswitch_right_r = $delta_regions{$top_delta_fr_fragment}{"n_r_right"};
				
				# Print strandswitch info
				my $ratio = sprintf("%.2f", $top_delta_fr / $scale ); # Get easily comprehendable ratio instead of deltafrvalue
				print F join("\t", $chromosome, $strandswitch_start, $strandswitch_stop, $ratio, $top_delta_fr, $strandswitch_left_f, $strandswitch_right_r), "\n";
				
				# Keep stats
				$delta_stats{$ratio}++;
				$i++;
				
				# Remove previous $scale-1 valid read's deltafrvalues and next $scale-1 valid read's deltafrvalues to prevent calling multiple strandswitches within the same strandswitch region
				foreach my $j ($top_delta_fr_fragment - (1 * $scale - 1) .. $top_delta_fr_fragment + ($scale - 1)) {
					next unless exists($deltas{$j});
					delete $deltas{$j};
				}    
			}
		}
		close F;
	
		# Write stats for BAM file to outputfile
		open F1, ">".$bam.'_'.$scale.$randomize.'deltarf_stats.txt'; # Open outputfile
		foreach my $delta_fr_ratio (sort {$a<=>$b} keys %delta_stats){
			print F1 join("\t", $delta_fr_ratio, $delta_stats{$delta_fr_ratio}, 100*$delta_stats{$delta_fr_ratio}/$i), "\n";
		}
		close F1;
	}
}
warn "Done \n"; 
