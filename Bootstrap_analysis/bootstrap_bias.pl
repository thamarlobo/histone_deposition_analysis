#!/usr/bin/perl -w
use strict;
use warnings;
use Parallel::ForkManager;
use POSIX ();

# Script bootstraps *_results.txt outputfiles of map_seqdata_around_oris.pl 
# Input: _results.txt outputfile(s) of map_seqdata_around_oris.pl
# Output: Outputfile(s) *_bootstrapped.txt contains the summed number of fragments from the forward strand and reverse strand per bin for each bootstrapround
# Output format: five tab-separated columns: $boot (=bootstrap round), $b (=bin), $n_f (=summed number of fragments from forward strand), $n_r (=summed number of fragments from reverse strand), $frvalue (=frvalue, not used downstream)
# Outputfile *_bootstrapped.txt can be used to determine the 95 % Confidence Interval by:
# 1. Calculating the average percentage of reads from forward strand ($n_f/($n_f + $n_r)) * 100) per bin for each bootstrapround
# 2. Sorted the resulting percentages of reads from forward strand per bin 
# 3. Selecting the 2.5th and the 97.5th percentiles per bin
#
# Author: Thamar Jessurun Lobo

# Paths and parameters
my @files = ("PATH(S)_TO_RESULTS.TXT"); # Insert path(s) to *_results.txt (map_seqdata_around_oris.pl)
my $n_rounds = 1000; # Number of bootstrap rounds. Default is 1000

# Get parallel sessions
#my $pm = Parallel::ForkManager->new(1);  # Uncomment to run in parallel
   
# For each file   
foreach my $f (@files){

    #$pm->start and next; # Uncomment to run in parallel
    
    # Open file
    open F1, $f;
    
    # Vars
    my $i =0;
    my %dict;
    my %origin_info;
    
    # Save origin info under to be assigned origin ID number
    while (<F1>){
        chomp;
        my ($chr, $ori_start, $ori_stop,  $deltafr, $dummy_pvalue, $dummy_fdr, $b, $n_f, $n_r, $frvalue) = split(/\t/);
        if (exists $dict{$chr."_".$ori_start."_".$ori_stop}){
            my $j =  $dict{$chr."_".$ori_start."_".$ori_stop};
            $origin_info{$j}{$b}{"n_f"}=$n_f;
            $origin_info{$j}{$b}{"n_r"}=$n_r;
            next; # Origin already got an ID number
        }
        $i +=1;
        $dict{$chr."_".$ori_start."_".$ori_stop}=$i; # Assign ID number to origin
        $origin_info{$i}{$b}{"n_f"}=$n_f;
        $origin_info{$i}{$b}{"n_r"}=$n_r;
    }
    close F1;
    
    # Open outputfile
    my $outfn = $f;
    $outfn =~ s/_results\.txt/_bootstrapped\.txt/; 
    warn "outfile is $outfn \n";
    open F2,'>'.$outfn;
    
    # Bootstrap:
    foreach my $boot (1..$n_rounds){ 
       
        # Get random origin IDs to use for bootstrap
        my $n_origins = scalar keys %dict;
        my @origins;
        foreach my $k (1..$n_origins){
            my $random_number = int(rand($n_origins)) + 1;
            push @origins, $random_number;
        }
        my $l = scalar @origins;
        warn "Number of origins to bootstrap is $l out of $n_origins origins in total. \n"; # Check: $l should be equal to $n_origins
    
        # Sum number of fragments per strand across selected origins
        my %summed_diad_counts;
        foreach my $origin (@origins){
            foreach my $b (keys %{$origin_info{$origin}}){
                $summed_diad_counts{$b}{"n_f"} += $origin_info{$origin}{$b}{"n_f"};
                $summed_diad_counts{$b}{"n_r"} += $origin_info{$origin}{$b}{"n_r"};
            }
        }
        # Write summed results to output
        foreach my $b (sort {$a <=> $b} keys %summed_diad_counts){
            my $n_f = 0;
            my $n_r = 0;
            $n_f += $summed_diad_counts{$b}{"n_f"} if exists $summed_diad_counts{$b}{"n_f"};
            $n_r += $summed_diad_counts{$b}{"n_r"} if exists $summed_diad_counts{$b}{"n_r"};
            
            # Calculate f-r/f+r and make it NA if both n_f and n_r are 0
            my $frvalue;
            if($n_f + $n_r == 0){
				$frvalue = "NA";
            }else{
                $frvalue = ($n_f-$n_r)/($n_f+$n_r);
            }
            
            # Write result of bootstrapround to outputfile
            print F2 join("\t", $boot, $b, $n_f, $n_r, $frvalue), "\n";
        }

    }
    close F2;
    #$pm-> finish; # Uncomment to run in parallel
}
#$pm->wait_all_children; # Uncomment to run in parallel
warn "Finished!\n";
