#!/usr/bin/perl -w
use strict;
use warnings;
use Parallel::ForkManager;
#use Math::Round;

# Script builds a matrix for input into pheatmap (R) 
# Input: Script expects inputfile *_results.txt as put out by map_histonseq_around_oris_sliding_windows.pl 
# Output: A matrix with the origin ID in the first column, the deltafr (represents origin strentgh) in the second column, and the percentage of reads from forward strand for each bin around the origin in the following columns.
#
# Author: Thamar Jessurun Lobo

# Paths and parameters
my @inputfiles = ("PATH(S)_TO_RESULTS.TXT"); # Insert path(s) to *_results.txt (map_histonseq_around_oris_sliding_windows.pl)
my $outpath = "./"; # Insert path to output directory
my $regionsize = 200000; # Define regionsize (in basepairs) to include to each side of the origin center. Default is 200000. MUST be the same as $region - ($windowsize - $slide)/2 in map_histonseq_around_oris_sliding_windows.pl
my $binsize = 1000; # Binsize in basepairs. Default is 1000. MUST be the same as $slide in map_histonseq_around_oris_sliding_windows.pl
my $addtobin = $regionsize/$binsize;

# Get parallel sessions
#my $pm = Parallel::ForkManager->new(1);  # Uncomment to run in parallel

# For each inputfile
foreach my $inputfile (@inputfiles){
    warn "Working on $inputfile\n";
    
    # Run in parallel
    #$pm->start and next;  # Uncomment to run in parallel
    
    # Open outputfile
    my $filen = $inputfile;
    $filen =~ s/^.*\///g;
    $filen =~ s/results.txt/matrix\.txt/;
    my $outputfile = $outpath.$filen;
    open F1, ">".$outputfile;
    
    # Define vars
    my %results;
    my @bins;
    my $current_id;
    my $current_deltafr;
    my $id;
    my $deltafr;
    
    open F, $inputfile;
    while( <F> ){
        chomp;
        my ($chrom, $ori_start, $ori_stop, $deltafr, $pvalue, $fdr, $bin, $n_f, $n_r, $frvalue) = split /\t/;

        # Create unique ID for origin
        my $id = join("_", $chrom, $ori_start, $ori_stop) ;
        if (not defined $current_id){
            $current_id = $id;
            $current_deltafr = $deltafr;
            @bins = (-200..200); # If $regionsize and/or $binsize are changed, replace (-)200 in (-200..200) by $regionsize/$binsize
            print F1 join("\t","id", "deltafr", @bins), "\n"; # Print header to outputfile 
        }elsif($current_id ne $id){  # Information belong to the next origin
            
            # Save percentage of reads from forward strand for empty bins for previous origin
            foreach my $empty_bin (@bins){
                next if $empty_bin eq "present";  # Bin was present, meaning it was not empty
                $results{$empty_bin}="NA";  # Bin was not present meaning it was empty, the percentage of reads from forward strand is "NA"
            }
            
            # Print percentage of reads from forward strands to outputfile for previous origin
            print F1 join( "\t", $current_id, $current_deltafr);
            foreach (sort { $a <=> $b } keys(%results) ){
                my $i = $results{$_};
                print F1 "\t", $i;
            }
            print F1 "\n";
            
            # Reset vars
            @bins = (-200..200);  # If $regionsize and/or $binsize are changed, replace (-)200 in (-200..200) by $regionsize/$binsize
            %results=();
            $current_id = $id;
            $current_deltafr = $deltafr;
        }
        
        # Save percentage of reads from forward strand
        $results{$bin}= ($n_f/($n_f + $n_r)) * 100; 
        
        # Set bin to "present" to indicate that it was not empty (since half the bins are negative but positive in the array index we need to use $addtobin)
        $bins[$bin + $addtobin] = "present";
        
        # Next if not eof
        next if not eof;
        
        warn "Writing info from last origin for $inputfile!\n";
        
        # If eof, print percentage of reads from forward strand for last origin to outputfile
        foreach my $empty_bin (@bins){
            next if $empty_bin eq "present"; # Bin was present, meaning it was not empty
            $results{$empty_bin}="NA"; # Bin was not present meaning it was empty, the percentage of reads from forward strand is "NA"
        }
        print F1 join("\t", $id, $deltafr);
        foreach (sort { $a <=> $b } keys(%results) ){
            my $i = $results{$_};
            print F1 "\t", $i;
        }
        print F1 "\n";
    }
    close F1;
    close F;
    #$pm->finish;  # Uncomment to run in parallel
}
#$pm->wait_all_children;  # Uncomment to run in parallel
warn "Done!\n";
