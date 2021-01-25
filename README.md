# histone_deposition_analysis
The scripts in this repo were used to analyze double-click-seq and OK-seq data for our manuscript: "Inversion of asymmetric histone deposition upon replication stress"

1. To prevent problems with OK-seq data due to repetitive sequences that are collapsed in the current genome reference (and hence showing extremely high coverage by NGS reads), we first divided the genome into bins of 10 kb and determined the median number of reads mapping to a bin (primary alignments with a minimum mapping quality of 20, PCR duplicates excluded). All reads that mapped to (part of) a bin with either less than a tenth or more than 10 times the median number of reads per bin, were excluded from further analysis. Blacklists for each OK-seq dataset with these bins and the number of reads mapping there can be found in the folder *Blacklist_RPE1_OKseq*.
2. The scripts in the folder *Identification_of_replication_origins_using_OKseq_data* were used to identify replication initiation zones from OK-seq data. First **detect_strandswitches.pl** was used to calculate strand-switching regions by processing BAM files from OK-seq. Next, the output of **detect_strandswitches.pl** was used with **get_significant_oris.R** to call estimated replication initiation zones (strand switched that had a benjamini hochberg adjusted pvalue < 0.01). The resulting initiation zones can be found in **SRX4036932_sorted_mdup_filtered.bam_100ndeltarfs_estimated_oris.txt**.
3. The scripts in the folder *Bias_calculation* were used to map double-click-seq data and OK-seq data at replication initiation zones. Specifically, the script **map_seqdata_around_oris.pl** was used to map double-click-seq data and OK-seq data at replication initiation zones and generate output files that were used to visualize the average histone bias and replication fork direction at replication initiation zones, as well as to generate inputfiles for **bootstrap_bias.pl** in the *Bootstrap_analysis* folder (see point 4). The script **map_seqdata_around_oris_gc.pl** was used to map double-click-seq data and OK-seq data at replication initiation zones that were grouped into quantiles based on GC-percentage or Replication Timing (RT)-score and generate output files that were used to visualize this average histone bias and replication fork direction at these groups of replication initiation zones, as well as to generate inputfiles for **bootstrap_bias_gcquantiles.pl** in the *Bootstrap_analysis* folder (see point 4). The script **map_seqdata_around_oris_sliding_windows.pl** was used to map double-click-seq data and OK-seq data at replication initiation zones with sliding windows to generate files that were used to build heatmaps (visualizing the histone bias and replication fork direction at each replication initiation zone).  The script **convert_to_matrix.pl** was used to build the matrix for input into pheatmap (R) from the output file from **map_seqdata_around_oris_sliding_windows.pl**.
4. The scripts in the folder *Bootstrap_analysis* were used to bootstrap outputfiles of the mapping of double-click-seq data and OK-seq, so that confidence intervals could be determined. Specifically, the script **bootstrap_bias.pl** was used to bootstrap results of **map_seqdata_around_oris.pl** and the script **bootstrap_bias_gcquantiles.pl** was used to bootstrap results of **map_seqdata_around_oris_gc.pl**.
5. The detected autosomal replication initiation zones were grouped into active or inactive regions based on hRPE-1 gene expression as downloaded from GEO accession GSE89413 (mapped to primary genome assembly GRCh38 using STAR aligner (v. 2.7.3) with --outSAMmapqUnique 50 ). For this the script **02a_cluster_expressed.pl** in folder *Active_vs_repressed_region_assignment* was used.  Genes with a minimal expression of 1 FPKM and clusters of neighboring genes with a minimum expression of 1 FPKM and maximum intergenic distance of 10 kb were regarded as active regions. The remaining sequences were marked as inactive regions. A bed file with the active and inactive regions, **20191126_RPE1_active_passive.bed**, can be found in in this folder as well.
