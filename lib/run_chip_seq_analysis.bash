#!/bin/bash

# Set the variables
bam_file_name="/data/illumina_pipeline/aligned_experiments/DM290/dm290.bam"
feature_file_name="/data/data2/jab112/2015_belsky_dissertation/chapter2/abf1_analysis/feature_files/abf1_chip_seq_dm290_macs_motif_midpoint_feature_file.csv"
output_file_name="/data/data2/jab112/2015_belsky_dissertation/chapter2/abf1_analysis/output_files/abf1_chip_seq_dm290_around_abf1_midpoint_motif.csv"
# output_file_name="${output_file_name}SRR1517820_h3_chip_exo_around_xu_tss_20_bw.csv"

# Set the parameters
program="ChIPSeqFeatureDensitySignal"
win=500
chip_shift=75
bw=20
span=1

# Set the jar file path
jar_file_path="/data/illumina_pipeline/scripts/java_scripts/jar_files/chip-seq-analysis.jar"

# Execute the script
java -jar $jar_file_path $program $bam_file_name $feature_file_name $output_file_name $win $chip_shift $bw $span
