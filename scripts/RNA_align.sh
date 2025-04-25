#!/bin/bash

# This script aligns RNA-seq reads to a masked genome using HISAT2, 
converts SAM files to BAM format,
# sorts the BAM files, and merges them into a single BAM file. Finally, it 
indexes the merged BAM file.

# Set the path to your masked genome index
GENOME_INDEX="abeoforma_masked_index"

# Set the number of threads for alignment
THREADS=8

# Create an output directory for the aligned reads
OUTPUT_DIR="RNAseq/aligned_RNA_reads"
mkdir -p $OUTPUT_DIR

# Loop through each pair of forward and reverse reads in RNAseq/trim
for forward in RNAseq/trim/*_forward_paired.fq.gz; do
    base=$(basename $forward _forward_paired.fq.gz)
    reverse="RNAseq/trim/${base}_reverse_paired.fq.gz"

    # Step 1: Align the RNA-seq reads to the masked genome
    hisat2 -p $THREADS -x $GENOME_INDEX -1 $forward -2 $reverse -S 
${OUTPUT_DIR}/${base}.sam
done

# Path to the aligned RNA-seq files
ALIGNED_DIR="RNAseq/aligned_RNA_reads"

# Convert SAM to BAM and sort
for samfile in $ALIGNED_DIR/*.sam; do
    # Convert .sam to .bam
    bamfile="${samfile%.sam}.bam"
    samtools view -bS $samfile > $bamfile

    # Sort the .bam file
    sorted_bamfile="${bamfile%.bam}.sorted.bam"
    samtools sort $bamfile -o $sorted_bamfile
done

# Merge sorted BAM files
samtools merge $ALIGNED_DIR/merged_RNA.bam $ALIGNED_DIR/*.sorted.bam

# Index the merged BAM file
samtools index $ALIGNED_DIR/merged_RNA.bam

