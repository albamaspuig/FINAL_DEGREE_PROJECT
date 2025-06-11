#!/bin/bash

# Set paths
GENOME_INDEX="/home/amas/annotation/index/abeoforma_masked_index"
THREADS=8
OUTPUT_DIR="/mnt/Franklin/amas/RNAseq/aligned_Vika2"
GENOME_FASTA="/home/amas/annotation/masked/Abeoforma_genome_v2.fasta.masked"
PORTCULLIS_OUT="/home/amas/annotation/portcullis_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Aligning RNA-seq reads with HISAT2..."
for f1 in /mnt/Franklin/amas/RNAseq/*_1_p.fastq; do
  f2="${f1/_1_p.fastq/_2_p.fastq}"
  base=$(basename "$f1" _1_p.fastq | sed 's/_1$//')
  sam_out="${OUTPUT_DIR}/${base}.sam"
  echo "  ↪ Processing $base"
  hisat2 -x "$GENOME_INDEX" -1 "$f1" -2 "$f2" -p "$THREADS" -S "$sam_out"
done

echo "Converting SAM to BAM..."
for sam_file in "$OUTPUT_DIR"/*.sam; do
  bam_out="${sam_file%.sam}.bam"
  samtools view -bS "$sam_file" > "$bam_out"
done

echo "Sorting BAM files..."
for bam in "$OUTPUT_DIR"/*.bam; do
  sorted_out="${bam%.bam}_sorted.bam"
  samtools sort "$bam" -o "$sorted_out"
done

echo "Merging all sorted BAM files..."
samtools merge "$OUTPUT_DIR/all_merged_sorted.bam" "$OUTPUT_DIR"/*_sorted.bam

echo "Indexing merged BAM..."
samtools index "$OUTPUT_DIR/all_merged_sorted.bam"

conda activate /home/shared/envs/portcullis
echo "Running Portcullis..."
portcullis full \
  -t "$THREADS" \
  -o "$PORTCULLIS_OUT" \
  "$GENOME_FASTA" \
  "$OUTPUT_DIR/all_merged_sorted.bam"

echo "Pipeline completed successfully!"
