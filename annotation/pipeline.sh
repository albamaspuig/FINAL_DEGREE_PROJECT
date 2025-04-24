#!/bin/bash

# Set number of threads
T=16

# Define paths
GENOME="/mnt/Franklin/Genomes/Abeoforma/Abeoforma_genome_v2.fasta"
READS_PATH="/home/amas/annotation/RNAseq/trimmed_reads"
PROTEIN_DB="/home/amas/annotation/abeoforma_db.fasta"
SPECIES="abeoforma"
WORKDIR="braker3_out"
BUSCO_LINEAGE="eukaryota_odb10"

mkdir -p ${WORKDIR}/masked
mkdir -p ${WORKDIR}/braker3_out
mkdir -p ${WORKDIR}/busco_final_results

# Mask repetitive sequences
BuildDatabase -name abeoforma_db -engine ncbi ${GENOME}

# Generate a de novo repeat library
RepeatModeler -database abeoforma_db -engine ncbi -pa ${T} -LTRStruct

# Mask the genome (soft masking)
RepeatMasker -pa ${T} -lib abeoforma_db-families.fa -gff -xsmall -dir ${WORKDIR}/masked ${GENOME}

# Map RNA-Seq reads to the genome
# Build a genome index
hisat2-build -p ${T} ${GENOME} abeoforma_index

# Align reads
# Create a list of forward and reverse read files
FORWARD_READS=$(ls ${READS_PATH}/*_p_trimmed.fastq | grep "_1_.*_trimmed.fastq" | tr '\n' ',')
REVERSE_READS=$(ls ${READS_PATH}/*_p_unpaired.fastq | grep "_2_.*_unpaired.fastq" | tr '\n' ',')

# Remove the last comma from the lists
FORWARD_READS=${FORWARD_READS%,}
REVERSE_READS=${REVERSE_READS%,}

# Run HISAT2 with the file variables
hisat2 -p ${T} -x abeoforma_index \
    -1 $FORWARD_READS \
    -2 $REVERSE_READS \
    -S abeoforma_rna.sam

# Convert SAM to BAM and sort
samtools view -@ ${T} -bS abeoforma_rna.sam | samtools sort -@ ${T} -o abeoforma_rna.sorted.bam

# Index BAM file
samtools index abeoforma_rna.sorted.bam

# Gene Prediction with BRAKER3 
braker.pl --species=abeoforma \
    --genome=${GENOME} \
    --prot_seq=${PROTEIN_DB} \
    --bam=abeoforma_rna.sorted.bam \
    --workingdir=${WORKDIR} \
    --threads=${T} \
    --busco_lineage=${BUSCO_LINEAGE} \
    --gff3  # Output in GFF3 format


