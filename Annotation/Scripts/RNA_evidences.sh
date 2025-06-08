#!/bin/bash

#Downloads, processes, aligns, and integrates SRA RNA-seq datasets for BRAKER

#Set the paths
THREADS=8
GENOME=/mnt/Franklin/amas/Abeoforma_genome_v2.fasta.masked #masked genome
WORKDIR=/mnt/Franklin/amas/RNAseq
ALIGNDIR=${WORKDIR}/aligned_Vika
TRIMDIR=${WORKDIR}/published
HISAT_INDEX=/home/amas/annotation/abeoforma_masked_index
PORTCULLIS_OUT=/home/amas/annotation/portcullis_output
PORTCULLIS_HINTS=${PORTCULLIS_OUT}/portcullis.hints.gff
MERGED_BAM=${ALIGNDIR}/all_merged_sorted.bam

source ~/miniconda3/bin/activate

cd ${TRIMDIR}
conda activate /home/shared/envs/pacbio

echo "Downloading SRA: ${SRA_ID}"
prefetch ${SRA_ID}
fasterq-dump --split-files ${SRA_ID} -O ${TRIMDIR}

#Split FASTC
#SRA_ID_1.fastq → forward reads
#SRA_ID_2.fastq → reverse reads

echo "Quality check before trimming"
fastqc ${SRA_ID}_1.fastq
fastqc ${SRA_ID}_2.fastq

echo "Trimming adapters and low-quality bases"
trimmomatic PE -threads ${THREADS} -phred33 \
  ${SRA_ID}_1.fastq ${SRA_ID}_2.fastq \
  ${SRA_ID}_forward_paired.fq.gz ${SRA_ID}_forward_unpaired.fq.gz \
  ${SRA_ID}_reverse_paired.fq.gz ${SRA_ID}_reverse_unpaired.fq.gz \
  ILLUMINACLIP:/home/shared/envs/atacseq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36

echo "Quality check after trimming"
fastqc ${SRA_ID}_forward_paired.fq.gz
fastqc ${SRA_ID}_reverse_paired.fq.gz

conda deactivate
conda activate /home/shared/envs/braker

echo "Aligning with HISAT2"
hisat2 -x ${HISAT_INDEX} \
  -1 ${SRA_ID}_forward_paired.fq.gz \
  -2 ${SRA_ID}_reverse_paired.fq.gz \
  -p ${THREADS} -S ${SRA_ID}.sam

echo "Converting and sorting BAM"
samtools view -@ ${THREADS} -bS ${SRA_ID}.sam > ${SRA_ID}.bam
samtools sort -@ ${THREADS} ${SRA_ID}.bam -o ${SRA_ID}_sorted.bam
cp ${SRA_ID}_sorted.bam ${ALIGNDIR}/

done

echo "Merging all BAMs..."
samtools merge -@ ${THREADS} ${MERGED_BAM} ${ALIGNDIR}/*_sorted.bam
samtools index ${MERGED_BAM}

conda deactivate
conda activate /home/shared/envs/portcullis/

echo "Running Portcullis on merged BAM"
portcullis full \
  -t ${THREADS} \
  -o ${PORTCULLIS_OUT} \
  ${GENOME} \
  ${MERGED_BAM}

echo "Converting Portcullis BED to GFF"
python /home/amas/annotation/convert_portcullis.py \
  ${PORTCULLIS_OUT}/3-filt/portcullis_filtered.pass.junctions.bed \
  ${PORTCULLIS_HINTS}

conda deactivate

echo "Done. Use ${PORTCULLIS_HINTS} as hint file for BRAKER."

