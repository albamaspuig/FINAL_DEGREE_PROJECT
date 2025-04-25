#!/bin/bash

# Set number of threads
T=16

# Define paths
GENOME="/mnt/Franklin/Genomes/Abeoforma/Abeoforma_genome_v2.fasta"
MASKED_GENOME="/mnt/Franklin/Genomes/Abeoforma/Abeoforma_genome_v2_masked.fasta"
READS_PATH="/home/amas/annotation/RNAseq/trimmed_reads"
PROTEIN_DB="/home/amas/annotation/abeoforma_db.fasta"
SPECIES="abeoforma"
WORKDIR="braker3"
BUSCO_LINEAGE="eukaryota_odb10"

mkdir -p ${WORKDIR}/masked
mkdir -p ${WORKDIR}/braker3_out
mkdir -p ${WORKDIR}/busco_final_results

echo "Initializing environment..."
source ~/miniconda3/bin/activate
conda activate /home/shared/envs/braker

echo "Step 1: Building repeat database..."
BuildDatabase -name abeoforma_db -engine ncbi ${GENOME}

echo "Step 2: Running RepeatModeler to identify repetitive elements..."
RepeatModeler -database abeoforma_db -engine ncbi -pa ${T} -LTRStruct
        # -dir repeats: the directory (repeats) where RepeatModeler store 
output. A set of consensus sequ>
        # -pa 8: number of threads
        # -LTRStruct: enables LTR (Long Terminal Repeat) structure 
detection. It allows RepeatModeler to >
        # -database: The input genome file for RepeatModeler.

gt suffixerator -db 
/mnt/Franklin/Genomes/Abeoforma/Abeoforma_genome_v2.fasta \
  -indexname abeoforma_db \
  -tis -sds -lcp -des -ssp -suf -dna

gt ltrharvest -index abeoforma_db -seqids yes -out ltrharvest_output.txt
        #output: ltrharvest_output.txt ← contains all the predicted LTR 
retrotransposons in GFF-like form>
        #Convert this file to .fa to merge with the RepeatModeler output 
and run it into RepeatMasker

cat abeoforma_db-families.fa ltrharvest_LTRs.fa > 
abeoforma_combined_repeats.fa

echo "Step 3: Masking genome using RepeatMasker..."
RepeatMasker -pa ${T} -lib abeoforma_db-families.fa -gff -xsmall -dir 
${WORKDIR}/masked ${GENOME}
        # The input can either be the original genome (before masking) or 
the genome after running Repeat>
        # -dir repeats: output where results stored as masked genome files 
and annotation files.
        # -gff: output the repetitive elements in GFF3 format (General 
Feature Format)(information about >
        # -species fungi: Specifies the species for repeat classification.
echo "Repeat masking completed."

echo "Step 4: Building HISAT2 genome index..."
hisat2-build -p ${T} ${MASKED_GENOME} abeoforma_masked_index

portcullis full \
  -t ${T} \
  -o portcullis_out \
  masked_genome.fa \
  aligned_reads.bam

        #masked_genome.fa ( RepeatMasker output)
        #aligned_reads.bam (sorted and indexed BAM from HISAT2)


echo "Step 5: Preparing RNA-Seq read files..."
FORWARD_READS=$(ls ${READS_PATH}/*_p_trimmed.fastq | grep 
"_1_.*_trimmed.fastq" | tr '\n' ',')
REVERSE_READS=$(ls ${READS_PATH}/*_p_unpaired.fastq | grep 
"_2_.*_unpaired.fastq" | tr '\n' ',')

FORWARD_READS=${FORWARD_READS%,}
REVERSE_READS=${REVERSE_READS%,}

echo "Step 6: Running HISAT2 for RNA-Seq alignment..."
unset PERL5LIB

hisat2 -p ${T} -x abeoforma_index \
    -1 $FORWARD_READS \
    -2 $REVERSE_READS \
    -S abeoforma_rna.sam

echo "Step 7: Converting and sorting SAM file..."
samtools view -@ ${T} -bS abeoforma_rna.sam | samtools sort -@ ${T} -o 
abeoforma_rna.sorted.bam

echo "Step 8: Indexing BAM file..."
samtools index abeoforma_rna.sorted.bam

echo "RNA-Seq alignment completed."

conda deactivate

echo "Switching to BRAKER3 environment..."
conda activate /home/shared/envs/braker3

export AUGUSTUS_BIN_PATH="/usr/bin/"
export AUGUSTUS_SCRIPTS_PATH="/usr/share/augustus/scripts/"
export 
PROTHINT_PATH="/home/shared/soft/GeneMark-ETP/bin/gmes/ProtHint/bin/"
export GENEMARK_PATH="/home/shared/soft/GeneMark-ETP/bin/"
export 
COMPLEASM_PATH="/home/shared/envs/braker3/lib/python3.7/site-packages"
export PATH=/home/shared/envs/braker/bin:$PATH
export 
PERL5LIB=/home/shared/envs/braker3/lib/perl5/5.32/site_perl:$PERL5LIB

echo "Step 9: Running BRAKER3 for gene prediction..."
braker.pl --species=abeoforma \
    --genome=${GENOME} \
    --prot_seq=${PROTEIN_DB} \
    --bam=abeoforma_rna.sorted.bam \
    --workingdir=${WORKDIR} \
    --threads=${T} \
    --busco_lineage=${BUSCO_LINEAGE} \
    --gff3

# Deactivate BRAKER environment when finished
conda deactivate

echo "Gene prediction completed."


