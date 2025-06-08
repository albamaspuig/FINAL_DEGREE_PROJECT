#!/bin/bash

# Set number of threads
T=8

# Define paths
GENOME="/mnt/Franklin/amas/Abeoforma/Abeoforma_genome_v2.fasta"
MASKED_GENOME="/mnt/Franklin/amas/Abeoforma_genome_v2.fasta.masked"
READS_PATH="/mnt/Franklin/amas/RNAseq/aligned_Vika/all_merged_sorted.bam"
PROTEIN_DB="/home/amas/annotation/proteins_clean.fa"
PORTCULLIS_HINTS="/home/amas/annotation/portcullis_output/portcullis.hints.gff"
SPECIES="abeoforma"
WORKDIR="braker3"
BUSCO_LINEAGE="eukaryota_odb10"o

mkdir -p ${WORKDIR}/masked

#Generate the proteomes database combining PASA with new proteomes database
cat /home/amas/annotation/combined_proteomes2.fasta /mnt/Franklin/amas/Abeoforma_PASA.agat.pep.fasta > prot_hint.fasta

#error, an unexpected symbol was found in the file proteins.fa in: NSR..QHARKCNR.QFKWTII.STKSKPTTKSKSTTKSKSTTKSKSTTKSK.LC.P.KPSGNY.WNDQSK
#error in protein file parsing: proteins.fa

awk '/^>/ {print; next} {gsub(/[^ACDEFGHIKLMNPQRSTVWYBXZ\*]/,""); print}' prot_hint.fasta > proteins_clean.fa

echo "Initializing environment..."
source ~/miniconda3/bin/activate
conda activate /home/shared/envs/braker

echo "Step 1: Building repeat database..."
BuildDatabase -name abeoforma_db -engine ncbi ${GENOME}

echo "Step 2: Running RepeatModeler to identify repetitive elements..."
RepeatModeler -database abeoforma_db -engine ncbi -pa ${T} -LTRStruct

# Suffixerator for repeat index
gt suffixerator -db  ${GENOME} \
  -indexname abeoforma_db \
  -tis -sds -lcp -des -ssp -suf -dna

# Harvest LTRs
gt ltrharvest -index abeoforma_db -seqids yes -out ltrharvest_output.txt

# Combine RepeatModeler and LTR harvest outputs
cat abeoforma_db-families.fa ltrharvest_LTRs.fa > abeoforma_combined_repeats.fa

echo "Step 3: Masking genome using RepeatMasker..."
RepeatMasker -pa ${T} -lib abeoforma_db-families.fa -gff -xsmall -dir ${WORKDIR}/masked ${GENOME}

echo "Repeat masking completed."

echo "Step 4: Building HISAT2 genome index..."
hisat2-build -p ${T} ${MASKED_GENOME} abeoforma_masked_index

echo "Step 5: Running HISAT2 for RNA-Seq alignment..."
# Running the external RNA alignment pipeline
./rna_aligning.sh  

# Deactivate HISAT2 environment
conda deactivate

# Activate Portcullis environment
conda activate /home/shared/envs/portcullis/

# Run Portcullis
echo "Step 6: Running Portcullis..."
portcullis full \
  -t 8 \
  -o ${WORKDIR}/portcullis_output \
  ${MASKED_GENOME}  \
  ${READS_PATH}

# Convert Portcullis junctions to GFF format
echo "Converting Portcullis output to GFF..."
python convert_portcullis.py \
  ${WORKDIR}/portcullis_output/3-filt/portcullis_filtered.pass.junctions.bed \
  ${WORKDIR}/portcullis_output/portcullis.hints.gff

conda deactivate

# Activate BRAKER3 environment
echo "Switching to BRAKER3 environment..."
conda activate /home/shared/envs/braker3

# Set necessary paths for AUGUSTUS and GeneMark
export AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config/"
export AUGUSTUS_BIN_PATH="/usr/bin/"
export AUGUSTUS_SCRIPTS_PATH="/usr/share/augustus/scripts/"
export PROTHINT_PATH="/home/shared/soft/GeneMark-ETP/bin/gmes/ProtHint/bin/"
export GENEMARK_PATH="/home/shared/soft/GeneMark-ETP/bin/"
export COMPLEASM_PATH="/home/shared/envs/braker3/lib/python3.7/site-packages"
export PATH=/home/shared/envs/braker/bin:$PATH
export PERL5LIB=/home/shared/envs/braker3/lib/perl5/5.32/site_perl

# Run BRAKER3 for gene prediction
echo "Step 7: Running BRAKER3 for gene prediction..."
braker.pl \
	--workingdir=/mnt/Franklin/amas/braker_output4 \
	--genome=${MASKED_GENOME} \
	--bam=${READS_PATH} \
	--prot_seq=${PROTEIN_DB} \
	--hints=${PORTCULLIS_HINTS} \
	--AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config/ \
	--AUGUSTUS_BIN_PATH=/usr/bin/ \
	--AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts/ \
	--PROTHINT_PATH=/home/shared/soft/GeneMark-ETP/bin/gmes/ProtHint/bin/ \
	--species=${SPECIES} \
	--etpmode \ #Doesn't exist in braker3
	--softmasking \
	--useexisting \
	--threads=8
---------------------------------------------------------------------------------------------
#Re run braker adding the proteins from PASA
braker.pl --workingdir=/home/amas/annotation/braker_output4 --genome=/home/amas/annotation/masked/Abeoforma_genome_v2.fasta.masked --bam=/mnt/Franklin/amas/RNAseq/aligned_Vika/all_merged_sorted.bam --prot_seq=/home/amas/annotation/proteins_clean_final.fa --hints=/home/amas/annotation/portcullis_output/portcullis.hints.gff --AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config/ --AUGUSTUS_BIN_PATH=/usr/bin/ --AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts/ --PROTHINT_PATH=/home/shared/soft/GeneMark-ETP/bin/gmes/ProtHint/bin/ --species="Abeoforma" --useexisting --threads=8


#If this doesn't work try adding published RNA or les restricive masking of the genome.

echo "Gene prediction completed."

