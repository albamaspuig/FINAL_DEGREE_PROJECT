#!/usr/bin/env bash
# Merge annotations

# Paths
GENOME_FA="/mnt/Franklin/amas/Abeoforma_genome_v2.fasta"
BRAKER_GFF="/home/amas/annotation/braker_output9/braker.gff3"
PASA_GFF="/mnt/Franklin/amas/Abeoforma_PASA.agat_v2.gff"
BUSCO_DB="/home/amas/annotation/eukaryota_odb10"

# Output directories
OUTDIR="$HOME/annotation"
MERGE_DIR="${OUTDIR}/merged"
EVM_DIR="${OUTDIR}/EVM"

mkdir -p "$MERGE_DIR" "$EVM_DIR"

#1. Merge annotations with AGAT
conda activate /home/shared/envs/agat

agat_sp_merge_annotations.pl --gff "$BRAKER_GFF" --gff "$PASA_GFF" --out "${MERGE_DIR}/merged_annotation.gff3"

'''
final result:
There is 75868 intron
There is 271083 exon
There is 10679 three_prime_utr
There is 9981 start_codon
There is 264855 cds
There is 19225 gene
There is 9996 stop_codon
There is 16890 five_prime_utr
There is 31482 mrna
'''

#Extract proteins from merged annotation

conda deactivate
conda activate /home/shared/envs/TEs

gffread "${MERGE_DIR}/merged_annotation.gff3" \
    -g "$GENOME_FA" \
    -y "${MERGE_DIR}/merged_proteins.fasta"


# Remove **exact** duplicate protein sequences
seqkit rmdup -s "${MERGE_DIR}/merged_proteins.fasta" > "${MERGE_DIR}/merged_proteins_nr.fasta"
#[INFO] 9681 duplicated records removed

grep -c "^>" merged_proteins.fasta
#31482

grep -c "^>" merged_proteins_nr.fasta
#21801

conda deactivate
conda activate /home/shared/envs/busco

# Run BUSCO on the non-redundant proteins to assess completeness
conda deactivate
conda activate /home/shared/envs/busco

busco -i "${MERGE_DIR}/merged_proteins_nr.fasta" -l "$BUSCO_DB" -o "${MERGE_DIR}/busco_merged_proteins_nr" -m protein -f --offline

'''
    ---------------------------------------------------
    |Results from dataset eukaryota_odb10              |
    ---------------------------------------------------
    |C:94.1%[S:80.8%,D:13.3%],F:2.0%,M:3.9%,n:255      |
    |240    Complete BUSCOs (C)                        |
    |206    Complete and single-copy BUSCOs (S)        |
    |34    Complete and duplicated BUSCOs (D)          |
    |5    Fragmented BUSCOs (F)                        |
    |10    Missing BUSCOs (M)                          |
    |255    Total BUSCO groups searched                |
    ---------------------------------------------------

'''

conda deactivate

# 2.  Integration with EVM (EvidenceModeler)
conda activate /home/shared/envs/agat

cd "$EVM_DIR"

agat_convert_sp_gxf2gxf.pl --gff "$PASA_GFF" -o pasa_evm.gff3
agat_convert_sp_gxf2gxf.pl --gff "$BRAKER_GFF" -o braker_evm.gff3

#Check name for weights
awk '$3 == "gene" {print $2}' pasa_evm.gff3 | sort | uniq
'''
.
AUGUSTUS
'''

awk '$3 == "gene" {print $2}' braker_evm.gff3 | sort | uniq
'''
GeneMark.hmm3   
gmst 
'''

# Prepare weights File
cat > weights.txt <<EOF
TRANSCRIPT            .               10
TRANSCRIPT            AUGUSTUS        5
ABINITIO_PREDICTION   GeneMark.hmm3   1
ABINITIO_PREDICTION   gmst            1
EOF

conda deactivate
conda activate /home/shared/envs/braker

#Partition the genome 
/home/shared/envs/braker/bin/EvmUtils/partition_EVM_inputs.pl \
    --partition_dir "$EVM_DIR" \
    --genome "$GENOME_FA" \
    --gene_predictions braker_evm.gff3 \
    --transcript_alignments pasa_evm.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000 \
    --partition_listing partitions_list.out

#Run EVM
/home/shared/envs/braker/bin/EvmUtils/write_EVM_commands.pl \
    --genome "$GENOME_FA" \
    --weights "$EVM_DIR/weights.txt" \
    --gene_predictions braker_evm.gff3 \
    --transcript_alignments pasa_evm.gff3 \
    --output_file_name evm.out \
    --partitions partitions_list.out > commands.list
 
# Run EVM parallelized
cat commands.list | parallel -j 4

# Recombine the results into final GFF3
/home/shared/envs/braker/bin/EvmUtils/recombine_EVM_partial_outputs.pl \
    --partitions partitions_list.out \
    --output_file_name evm.out


#Convert to gff3
/home/shared/envs/braker/bin/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
    --partitions partitions_list.out \
    --output evm.out \
    --genome "$GENOME_FA"
    
# Merge GFF files
cat Abeoforma_whisleri*/evm.out > evm.all.out
cat Abeoforma_whisleri*/evm.out.gff3 > evm.out.all.gff

#Extract proteins
gffread evm.out.all.gff3 -g "$GENOME_FA" -y "${EVM_DIR}/evm_proteins.fasta"

# Run BUSCO to assess completeness

conda activate /home/shared/envs/busco

busco -i "${EVM_DIR}/evm_proteins.fasta" \
    -l "$BUSCO_DB" \
    -o "${EVM_DIR}/evm_busco" \
    -m protein -f --offline

conda deactivate

'''
       C:83.1%[S:82.7%,D:0.4%],F:5.9%,M:11.0%,n:255       
        212     Complete BUSCOs (C)                        
        211     Complete and single-copy BUSCOs (S)        
        1       Complete and duplicated BUSCOs (D)         
        15      Fragmented BUSCOs (F)                      
        28      Missing BUSCOs (M)                         
        255     Total BUSCO groups searched  

'''

