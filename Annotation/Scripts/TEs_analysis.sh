#!/usr/bin/env bash
#Genomic Distribution of Repeats
#Identify where TEs are located relative to genes (exonic, intronic, intergenic).

# TE annotations → from RepeatMasker (.gff file)
# Gene annotations → from BRAKER or PASA (.gff3)

# Paths to input files
TE_GFF="/home/amas/annotation/masked2/Abeoforma_genome_v2.fasta.out.gff"
GENE_GFF="/home/amas/annotation/braker_output8/braker.gff3"
GENOME_FA="/mnt/Franklin/amas/Abeoforma_genome_v2.fasta"

# Output filenames
TE_BED="TEs.bed"
GENES_BED="genes.bed"
TE_GENE_OVERLAPS="TE_gene_overlaps.bed"

conda activate /home/shared/envs/agat

#Convert to BED format first for easier handling
echo "Converting GFFs to BED format..."
agat_convert_sp_gff2bed.pl --gff "$GENE_GFF" -o "$GENES_BED"
agat_convert_sp_gff2bed.pl --gff "$TE_GFF" -o "$TE_BED"

#cp masked2/TEs.bed TE_analysis/

conda deactivate

conda activate /home/shared/envs/braker

#Overlap with genes:  This file will show which TEs overlap with which genes.
echo "Finding overlaps between TEs and genes..."
bedtools intersect -a "$TE_BED" -b "$GENES_BED" -wa -wb > "$TE_GENE_OVERLAPS"

#Summary stats
OVERLAPPING_TES=$(cut -f1,2,3,4 "$TE_GENE_OVERLAPS" | sort | uniq | wc -l)
TOTAL_TES=$(wc -l < "$TE_BED")                                           

#Proportion of TEs overlapping genes = overlapping TEs / total TEs
echo "Number of overlapping TEs: $OVERLAPPING_TES"
echo "Total number of TEs: $TOTAL_TES"
echo "Proportion overlapping with genes: $(awk -v a=$OVERLAPPING_TES -v b=$TOTAL_TES 'BEGIN{print (a/b)*100 "%"}')"


conda deactivate
conda activate /home/shared/envs/agat

#Extract introns from genes:
echo "Extracting introns from gene annotation..."

agat_sp_add_introns.pl --gff "$GENE_GFF" -o genes_with_introns.gff3 # introns may be stored here

grep -P "\tintron\t" genes_with_introns.gff3 > introns_only.gff3
agat_convert_sp_gff2bed.pl --gff introns_only.gff3 -o introns.bed  # Convert extracted introns to BED for overlap analysis

conda deactivate
conda activate /home/shared/envs/braker

echo "Calculating number of TEs in introns..."
INTRONIC_TES=$(bedtools intersect -a "$TE_BED" -b introns.bed -u | wc -l)

echo "Calculating number of TEs in intergenic regions..."
# Get genome length BED for subtraction
cut -f1,2 "$GENOME_FA.fai" > genome_file.genome

bedtools makewindows -g genome_file.genome -w 1000000 > genome_windows.bed # chunked genome
bedtools subtract -a genome_windows.bed -b "$GENES_BED" > intergenic_regions.bed

INTERGENIC_TES=$(bedtools intersect -a "$TE_BED" -b intergenic_regions.bed -u | wc -l)

echo "Intronic TEs: $INTRONIC_TES"

echo "Intergenic TEs: $INTERGENIC_TES"

echo "Proportion intronic: $(awk -v a=$INTRONIC_TES -v b=$TOTAL_TES 'BEGIN{print (a/b)*100 "%"}')"

echo "Proportion intergenic: $(awk -v a=$INTERGENIC_TES -v b=$TOTAL_TES 'BEGIN{print (a/b)*100 "%"}')"


conda deactivate


