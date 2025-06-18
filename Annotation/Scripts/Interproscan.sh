#!/bin/bash
# functional_annotation_interproscan.sh
#
# Functional annotation of predicted proteins using InterProScan.
# Predicts protein domains and families, assigns Gene Ontology (GO) terms, and provides detailed domain descriptions.

set -euo pipefail

# Arguments
INPUT_FASTA=${1:-"merged_proteins_clean.fasta"}
OUTPUT_DIR=${2:-"interpro_output"}
NUM_CHUNKS=${3:-10}
THREADS=${4:-8}

# Activate environment
source ~/miniconda3/bin/activate
conda activate /home/shared/envs/interpro

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "Running InterProScan on: $INPUT_FASTA"
echo "Splitting FASTA into $NUM_CHUNKS parts..."

# Split FASTA
seqkit split -p "$NUM_CHUNKS" "$INPUT_FASTA" -O split_fastas

# Run InterProScan on each chunk
for file in split_fastas/*.fasta
do
    part=$(basename "$file" .fasta)
    echo "Processing $part ..." | tee -a interproscan_all.log
    interproscan.sh -i "$file" -f tsv -goterms -pa -dp -cpu "$THREADS" > "${part}.tsv" 2>> interproscan_all.log
done

echo "Merging output files..."
cat *.tsv > merged_interproscan.tsv

echo "Generating summary statistics..."

# Top 20 domains (InterPro accessions)
awk -F'\t' '($13 != "-") {print $13}' merged_interproscan.tsv | sort | uniq -c | sort -nr | head -20 > top20_domains.txt

# Counts per domain type (e.g., PFAM, SMART)
awk -F'\t' '($5 != "-") {print $5}' merged_interproscan.tsv | sort | uniq -c | sort -nr > domain_type_counts.txt

# List each protein with its InterPro domain and description
awk -F'\t' '($13 != "-") {print $1"\t"$5"\t"$13}' merged_interproscan.tsv | sort | uniq > gene_domain_descriptions.tsv

echo "Functional annotation with InterProScan completed successfully."
echo "Output directory: $OUTPUT_DIR"

conda deactivate

