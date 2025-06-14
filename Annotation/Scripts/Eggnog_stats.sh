#!/bin/bash
# Extract and summarize key functional annotation (GO terms, EC numbers, KEGG pathways, PFAM domains)

conda activate /home/shared/envs/eggnog

mkdir -p /home/amas/annotation/eggnog_data
export EGGNOG_DATA_DIR=/home/amas/annotation/eggnog_data
python download_eggnog_data.py

emapper.py \
  -i cleaned.fasta \
  -o merged_eggnog \
  --output_dir /home/amas/annotation/eggnog_output4 \
  --itype proteins \
  --override \
  --cpu 8 
  
conda deactivate 

ANNOT="merged_eggnog.emapper.annotations"

# Remove comment lines and header
grep -v "^#" $ANNOT > cleaned_annotations.tsv

# Extract total number of genes annotated
total_genes=$(wc -l < cleaned_annotations.tsv)
echo "Total genes annotated: $total_genes"

# Extract and count unique GO terms
echo "Extracting GO terms..."
grep -v "^#" $ANNOT | awk -F"\t" '$10 != "-" {print $10}' | tr ',' '\n' | sort | uniq -c | sort -nr > go_term_counts.txt
echo "Top 10 GO terms:"
head -10 go_term_counts.txt

# Extract and count unique EC numbers
echo "Extracting EC numbers..."
grep -v "^#" $ANNOT | awk -F"\t" '$11 != "-" {print $11}' | tr ',' '\n' | sort | uniq -c | sort -nr > ec_counts.txt
echo "Top 10 EC numbers:"
head -10 ec_counts.txt

# Extract and count unique KEGG Orthologs (KO)
echo "Extracting KEGG Orthologs (KO)..."
grep -v "^#" $ANNOT | awk -F"\t" '$12 != "-" {print $12}' | tr ',' '\n' | sort | uniq -c | sort -nr > kegg_ko_counts.txt
echo "Top 10 KEGG Orthologs:"
head -10 kegg_ko_counts.txt

# Extract and count unique KEGG Pathways
echo "Extracting KEGG Pathways..."
grep -v "^#" $ANNOT | awk -F"\t" '$13 != "-" {print $13}' | tr ',' '\n' | sort | uniq -c | sort -nr > kegg_pathway_counts.txt
echo "Top 10 KEGG Pathways:"
head -10 kegg_pathway_counts.txt

# Extract and count PFAM domains
echo "Extracting PFAM domains..."
grep -v "^#" $ANNOT | awk -F"\t" '$21 != "-" {print $21}' | tr ',' '\n' | sort | uniq -c | sort -nr > pfam_counts.txt
echo "Top 10 PFAM domains:"
head -10 pfam_counts.txt

#predicted proteins across COG functional categories.
awk -F"\t" '$7 != "-" {print $7}' cleaned_annotations.tsv | tr -d ' ' | fold -w1 | sort | uniq -c | sort -nr > cog_function_counts.txt



cut -f10 cleaned_annotations.tsv | tr ',' '\n' | grep '^GO:' | sort | uniq -c | sort -nr | awk '{print $2 "\t" $1}' > REVIGO_input.txt

head REVIGO_input.txt
'''
GO:0008150	6711
GO:0005575	6650
GO:0044464	6549
GO:0005623	6549
GO:0005622	6307
GO:0044424	6284
GO:0009987	6009
GO:0003674	6000
GO:0043226	5645
GO:0043229	5599
'''
