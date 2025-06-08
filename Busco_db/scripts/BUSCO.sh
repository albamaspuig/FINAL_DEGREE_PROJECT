#!/bin/bash

# Genome BUSCO Custom Dataset Creation Script
# Build a custom BUSCO database for unicellular holozoans

# Activate OrthoFinder environment
source ~/miniconda3/bin/activate
conda activate /home/shared/envs/orthofinder/

# Step 1: Run OrthoFinder to find orthogroups
orthofinder -f /mnt/Franklin/amas/busco_input/ -t 4 -a 4

# Expected outputs:
#   OrthoFinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
#   OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/*.fa

conda deactivate

# Step 2: Extract, align, and build HMMs
conda activate /home/shared/envs/evosearch

BUSCO_DIR="/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10"
HMM_DIR="$BUSCO_DIR/hmms"
mkdir -p "$HMM_DIR"

cp /mnt/Franklin/amas/busco_input/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/*.fa "$HMM_DIR"

# Align and build HMMs
for file in "$HMM_DIR"/*.fa; do
    base=$(basename "$file" .fa)
    mafft --auto "$file" > "$HMM_DIR/${base}.aln"
    hmmbuild "$HMM_DIR/${base}.hmm" "$HMM_DIR/${base}.aln"
done

# Step 3: Map HOG IDs to integers
python3 mapping.py

# Rename files to use numeric BUSCO IDs
while IFS=$'\t' read -r new_id old_id; do
    for ext in hmm fa aln; do
        mv "$HMM_DIR/${old_id}.${ext}" "$HMM_DIR/${new_id}.${ext}" 2>/dev/null || true
    done
done < "$BUSCO_DIR/id_mapping.tsv"

# Step 4: Concatenate HMMs and build index
cat "$HMM_DIR"/*.hmm > "$BUSCO_DIR/myholozoa.hmm"
hmmpress "$BUSCO_DIR/myholozoa.hmm"

# Step 5: Create dataset config file
cat << EOF > "$BUSCO_DIR/dataset.cfg"
name=my_busco_dataset_odb10
species=custom
domain=eukaryota
creation_date=$(date +%Y-%m-%d)
number_of_BUSCOs=$(ls "$HMM_DIR"/*.hmm | wc -l)
number_of_species=1
max_intron=100000
max_seq_len=100000
OrthoDB_version=10.1
EOF

# Step 6: Generate lengths_cutoff file
python3 lengthcutoff.py

# Step 7: Generate scores_cutoff file
./scorescutoff.sh

# Step 8: Test the custom BUSCO lineage
hmmsearch --tblout "$BUSCO_DIR/test_results.tbl" "$BUSCO_DIR/myholozoa.hmm" \
    /mnt/Franklin/amas/busco_input/Capsaspora_owczarzaki.faa

# Step 9: Run BUSCO test
busco -i /mnt/Franklin/amas/busco_input/Capsaspora_owczarzaki.faa \
      -o busco_test_output \
      -l "$BUSCO_DIR" \
      -m prot --cpu 4 --offline -f



