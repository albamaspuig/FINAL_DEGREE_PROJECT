#!/bin/bash
import os

HMM_DIR="/mnt/Franklin/amas/busco_input/uni_holozoans_odb10/hmms"
OUTFILE="/mnt/Franklin/amas/busco_input/uni_holozoans_odb10/scores_cutoff"

echo "Generating scores_cutoff..."
> "$OUTFILE"  # clear or create output file

# Loop over each .hmm file in the HMM directory
for hmm_file in "$HMM_DIR"/*.hmm; do
    base=$(basename "$hmm_file" .hmm)  # Extract the base name of the hmm file (remove directory and extension)
    ref_fasta="$HMM_DIR/$base.fa"  # Define the corresponding reference fasta file path (same basename, .fa extension)

    # Check if the reference fasta file exists; if not, skip this iteration
    if [[ ! -f "$ref_fasta" ]]; then
        echo "Missing reference fasta for $base, skipping..."
        continue
    fi

    # Run hmmsearch of the HMM profile against the reference fasta sequences
    # --tblout outputs results in a tabular format to stdout
    # --noali suppresses alignment output (only summary info)
    # Parse output to:
    #   - exclude comment lines (grep -v '^#')
    #   - extract the 13th column (the full sequence bit score)
    #   - sort scores in descending order
    #   - take the highest score (head -1)
    
    best_score=$(hmmsearch --tblout /dev/stdout --noali "$hmm_file" "$ref_fasta" \
                 | grep -v '^#' | awk '{print $13}' | sort -nr | head -1)

    # If no score was found (empty result), set a fallback score of 50
    if [[ -z "$best_score" ]]; then
        best_score=50  # fallback score if no hit
    fi

    # Append the base HMM ID and its best score to the output file separated by a tab
    echo -e "${base}\t${best_score}" >> "$OUTFILE"
done

echo "Finished writing scores_cutoff to $OUTFILE"


