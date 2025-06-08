#!/usr/bin/env python3
import os

aln_dir = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/hmms"
output_file = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/lengths_cutoff"

#check if a line in the file is a FASTA header (starts with ">")
def is_fasta_header(line):
    return line.startswith(">")

with open(output_file, "w") as out:
    # Get a sorted list of all alignment files in the directory ending with .aln
    aln_files = sorted(f for f in os.listdir(aln_dir) if f.endswith(".aln"))

    for fname in aln_files:
    # Extract BUSCO gene ID by removing the '.aln' extension from filename
        busco_id = fname.replace(".aln", "")
        aln_path = os.path.join(aln_dir, fname)

        try:
            with open(aln_path) as f:
                # Read all lines that are not FASTA headers, stripping whitespace
                seqs = [line.strip() for line in f if not is_fasta_header(line)]
                
                # Join sequences together into a single string (concatenated alignment)
                aligned = "".join(seqs)
                
                # Calculate unaligned sequence length by removing gap characters ('-')
                unaligned_length = len(aligned.replace("-", ""))

                # If no sequence content, skip
                if unaligned_length == 0:
                    print(f"Warning: Empty alignment in {fname}, skipping...")
                    continue

                # Define minimum length cutoff as 40% of the unaligned length
                min_len = int(unaligned_length * 0.4)
                # Define maximum length cutoff as 160% of the unaligned length
                max_len = int(unaligned_length * 1.6)

		# Write the BUSCO ID, a zero placeholder, and length cutoffs to the output file
                out.write(f"{busco_id}\t0\t{min_len}\t{max_len}\n")

        except Exception as e:
        # Print error message if the alignment file cannot be processed
            print(f"Error reading {fname}: {e}")
            
            
            
            
