#!/usr/bin/env python3
import os

aln_dir = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/hmms"
output_file = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/lengths_cutoff"

def is_fasta_header(line):
    return line.startswith(">")

with open(output_file, "w") as out:
    aln_files = sorted(f for f in os.listdir(aln_dir) if f.endswith(".aln"))

    for fname in aln_files:
        busco_id = fname.replace(".aln", "")
        aln_path = os.path.join(aln_dir, fname)

        try:
            with open(aln_path) as f:
                seqs = [line.strip() for line in f if not is_fasta_header(line)]
                aligned = "".join(seqs)
                unaligned_length = len(aligned.replace("-", ""))

                # If no sequence content, skip
                if unaligned_length == 0:
                    print(f"Warning: Empty alignment in {fname}, skipping...")
                    continue

                min_len = int(unaligned_length * 0.4)
                max_len = int(unaligned_length * 1.6)

                out.write(f"{busco_id}\t0\t{min_len}\t{max_len}\n")

        except Exception as e:
            print(f"Error reading {fname}: {e}")
