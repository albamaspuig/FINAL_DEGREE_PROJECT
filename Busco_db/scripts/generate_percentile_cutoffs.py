import os
import numpy as np
from Bio import AlignIO

aln_dir = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/hmms"
output_file = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/lengths_cutoff"

with open(output_file, "w") as out:
    for fname in os.listdir(aln_dir):
        if fname.endswith(".aln"):
            busco_id = fname.replace(".aln", "")
            aln_path = os.path.join(aln_dir, fname)
            try:
                alignment = AlignIO.read(aln_path, "fasta")
                lengths = [len(str(rec.seq).replace("-", "")) for rec in alignment]
                min_len = int(np.percentile(lengths, 10))
                max_len = int(np.percentile(lengths, 90))
                out.write(f"{busco_id}\t0\t{min_len}\t{max_len}\n")
            except Exception as e:
                print(f"Failed to process {fname}: {e}")

