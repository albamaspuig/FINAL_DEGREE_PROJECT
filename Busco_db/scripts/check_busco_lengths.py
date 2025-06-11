import os
from Bio import AlignIO

aln_dir = "/mnt/Franklin/amas/busco_input/uni_holozoans_odb10/hmms"

for fname in os.listdir(aln_dir):
    if fname.endswith(".aln"):
        busco_id = fname.replace(".aln", "")
        aln_path = os.path.join(aln_dir, fname)
        alignment = AlignIO.read(aln_path, "fasta")
        lengths = [len(str(record.seq).replace("-", "")) for record in alignment]
        max_len = max(lengths)
        min_len = min(lengths)
        print(f"BUSCO {busco_id}: min={min_len}, max={max_len}, n={len(lengths)}")

