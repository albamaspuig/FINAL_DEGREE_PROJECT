#!/usr/bin/env python3
import os

hmm_dir = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/hmms"
output_file = "/mnt/Franklin/amas/busco_input/my_busco_dataset_odb10/scores_cutoff"

default_score = 50  # You can adjust this

with open(output_file, "w") as out:
    for fname in os.listdir(hmm_dir):
        if fname.endswith(".hmm"):
           busco_id = fname.replace(".hmm", "")  # Example: N0.HOG0005658
           out.write(f"{busco_id}\t{default_score}\n") 
