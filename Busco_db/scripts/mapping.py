#!/usr/bin/env python3

import os

input_dir = "/mnt/Franklin/amas/busco_input/uni_holozoans_odb10/hmms"
id_map_file = "/mnt/Franklin/amas/busco_input/uni_holozoans_odb10/id_mapping.tsv"

hog_ids = sorted([f.split(".")[0] for f in os.listdir(input_dir) if f.endswith(".hmm")])
with open(id_map_file, "w") as out:
    for i, hog in enumerate(hog_ids, 1):
        out.write(f"{i}\t{hog}\n")
