print("Converting Portcullis output to GFF...")

import sys

input_bed = sys.argv[1]
output_gff = sys.argv[2]

with open(input_bed, 'r') as infile, open(output_gff, 'w') as outfile:
    for line in infile:
        fields = line.strip().split('\t')
        if len(fields) < 6:
            continue
        chrom = fields[0]
        start = int(fields[1]) + 1  # BED is 0-based, GFF is 1-based
        end = int(fields[2])
        score = fields[4]
        strand = fields[5]

        # Creating the attributes separately
        attributes = f"src=E;pri=4;"

        # Writing the GFF line
        outfile.write(f"{chrom}\tPortcullis\tintron\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}\n")

