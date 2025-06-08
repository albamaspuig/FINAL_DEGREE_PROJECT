# Final Degree Project: Genome Annotation Pipeline
This project implements a pipeline for annotating the genome of *Abeoforma*. The pipeline includes 
RNA-seq alignment, genome masking, and gene annotation using tools like HISAT2, Portcullis, and 
BRAKER...

This repository includes:
1. A structural and functional reannotation pipeline for *Abeoforma whisleri*
2. The creation of a custom BUSCO database for unicellular holozoans

### 1. annotation/
Contains all scripts and results related to genome annotation.

- `scripts/`: Annotation pipeline scripts, including genome annotation, RNA-seq processing, and external evidence integration
- `data/`: Input RNAseq and protein data
- `results/`: BRAKER output (with and without external evidence)

### 2. Busco_db/
Contains scripts and files to create a custom BUSCO database.

- `scripts/`: Scripts to process scores and generate dataset files
- `my_busco_dataset_odb10/`: Final BUSCO dataset
- `busco_input/`: Input orthologs used to create the DB
