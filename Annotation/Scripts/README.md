## Scripts for Genome Annotation of *Abeoforma Whisleri*

This folder contains scripts for genome annotation using BRAKER3 with evidence from RNA-seq, protein homology, and repeat masking.

### 1. `braker3_annotation.sh`
**Purpose:** Pipeline for the full annotation including repeat modeling, masking, RNA-seq alignment, Portcullis hints, and BRAKER3 execution.

### 2. `rna_aligning.sh`
**Purpose:** Aligns RNA-seq reads using HISAT2, sorts and merges alignments, and prepares input files for Portcullis.

### 3. `RNA_evidences.sh`
**Purpose:** Downloads and processes external SRA RNA-seq datasets, performs trimming and alignment, and integrates results into the merged BAM for annotation.

### 4. `convert_portcullis.py`
**Purpose:** Converts BED intron junctions output from Portcullis into GFF format required by BRAKER3.

### 5. `compare.sh`
**Purpose:** Compare both PASA and Braker3 annotations to see which one produced more accurate results. Use BUSCO, gffcompare and AGAT.

### 6. `TEs_analysis.sh`
**Purpose:** Genomic Distribution of Repeats. Identify where TEs are located relative to genes (exonic, intronic, intergenic).

### 7. `merge.sh`
**Purpose:** Pipeline to merge BRAKER and PASA annotations using AGAT and EVM, extract proteins, and assess completeness with BUSCO.

<img src="/images/merging_workflow.png" alt="braker3 Pipeline" width="500"/>

*Diagram illustrating the merging workflow for the two methods.*

### 8. `Revigo_filter.py`
**Purpose:** This script filters REVIGO output tables to retain only the most representative GO terms based on dispensability. It ensures clean, reproducible preparation of enriched GO terms (BP, CC, MF) for visualization, reporting, and supplementary materials.


