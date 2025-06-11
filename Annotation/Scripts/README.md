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

### 4. `compare.sh`
**Purpose:** Compare both PASA and Braker3 annotations to see which one produced more accurate results. Use BUSCO, gffcompare and AGAT.

### 4. `TEs_analysis.sh`
**Purpose:** Genomic Distribution of Repeats. Identify where TEs are located relative to genes (exonic, intronic, intergenic).

### 4. `merge.sh`
**Purpose:** Pipeline to merge BRAKER and PASA annotations using AGAT and EVM, extract proteins, and assess completeness with BUSCO.
