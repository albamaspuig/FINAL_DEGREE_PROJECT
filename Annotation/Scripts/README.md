# Annotation Pipeline Scripts

Scripts for genome annotation of *Abeoforma* using BRAKER3 with evidence from RNA-seq, protein homology, and repeat masking.

## Main Pipeline

- `braker3_annotation_pipeline.sh`: Orchestrates the full pipeline: repeat modeling, masking, alignment, Portcullis hints, BRAKER3 execution.

## RNA-Seq Processing

- `rna_aligning.sh`: Aligns RNA-seq reads using HISAT2, sorts, merges, and prepares input for Portcullis.

## External RNA Evidence

- `RNA_evidences.sh`: Downloads and processes external SRA datasets, performs trimming, alignment, and integrates into the merged BAM.

## Utilities

- `convert_portcullis.py`: Converts BED intron junctions from Portcullis to GFF format (required for BRAKER3).
