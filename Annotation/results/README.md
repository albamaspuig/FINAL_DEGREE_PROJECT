## Results for Genome Annotation of *Abeoforma whisleri*

This folder contains the results of the genome annotation analysis, including **structural** and **functional** assessments, as well as **transposable element (TE)** analysis.

### 1. `comparison/`
Contains files related to the comparison between the **previous annotation** (PASA) and the **new annotation** (BRAKER3). Includes:
- **GFFCompare** outputs
- **AGAT** statistics and summaries
- **BUSCO** completeness assessments for both annotations

### 2. `Repeat_analysis/`
Contains results from **transposable element (TE)** analysis, including:
-   **`EDTA/`**
    -   `TE_classification_summary.csv`: Comprehensive summary of TE classification and content.
    -   `EDTA_plots.pdf/`: PDF document containing plots related to TEs characteristics.

-   **`LTRharvest/`**
    -   `ltr_summary.txt`: Summary statistics from LTRharvest.
    -   `ltr_library.fa`: Library of identified LTR retrotransposons.
    -   `ltr_plots.pdf`: PDF document containing plots related to LTR element characteristics.

-   **`RepeatMasker/`**
    -   `Abeoforma_genome_v2.fasta.out`: The main RepeatMasker output file detailing identified repeats.
    -   `Abeoforma_genome_v2.fasta.tbl`: A table summarizing the RepeatMasker results.
    -   `RepeatMasker_plot.png`: Visual representation of the genome's repeat content.

### 3. `eggNOG/`
This directory contains the functional annotation results generated using eggNOG-mapper. 