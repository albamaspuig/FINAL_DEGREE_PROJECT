## Scripts for Custom BUSCO Database Creation

This folder contains scripts used to build a custom BUSCO lineage dataset for unicellular holozoans. These tools allow generating a new BUSCO dataset from a set of single-copy orthologs obtained via OrthoFinder.

### 1. `mapping.py`
**Purpose:** Maps BUSCO gene IDs (e.g., `HOG0001234`) to numeric identifiers (`1`, `2`, `3`, ...), as required by BUSCO naming conventions.  
**Output:** `id_mapping.tsv`

### 2. `check_busco_length.py`
**Purpose:** To inspect the distribution of unaligned sequence lengths per BUSCO gene. Prints the minimum, maximum, and total number of sequences per alignment to help validate cutoff ranges.
**Output:** Console output summarizing length stats per BUSCO gene.

### 3. `lengthcutoff.py`
**Purpose:** Calculates the acceptable aligned sequence length range for each BUSCO gene. Uses ±20% of the ungapped alignment length to define the valid length window.  
**Output:** `lengths_cutoff`

### 4. `scorescutoff.sh`
**Purpose:** Uses `hmmsearch` on reference BUSCO sequences to compute the minimum score for each HMM profile to be considered valid.  
**Output:** `scores_cutoff`

### 5. `BUSCO.sh`
**Purpose:** Complete, reproducible pipeline to:
- Run OrthoFinder to detect single-copy orthologs
- Align and build HMM profiles
- Map BUSCO IDs to integers
- Concatenate and index HMMs
- Generate length and score cutoffs
- Validate the resulting lineage with `hmmsearch` and `busco`
---
