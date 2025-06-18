# Busco_db Directory

This directory contains the custom BUSCO dataset and associated scripts developed for genome completeness assessment in unicellular holozoans.

## Folder Structure

### `uni_holozoans_odb10/`

This subdirectory holds the custom BUSCO dataset for unicellular holozoans, including all necessary files for running BUSCO with this lineage-specific dataset.

- `dataset.cfg` — Configuration file specifying dataset parameters (name, species, number of BUSCOs, etc.).
- `hmms/` — Folder containing the individual HMM profiles for each BUSCO ortholog group.
- `myholozoa.hmm` and related files (`.h3f`, `.h3i`, `.h3m`, `.h3p`) — Concatenated and compressed HMM database files for efficient searching by BUSCO.
- `Orthogroups_SingleCopyOrthologues.txt` — List of single-copy ortholog groups used in this dataset.
- `id_mapping.tsv` — Mapping file linking original orthogroup IDs to numeric IDs used in the dataset.
- `lengths_cutoff` — File containing length cutoffs per BUSCO group, used for filtering hits.
- `scores_cutoff` — File specifying bit score cutoffs per BUSCO group for reliable detection.

### Species Used for OrthoFinder Analysis

| #  | Species                                      | Lineage           | Source | Number of Proteins |
|----|----------------------------------------------|-------------------|--------|--------------------|
| 1  | *Capsaspora owczarzaki*                      | Filozoa           | MCG    | 10,186             |
| 2  | *Chromosphaera perkinsii*                    | Corallochytrea    | MCG    | 8,500              |
| 3  | *Corallochytrium limacisporum* (Hawaii)      | Corallochytrea    | MCG    | 7,535              |
| 4  | *Corallochytrium limacisporum* (India)       | Corallochytrea    | MCG    | 7,667              |
| 5  | *Creolimax fragrantissima*                   | Ichthyosporea     | NCBI   | 8,694              |
| 6  | *Ichthyophonus hoferi*                       | Ichthyosporea     | NCBI   | 6,351              |
| 7  | *Ministeria vibrans*                         | Filasterea        | MCG    | 12,742             |
| 8  | *Monosiga brevicollis*                       | Choanoflagellata  | NCBI   | 9,203              |
| 9  | *Pigoraptor chileana*                        | Aphelidea         | MCG    | 14,618             |
| 10 | *Pigoraptor vietnamica*                      | Aphelidea         | MCG    | 14,922             |
| 11 | *Salpingoeca rosetta*                        | Choanoflagellata  | NCBI   | 11,731             |
| 12 | *Sphaeroforma arctica*                       | Ichthyosporea     | NCBI   | 18,730             |
| 13 | *Sphaerothecum destruens*                    | Ichthyosporea     | NCBI   | 15,930             |

**Table:** Species used for the construction of the *Unicellular Holozoans* custom BUSCO dataset and OrthoFinder analysis, showing their taxonomic lineage, data source, and number of predicted proteins.

### `scripts/`

### 1. `mapping.py`  
**Purpose:** Maps BUSCO gene IDs (e.g., `HOG0001234`) to numeric identifiers (`1`, `2`, `3`, ...), as required by BUSCO naming conventions.  
**Output:** `id_mapping.tsv`

### 2. `check_busco_length.py`  
**Purpose:** Inspects the distribution of unaligned sequence lengths per BUSCO gene. Prints the minimum, maximum, and total number of sequences per alignment to help validate cutoff ranges.  
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

## Example Usage

To run BUSCO using this custom unicellular holozoans dataset on your input proteome:

```bash
busco -i your_proteome.faa -l path/to/Busco_db/uni_holozoans_odb10 -o busco_output -m prot --offline -f
