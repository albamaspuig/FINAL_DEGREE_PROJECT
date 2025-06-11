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

#### Species used for Orthofinder analysis:

| #  | Species Name                         |
|----|------------------------------------|
| 1  | Capsaspora_owczarzaki.faa           |
| 2  | Chromosphaera_perkinsii.faa         |
| 3  | Corallochytrium_limacisporum_HI.faa|
| 4  | Corallochytrium_limacisporum_IN.faa|
| 5  | Creolimax_fragrantissima.faa        |
| 6  | Ichthyophonus_hoferi.faa            |
| 7  | Ministeria_vibrans.faa              |
| 8  | Monosiga_brevicollis.faa            |
| 9  | Pigoraptor_chileana.faa             |
| 10 | Pigoraptor_vietnamica.faa           |
| 11 | Salpingoeca_rosetta.faa             |
| 12 | Sphaeroforma_arctica.faa            |
| 13 | Sphaerothecum_destruens.faa         |

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
busco -i your_proteome.faa -l path/to/Busco_db/uni_holozoans_odb10 -o busco_output -m prot
