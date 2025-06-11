#!/bin/bash

source ~/miniconda3/bin/activate

# Functional annotation
# 1. Quality Assessment 
# Run BUSCO

conda activate /home/shared/envs/busco

#Run BUSCO in proteome mode
# PASA busco
busco -i /mnt/Franklin/amas/Abeoforma_PASA.agat.pep.fasta -l /home/amas/annotation/eukaryota_odb10 -o /mnt/Franklin/amas/busco_pasa_proteins -m protein -f --offline
#92.9%

# braker3 output without published RNA or PASA proteomes as evidences
busco -i /home/amas/annotation/braker_output8/braker.aa -l /home/amas/annotation/eukaryota_odb10 -o /mnt/Franklin/amas/busco_braker_proteins3 -m protein -f --offline
#72.5% --> missing 66

# braker3 output with all evidences: published RNA and PASA proteomes
busco -i /home/amas/annotation/braker_output9/braker.aa -l /home/amas/annotation/eukaryota_odb10 -o /home/amas/annotation/busco_braker_proteins9 -m protein -f --offline
#80.00%

#BUSCO on the assembly
busco -i /mnt/Franklin/amas/Abeoforma_genome_v2.fasta -l /home/amas/annotation/eukaryota_odb10 -o /mnt/Franklin/amas/busco_genome_assembly -m genome -f --offline
#84%

#Generate plots
cp /home/amas/annotation/home/amas/annotation/busco_braker_proteins9/short_summary.* ~/annotation/BUSCO_summaries/
cp /mnt/Franklin/amas/busco_pasa_proteins/short_summary.* ~/annotation/BUSCO_summaries/

#BUSCO plot generation tool.

python3 /home/shared/envs/busco/bin/generate_plot.py -wd ~/annotation/BUSCO_summaries

conda deactivate

# Structural Comparison of Annotations
# gff compare --> how many genes and transcripts from BRAKER overlap or are supported by PASA,
#https://github.com/gpertea/gffcompare
#gffcompare -r <reference.gtf/gff> -o <output_prefix> <query.gtf/gff>

conda activate /home/shared/envs/TEs/

#using PASA as reference
gffcompare -r /mnt/Franklin/amas/Abeoforma_PASA.agat_v2.renamed.gff -o /mnt/Franklin/amas/pasa_vs_braker3 /home/amas/annotation/braker_output9/braker.gff3 

cat /mnt/Franklin/amas/pasa_vs_braker3.stats
'''
#= Summary for dataset: /home/amas/annotation/braker_output9/braker.gff3 
#     Query mRNAs :   10005 in    9746 loci  (7839 multi-exon transcripts)
#            (247 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   21477 in   18808 loci  (17355 multi-exon)
# Super-loci w/ reference transcripts:     3636
#-----------------| Sensitivity | Precision  |
        Base level:    17.0     |    37.8    |
        Exon level:    18.4     |    34.5    |
      Intron level:    20.8     |    38.9    |
Intron chain level:     8.3     |    18.4    |
  Transcript level:     6.7     |    14.4    |
       Locus level:     7.6     |    14.6    |

     Matching intron chains:    1443
       Matching transcripts:    1443
              Matching loci:    1423

          Missed exons:  123188/156115	( 78.9%)
           Novel exons:   50766/83199	( 61.0%)
        Missed introns:  107776/137191	( 78.6%)
         Novel introns:   44581/73410	( 60.7%)
           Missed loci:   15122/18808	( 80.4%)
            Novel loci:    6078/9746	( 62.4%)

'''

# Gene Feature Summary (AGAT)
#Quick overview of features in both annotations — which predicts more genes, exons per transcript, UTRs (if present), etc.

conda deactivate
conda activate /home/shared/envs/agat

agat_sp_statistics.pl --gff /mnt/Franklin/amas/Abeoforma_PASA.agat_v2.renamed.gff -o AGAT/pasa_stats.tsv
agat_sp_statistics.pl --gff /home/amas/annotation/braker_output9/braker.gff3  -o AGAT/braker3_stats.tsv

agat_sp_compare_two_annotations.pl -gff1 /home/amas/annotation/braker_output9/braker.gff3 -gff2 /mnt/Franklin/amas/Abeoforma_PASA.agat_v2.renamed.gff -o AGAT/braker3_vs_pasa_compare.tsv

'''
----------------------------------------------------------------------------------------------
|                                       gene with mrna                                       |
----------------------------------------------------------------------------------------------
|            braker            |Abeoforma_PASA.agat_v2.renamed|        Number of cases       |
----------------------------------------------------------------------------------------------
|              0               |              1               |            15197             |
|              1               |              0               |             6089             |
|              1               |              1               |             3508             |
|              1               |              2               |              72              |
|              1               |              3               |              11              |
|              1               |              4               |              6               |
|              1               |              5               |              3               |
|              2               |              0               |              1               |
|              2               |              1               |              25              |
|              2               |              2               |              1               |
|              2               |              3               |              1               |
|              2               |              5               |              1               |
----------------------------------------------------------------------------------------------
Number gene in braker: 9747
Number gene in Abeoforma_PASA.agat_v2.renamed: 18956
'''

'''
~802% of PASA genes (15,197/18,956) are unique to PASA → indicates substantial additional content captured by the transcriptome-based assembly.


~63% of BRAKER genes (6,089/9,747) are unique to BRAKER, suggesting that BRAKER may capture some ab initio-predicted genes not well-represented in the transcriptome assembly (lowly expressed or condition-specific genes).


~36% of BRAKER genes (3,508/9,747) directly overlap with PASA genes, indicating that many of the high-confidence BRAKER genes are already captured by PASA.
'''

#Run TE_analysis for each annotation to compare where TEs are annotated.

