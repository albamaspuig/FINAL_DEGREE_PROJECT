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
# Reference mRNAs :   21477 in   18956 loci  (17355 multi-exon)
# Super-loci w/ reference transcripts:     9455
#-----------------| Sensitivity | Precision  |
        Base level:    43.8     |    97.5    |
        Exon level:    46.6     |    87.4    |
      Intron level:    52.7     |    98.4    |
Intron chain level:    21.4     |    47.4    |
  Transcript level:    17.3     |    37.2    |
       Locus level:    19.4     |    37.7    |

     Matching intron chains:    3717
       Matching transcripts:    3717
              Matching loci:    3673

          Missed exons:   72937/156115  ( 46.7%)
           Novel exons:     884/83199   (  1.1%)
        Missed introns:   63702/137191  ( 46.4%)
         Novel introns:     484/73410   (  0.7%)
           Missed loci:    9380/18956   ( 49.5%)
            Novel loci:     212/9746    (  2.2%)

 Total union super-loci across all input datasets: 9667
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
|            braker            |    Abeoforma_PASA.agat_v2    |        Number of cases       |
----------------------------------------------------------------------------------------------
|              0               |              0               |              2               |
|              0               |              1               |             9390             |
|              1               |              0               |              233             |
|              1               |              1               |             9262             |
|              1               |              2               |              110             |
|              1               |              3               |              2               |
|              2               |              1               |              61              |
|              2               |              2               |              5               |
|              2               |              3               |              1               |
|              3               |              1               |              2               |
----------------------------------------------------------------------------------------------
Number gene in braker: 9747
Number gene in Abeoforma_PASA.agat_v2: 18954

'''



