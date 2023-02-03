# Project Description

This project attempts to replicate results, specifically Figures 2A and 3B+C from the paper by Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. **“A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data”** _Nature Biotechnology_ 32 (9): 926–32. PMID: 4243706

*Final Report: `BF528_tinman_Project_3.pdf`*

# Contributors
  
1.  Data Curator: Aneeq Husain - @aneeqh
2.  Programmer: Reshma Ramaiah - @ReshmaRamaiah10
3.  Analyst: Allison Choy - @AllisonWChoy
4.  Biologist: Rizky Kafrawi -@rkafrawi

# Repository Contents

## DataCurator

Generates a total of 9 bam files with alignment measures for each sample using the following tools:
- ``STAR.qsub`` - Script for the STAR aligner
- ``multiqc.qsub`` - Script for the multiqc aligner

## Programmer

- ``featureCounts.qsub`` - This script calculates feature counts from bam files
- ``multiqc.qsub`` - This script conducts a multiqc analysis on the feature counts results
- ``combine_csv.R``
  - This script combines featureCounts output files into one csv file
  - Additionally, it creates boxplots of distribution of counts
- ``run_deseq.R`` - This script does a DeSeq2 analysis

## Analyst

- `` Analyst1.R`` 
- Script to generate microarray analysis with limma for tox-group dataset
  - Writes out .csv files for significant genes sorted by adj. p-values and top 10 genes for each treatment group
  - Generates histogram and volcano plots for visualizing fold change and p-values for each treatment analysis

- ``Analyst2.R``
- Generates concordance calculations for overall, above-, and below-median genes for all analyses
- Generates the following plots:
  - concordance against # DE genes for treatment effects for all analyses
  - concordance for all subsets for each chemical treatment
  - challenge plot with increasing number of DE genes for above- and below-median subsets ranked by FC
  
## Biologist

- ``heatmap_project_3_biologist.R``
  - Generates heatmap for part 7.2 in project