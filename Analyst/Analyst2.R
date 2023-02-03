# Author: Allison Choy
# BF528 Project 3
# Role: Analyst

# This script is meant to run a concordance calculation of DE analysis between
# RNA-Seq and Microarray analysis
# for overall, above-, and below- median subsets of data
# arranged by decreasing fold change values
# Also included are calculations with probes mapped to a rat genome database rat2302.db
# and it subsequent studies to see if there is a difference
# Lastly, the challenge plot is also included at the end

# Set directory
setwd('/projectnb2/bf528/users/tinman_2022/project_3/Analyst')

# Installing (if necessary) and loading libraries
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.14")
  BiocManager::install("AnnotationDbi")
  BiocManager::install("rat2302.db")
}

if (!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
  install.packages("ggpubr")
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(rat2302.db))

# Reading files into dataframes
lef_rna <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Programmer/de_seq_results/AhR_deseq_results1.csv', as.is = TRUE)
flu_rna <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Programmer/de_seq_results/CARPXR_deseq_results1.csv', as.is = TRUE)
ifo_rna <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Programmer/de_seq_results/DNADamage_deseq_results1.csv', as.is = TRUE)

lef_mic <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Analyst/AhR_padj_limma.csv', as.is = TRUE)
flu_mic <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Analyst/CARPXR_padj_limma.csv', as.is = TRUE)
ifo_mic <- read.csv('/projectnb/bf528/users/tinman_2022/project_3/Analyst/DNADamage_padj_limma.csv', as.is = TRUE)

affy_map <- read_csv('/project/bf528/project_3/refseq_affy_map.csv')

#############################################################################################
# Function to make deseq files into filtered tibble for FC and PV
deseq_tib <- function(dataframe){
  as_tibble(dataframe) %>% 
    rename('X' = 'geneID') %>% 
    filter((pvalue < 0.05) & (abs(log2FoldChange) > log2(1.5)))
}

# Function to make limma files into filtered tibble for FC and PV
limma_tib <- function(dataframe){
  as_tibble(dataframe) %>% 
    rename('X' = 'probeid') %>% 
    filter((P.Value < 0.05) & (abs(logFC) > log2(1.5)))
}

#############################################################################################
# pulls number of probeids in ref matrix of genes
affy_limma<-function(limma_set){
  n_lim <- affy_map %>% 
    filter(PROBEID %in% limma_set$probeid) #%>% 
  return(length(unique(n_lim$SYMBOL)))
}

# pulls number of unique genes from deseq gene IDs
affy_deseq<-function(deseq_set){
  n_lim <- affy_map %>% 
    filter(REFSEQ %in% deseq_set$geneID) #%>% 
  return(length(unique(n_lim$SYMBOL)))
}

#############################################################################################
# Obtains list of genes from limma results that overlap with ref matrix
limma_gene<-function(limma_set){
  lim_gene <- affy_map %>% 
    inner_join(limma_set, c('PROBEID' = 'probeid'))
  return(unique(lim_gene$SYMBOL, .keep_all = TRUE))
}

# Obtains list of unique genes from DESeq results that overlap with ref matrix
deseq_gene <- function(deseq_set){
  deseq_gene <- affy_map %>% 
    inner_join(deseq_set, c('REFSEQ' = 'geneID'))
  return(unique(deseq_gene$SYMBOL, .keep_all = TRUE))
}

#############################################################################################
# Find # concordant genes in both methods
#n0 <-length(intersect(lim_gene1, deseq_gene1))

# Calculating Concordance, background of corrected number of 'true' overlaps
# Equation: nx=Nn0−n1n2/n0+N−n1−n2
# to be consistent, n1 always refers to microarray data and n2 RNASeq data

# Total genes = 13079 from paper, can also use total from affy_map
N <- length(unique(affy_map$SYMBOL)) #12748
N_affy <-13079 # from paper

# Concordance with adjusted 'true' overlaps
concordance <- function(N, n0, n1, n2){
  # Adjusted 'true' number of overlaps
  nx = ((N*n0) - n1*n2)/(n0+N-n1-n2)
  conc <- (2*nx)/(n1+n2)
  return(conc)
}

# For example: conc1 <- concordance(N, n0, n1, n2)

#############################################################################################
# Converting each to a tibble, mic for microarray, rna for RNA-seq
# Also filtered by p-value < 0.05 and |log2FC| > 1.5
lef_mic_tib <- limma_tib(lef_mic)
flu_mic_tib <- limma_tib(flu_mic)
ifo_mic_tib <- limma_tib(ifo_mic)

lef_rna_tib <- deseq_tib(lef_rna)
flu_rna_tib <- deseq_tib(flu_rna)
ifo_rna_tib <- deseq_tib(ifo_rna)

# Calculating gene size in mic set (n1) and rna set (n2)
n1_lef <- affy_limma(lef_mic_tib)
n1_flu <- affy_limma(flu_mic_tib)
n1_ifo <- affy_limma(ifo_mic_tib)

n2_lef <- affy_deseq(lef_rna_tib)
n2_flu <- affy_deseq(flu_rna_tib)
n2_ifo <- affy_deseq(ifo_rna_tib)

# Getting the gene list of DE genes for each drug and method
lef_mic_gene <- limma_gene(lef_mic_tib)
flu_mic_gene <- limma_gene(flu_mic_tib)
ifo_mic_gene <- limma_gene(ifo_mic_tib)

lef_rna_gene <- deseq_gene(lef_rna_tib)
flu_rna_gene <- deseq_gene(flu_rna_tib)
ifo_rna_gene <- deseq_gene(ifo_rna_tib)

# Determining the raw number of overlaps between two sets (n0)
n0_lef <- length(intersect(lef_mic_gene, lef_rna_gene))
n0_flu <- length(intersect(flu_mic_gene, flu_rna_gene))
n0_ifo <- length(intersect(ifo_mic_gene, ifo_rna_gene))

n0_all <-rbind(n0_lef, n0_flu, n0_ifo)
write.csv(n0_all, 'counts.csv')

# Concordance determination between methods for all 3 treatments
conc_lef <- concordance(N, n0_lef, n1_lef, n2_lef)
conc_flu <- concordance(N, n0_flu, n1_flu, n2_flu)
conc_ifo <- concordance(N, n0_ifo, n1_ifo, n2_ifo)

#############################################################################################
# Selects genes that were mapped for 8 or more times with probesets
# Only yielded one result from Leflunomide probeset: Sorl1
deseq_gene_test<-affy_map %>% 
  inner_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n()) %>% 
  filter(count >=8)
# Yielded three results from Fluconazole probeset: Sorl1, Sqstm1, Tpm3
deseq_gene_test2<-affy_map %>% 
  inner_join(flu_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n()) %>% 
  filter(count >=8)
# Yielded no results from Ifosfamide probeset
deseq_gene_test3<-affy_map %>% 
  inner_join(ifo_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n()) %>% 
  filter(count >=8)

# One way to check that there are no duplicates
length(lef_mic_gene) == n1_lef # TRUE

# Checks to make sure duplicates are gone
length(lef_rna_gene) == n2_lef # TRUE

#############################################################################################
# Making table of points for plotting
# Pulling concordance data from microarray
mic_conc <- tibble(treatment = c('LEF', 'FLU', 'IFO'),
              DEGs = c(n1_lef, n1_flu, n1_ifo),
              concordancev = c(conc_lef, conc_flu, conc_ifo))

# Plot for microarray concordance
(mic_conc_plot <- ggplot(mic_conc, aes(x = DEGs, y = concordancev)) + 
    # in case want to add colors: fill = c('#619CFF', 'white', '#F8766D'), 
    geom_point(shape = 1, size = 2) + 
    geom_smooth(aes(x = DEGs, y = concordancev), method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') +
    theme_classic() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.36), labels = scales::percent_format(accuracy = 1)) +
    stat_cor(aes(label = paste(..rr.label..)), label.x = 50) +
    geom_label(aes(label = treatment), vjust = 1.5, nudge_x = 0.15) +
    labs(x = 'Treatment Effect\n(Number of DEGs from Microarray)',
         y = 'Concordance of DEG',
         title = 'Microarray Analysis'))

# Pulling concordance data from RNA-Seq
rna_conc <- tibble(treatment = c('LEF', 'FLU', 'IFO'),
                   DEGs = c(n2_lef, n2_flu, n2_ifo),
                   concordancev = c(conc_lef, conc_flu, conc_ifo))

# Plot for rna concordance
(rna_conc_plot <- ggplot(rna_conc, aes(x = DEGs, y = concordancev)) + 
    geom_point(shape = 1, size = 2) + 
    geom_smooth(aes(x = DEGs, y = concordancev), method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') +
    theme_classic() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.36), labels = scales::percent_format(accuracy = 1)) +
    stat_cor(aes(label = paste(..rr.label..)), label.x = 50) +
    geom_label(aes(label = treatment), vjust = 1.5, nudge_x = 0.15) +
    labs(x = 'Treatment Effect\n(Number of DEGs from RNASeq)',
         y = 'Concordance of DEG',
         title = 'RNA-Seq Analysis'))

Concordance <- ggarrange(rna_conc_plot, mic_conc_plot, labels = 'AUTO')
annotate_figure(Concordance, top = text_grob('Overall Concordance vs Number of DE Genes',
                                         size = 16))
ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/Concordance.png')

#############################################################################################
## WARNING: RABBIT HOLE
# A look at the concordance level with the ratdb Bioconductor package
# due to duplications and/or missing matches

# Total number of unique genes from database - not used for calculations, see Report
N_key <- keys(rat2302.db, keytype = 'SYMBOL')
total_N <- select(rat2302.db, key = N_key, columns = 'SYMBOL', keytype = 'SYMBOL')
length(unique(total_N$SYMBOL)) #42380

# Pulling microarray probeids
mic_key <- keys(rat2302.db, keytype = 'PROBEID')
mic_arr <- select(rat2302.db, key = mic_key, columns = c('SYMBOL', 'GENENAME'), keytype = 'PROBEID')

# Redo calculation of subset mapping to this db - mic set
mic_lef_probes <- tibble(mic_arr) %>% filter(PROBEID %in% lef_mic_tib$probeid)
n1_rat_lef <- length(unique(mic_lef_probes$SYMBOL)) #340 vs 288
mic_flu_probes <- tibble(mic_arr) %>% filter(PROBEID %in% flu_mic_tib$probeid)
n1_rat_flu <- length(unique(mic_flu_probes$SYMBOL)) #721 vs 619
mic_ifo_probes <- tibble(mic_arr) %>% filter(PROBEID %in% ifo_mic_tib$probeid)
n1_rat_ifo <- length(unique(mic_ifo_probes$SYMBOL)) #54 vs 46

# Pulling rna dataset probes (REFSEQ)
rna_key <- keys(rat2302.db, keytype = 'REFSEQ')
rna_arr <- select(rat2302.db, key = rna_key, columns = c('SYMBOL', 'GENENAME'), keytype = 'REFSEQ')

# Redo calculation of subset mapping to this db - mic set
rna_lef_probes <- tibble(rna_arr) %>% filter(REFSEQ %in% lef_rna_tib$geneID) 
#%>% group_by(REFSEQ) %>% summarize(count = n()) %>% filter(count >=2) # resulted in 1 at count of 2
n2_rat_lef <- length(unique(rna_lef_probes$SYMBOL)) #1531 vs 1385
rna_flu_probes <- tibble(rna_arr) %>% filter(REFSEQ %in% flu_rna_tib$geneID) 
n2_rat_flu <- length(unique(rna_flu_probes$SYMBOL)) #2516 vs 2282
rna_ifo_probes <- tibble(rna_arr) %>% filter(REFSEQ %in% ifo_rna_tib$geneID) 
n2_rat_ifo <- length(unique(rna_ifo_probes$SYMBOL)) #129 vs 119

# Results in report, test to see if rat db is worth using
lef_rna_join <- tibble(rna_arr) %>% left_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(padj)) %>% group_by(REFSEQ) %>% summarize(count = n()) %>% filter(count >=2) # 0 results
lef_rna_sym0 <- tibble(rna_arr) %>% left_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(padj)) %>% group_by(SYMBOL) %>% summarize(count = n()) %>% filter(count >=2)
lef_rna_affy_join <- tibble(affy_map) %>% left_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(padj)) %>% group_by(REFSEQ) %>% summarize(count = n()) %>% filter(count >=2) # 462
lef_rna_sym <- tibble(affy_map) %>% left_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(padj)) %>% group_by(SYMBOL) %>% summarize(count = n()) %>% filter(count >=2) # 462 same as REFSEQ
lef_rna_sym2 <- tibble(affy_map) %>% left_join(lef_rna_tib, c('REFSEQ' = 'geneID')) %>% 
  filter(!is.na(padj)) %>% group_by(SYMBOL) %>% summarize(count = n()) %>% filter(count >=8) # 2, 1 is NA with count of 150

# recalculating concordance with rat db mapping
n0_rat_lef <- length(intersect(mic_lef_probes$SYMBOL, rna_lef_probes$SYMBOL))
conc_rat_lef <- concordance(N, n0_rat_lef, n1_rat_lef, n2_rat_lef)
n0_rat_flu <- length(intersect(mic_flu_probes$SYMBOL, rna_flu_probes$SYMBOL))
conc_rat_flu <- concordance(N, n0_rat_flu, n1_rat_flu, n2_rat_flu)
n0_rat_ifo <- length(intersect(mic_ifo_probes$SYMBOL, rna_ifo_probes$SYMBOL))
conc_rat_ifo <- concordance(N, n0_rat_ifo, n1_rat_ifo, n2_rat_ifo)

# Plot of concordance vs mic analysis
mic_rat_conc <- tibble(treatment = c('LEF', 'FLU', 'IFO'),
                   rat_DEGs = c(n1_rat_lef, n1_rat_flu, n1_rat_ifo),
                   rat_conc = c(conc_rat_lef, conc_rat_flu, conc_rat_ifo))

mic_overall <- mic_conc %>% left_join(mic_rat_conc, 'treatment')
#write.csv(mic_overall, 'Mic_overall_conc.csv')

(mic_rat_conc_plot <- ggplot(mic_rat_conc, aes(x = rat_DEGs, y = rat_conc)) + 
    # in case want to add colors: fill = c('#619CFF', 'white', '#F8766D'), 
    geom_point(shape = 1, size = 2) + 
    geom_smooth(aes(x = rat_DEGs, y = rat_conc), method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') +
    theme_classic() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.36), labels = scales::percent_format(accuracy = 1)) +
    stat_cor(aes(label = paste(..rr.label..)), label.x = 50) +
    geom_label(aes(label = treatment), vjust = 1.5, nudge_x = 0.15) +
    labs(x = 'Treatment Effect\n(Number of DEGs from Microarray)',
         y = 'Concordance of DEG',
         title = 'Microarray Analysis'))

# Plot of concordance vs rna analysis
rna_rat_conc <- tibble(treatment = c('LEF', 'FLU', 'IFO'),
                   rat_DEGs = c(n2_rat_lef, n2_rat_flu, n2_rat_ifo),
                   rat_conc = c(conc_rat_lef, conc_rat_flu, conc_rat_ifo))

rna_overall <- rna_conc %>% left_join(rna_rat_conc, 'treatment')
#write.csv(rna_overall, 'RNA-seq_overall_conc.csv')

# Plot for rna concordance
(rna_rat_conc_plot <- ggplot(rna_rat_conc, aes(x = rat_DEGs, y = rat_conc)) + 
    geom_point(shape = 1, size = 2) + 
    geom_smooth(aes(x = rat_DEGs, y = rat_conc), method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') +
    theme_classic() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.36), labels = scales::percent_format(accuracy = 1)) +
    stat_cor(aes(label = paste(..rr.label..)), label.x = 50) +
    geom_label(aes(label = treatment), vjust = 1.5, nudge_x = 0.15) +
    labs(x = 'Treatment Effect\n(Number of DEGs from RNASeq)',
         y = 'Concordance of DEG',
         title = 'RNA-Seq Analysis'))

# pulling plots together
Rat_concordance <- ggarrange(rna_rat_conc_plot, mic_rat_conc_plot, labels = 'AUTO')
annotate_figure(Rat_concordance, top = text_grob('Overall Concordance of DE Genes vs. DEG counts with ratDB',
                                         size = 16))
#ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/Concordance_rat.png')

#############################################################################################
# Subdividing DE genes
# RNA-Seq
lef_rna_above <- lef_rna_tib %>% 
  filter(baseMean > median(baseMean))
lef_rna_below <- lef_rna_tib %>% 
  filter(baseMean < median(baseMean))
flu_rna_above <- flu_rna_tib %>% 
  filter(baseMean > median(baseMean))
flu_rna_below <- flu_rna_tib %>% 
  filter(baseMean < median(baseMean))
ifo_rna_above <- ifo_rna_tib %>% 
  filter(baseMean > median(baseMean))
ifo_rna_below <- ifo_rna_tib %>% 
  filter(baseMean < median(baseMean))

# Microarray
lef_mic_above <- lef_mic_tib %>% 
  filter(AveExpr > median (AveExpr))
lef_mic_below <- lef_mic_tib %>% 
  filter(AveExpr < median (AveExpr))
flu_mic_above <- flu_mic_tib %>% 
  filter(AveExpr > median (AveExpr))
flu_mic_below <- flu_mic_tib %>% 
  filter(AveExpr < median (AveExpr))
ifo_mic_above <- ifo_mic_tib %>% 
  filter(AveExpr > median (AveExpr))
ifo_mic_below <- ifo_mic_tib %>% 
  filter(AveExpr < median (AveExpr))

#############################################################################################
# Computing concordance of each subset
# Above and Below Median Concordance function
med_conc <- function(ma_med, deseq_med){
  n_mic <- affy_limma(ma_med)
  n_rna <- affy_deseq(deseq_med)
  gene_rna <- deseq_gene(deseq_med)
  gene_mic <- limma_gene(ma_med)
  n0_med <- length(intersect(gene_rna, gene_mic))
  med_conc <- concordance(N, n0_med, n_mic, n_rna)
  return(med_conc)
}

med_rat_conc <- function(ma_med, deseq_med){
  mic_probes <- tibble(mic_arr) %>% filter(PROBEID %in% ma_med$probeid)
  n1_rat <- length(unique(mic_probes$SYMBOL)) 
  rna_probes <- tibble(rna_arr) %>% filter(REFSEQ %in% deseq_med$geneID) 
  n2_rat <- length(unique(rna_probes$SYMBOL))
  n0_rat <- length(intersect(mic_probes$SYMBOL, rna_probes$SYMBOL))
  conc_rat <- concordance(N, n0_rat, n1_rat, n2_rat)
  return(conc_rat)
}

above_conc_lef <- med_conc(lef_mic_above, lef_rna_above)
above_conc_flu <- med_conc(flu_mic_above, flu_rna_above)
above_conc_ifo <- med_conc(ifo_mic_above, ifo_rna_above)
below_conc_lef <- med_conc(lef_mic_below, lef_rna_below)
below_conc_flu <- med_conc(flu_mic_below, flu_rna_below)
below_conc_ifo <- med_conc(ifo_mic_below, ifo_rna_below)

above_conc_rat_lef <- med_rat_conc(lef_mic_above, lef_rna_above)
above_conc_rat_flu <- med_rat_conc(flu_mic_above, flu_rna_above)
above_conc_rat_ifo <- med_rat_conc(ifo_mic_above, ifo_rna_above)
below_conc_rat_lef <- med_rat_conc(lef_mic_below, lef_rna_below)
below_conc_rat_flu <- med_rat_conc(flu_mic_below, flu_rna_below)
below_conc_rat_ifo <- med_rat_conc(ifo_mic_below, ifo_rna_below)

#############################################################################################
# Pulling all data together into a tibble for plotting
df_lef <- tibble(treatment = rep('Leflunomide'),
              measure = c('Above Median', 'Below Median', 'Overall'),
              Conc = c(above_conc_lef, below_conc_lef, conc_lef))
df_lef$measure <- factor(df_lef$measure)

df_flu <- tibble(treatment = rep('Fluconazole'),
              measure = c('Above Median', 'Below Median', 'Overall'),
              Conc = c(above_conc_flu, below_conc_flu, conc_flu))
df_flu$measure <- factor(df_flu$measure)

df_ifo <- tibble(treatment = rep('Ifosfamide'),
              measure = c('Above Median', 'Below Median', 'Overall'),
              Conc = c(above_conc_ifo, below_conc_ifo, conc_ifo))
df_ifo$measure <- factor(df_ifo$measure)

df_all <- rbind(df_lef, df_flu, df_ifo)
#######################################################
# Making tables for plotting for ratdb mapping
df_rat_lef <- tibble(treatment = rep('Leflunomide'),
                 measure = c('Above Median', 'Below Median', 'Overall'),
                 Conc = c(above_conc_rat_lef, below_conc_rat_lef, conc_rat_lef))
df_rat_lef$measure <- factor(df_rat_lef$measure)

df_rat_flu <- tibble(treatment = rep('Fluconazole'),
                 measure = c('Above Median', 'Below Median', 'Overall'),
                 Conc = c(above_conc_rat_flu, below_conc_rat_flu, conc_rat_flu))
df_rat_flu$measure <- factor(df_rat_flu$measure)

df_rat_ifo <- tibble(treatment = rep('Ifosfamide'),
                 measure = c('Above Median', 'Below Median', 'Overall'),
                 Conc = c(above_conc_rat_ifo, below_conc_rat_ifo, conc_rat_ifo))
df_rat_ifo$measure <- factor(df_rat_ifo$measure)

df_rat_all <- rbind(df_rat_lef, df_rat_flu, df_rat_ifo)

#############################################################################################
# Bar plot for all three measurements

(plot_all <- ggplot(df_all, aes(x = treatment, y = Conc, group = measure)) + 
    geom_col(aes(fill = measure), position = 'dodge') + 
    geom_text(aes(label = paste0(round(Conc*100, digits = 2), '%'), y = Conc), position = position_dodge(0.9), vjust = -0.5) + 
    theme_bw() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.4),labels = scales::percent_format(accuracy = 1)) +
    labs(x = 'Chemical for Treatment',
         y = 'Concordance of DEGs',
         fill = 'Measure',
         title = 'Concordance Measures of Overall, Above-, and Below-Median\nsubsets for each Chemical Treatment') +
    theme(legend.position = 'bottom'))

ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/Concordance_barplot.png')

(plot_rat_all <- ggplot(df_rat_all, aes(x = treatment, y = Conc, group = measure)) + 
    geom_col(aes(fill = measure), position = 'dodge') + 
    geom_text(aes(label = paste0(round(Conc*100, digits = 2), '%'), y = Conc), position = position_dodge(0.9), vjust = -0.5) + 
    theme_bw() +
    # formatted to include range of % and y limits
    scale_y_continuous(limits = c(0.0,0.4),labels = scales::percent_format(accuracy = 1)) +
    labs(x = 'Chemical for Treatment',
         y = 'Concordance of DEGs',
         fill = 'Measure',
         title = 'Concordance Measures of Overall, Above-, and Below-Median\nsubsets from ratDB for each Chemical Treatment') +
    theme(legend.position = 'bottom'))

#ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/Concordance_rat_barplot.png')

#############################################################################################
# Challenge plot
lef_rna_above_FC <- lef_rna_above %>% arrange(desc(log2FoldChange))
lef_mic_above_FC <- lef_mic_above %>% arrange(desc(logFC))
flu_rna_above_FC <- flu_rna_above %>% arrange(desc(log2FoldChange))
flu_mic_above_FC <- flu_mic_above %>% arrange(desc(logFC))
ifo_rna_above_FC <- ifo_rna_above %>% arrange(desc(log2FoldChange))
ifo_mic_above_FC <- ifo_mic_above %>% arrange(desc(logFC))

lef_rna_below_FC <- lef_rna_below %>% arrange(desc(log2FoldChange))
lef_mic_below_FC <- lef_mic_below %>% arrange(desc(logFC))
flu_rna_below_FC <- flu_rna_below %>% arrange(desc(log2FoldChange))
flu_mic_below_FC <- flu_mic_below %>% arrange(desc(logFC))
ifo_rna_below_FC <- ifo_rna_below %>% arrange(desc(log2FoldChange))
ifo_mic_below_FC <- ifo_mic_below %>% arrange(desc(logFC))

challenge <- function(mic_FC, rna_FC, trt){
  count <- max(c(length(mic_FC$probeid), length(rna_FC$geneID)))
  count2 <- min(c(length(mic_FC$probeid), length(rna_FC$geneID)))
  counter = 0
  while (count != 0){
    if (count == max(c(length(mic_FC$probeid), length(rna_FC$geneID)))){
      counter = counter +1
      mic1 <- mic_FC[1:counter,]
      rna1 <- rna_FC[1:counter,]
      first_pt <- med_conc(mic1, rna1)
      ot_conc <- first_pt
      count_list <- counter
      count <- count-1
    }else if (count < min(c(length(mic_FC$probeid), length(rna_FC$geneID)))){
      counter = counter +1
      mic1 <- mic_FC[1:count2,]
      rna1 <- rna_FC[1:counter,]
      next_pt <- med_conc(mic1, rna1)
      ot_conc <- append(ot_conc, next_pt)
      count_list <- append(count_list, counter)
      count <- count -1
    }else{
      counter = counter + 1
      mic1 <- mic_FC[1:counter,]
      rna1 <- rna_FC[1:counter,]
      next_pt <- med_conc(mic1, rna1)
      ot_conc <- append(ot_conc, next_pt)
      count_list <- append(count_list, counter)
      count <- count -1
    }
  }
  test_tib <- tibble(trt,
                     DEGs = count_list,
                     Concordance = ot_conc)
  return(test_tib)
}

lef_conc_above <-challenge(lef_mic_above_FC, lef_rna_above_FC, 'LEF')
flu_conc_above <-challenge(flu_mic_above_FC, flu_rna_above_FC, 'FLU')
ifo_conc_above <-challenge(ifo_mic_above_FC, ifo_rna_above_FC, 'IFO')

lef_conc_below <-challenge(lef_mic_below_FC, lef_rna_below_FC, 'LEF')
flu_conc_below <-challenge(flu_mic_below_FC, flu_rna_below_FC, 'FLU')
ifo_conc_below <-challenge(ifo_mic_below_FC, ifo_rna_below_FC, 'IFO')

conc_above <- rbind(lef_conc_above, flu_conc_above, ifo_conc_above)
conc_below <- rbind(lef_conc_below, flu_conc_below, ifo_conc_below)

conc_above_plot <- ggplot(conc_above, aes(x = DEGs, y = Concordance, group = trt)) +
  geom_line(aes(color = trt)) +
  scale_y_continuous(limits = c(-0.001,0.35),labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  labs(x = 'Number of DEGs ordered by Fold Change',
       y = 'DEG Concordance (%)',
       color = 'Treatment',
       title = 'Above-Median Expression')

conc_below_plot <- ggplot(conc_below, aes(x = DEGs, y = Concordance, group = trt)) +
  geom_line(aes(color = trt)) +
  scale_y_continuous(limits = c(-0.001,0.35),labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  labs(x = 'Number of DEGs ordered by Fold Change',
       y = 'DEG Concordance (%)',
       color = 'Treatment',
       title = 'Below-Median Expression')

Challenge_plot <- ggarrange(conc_above_plot, conc_below_plot, labels = 'AUTO',
                            common.legend = TRUE, legend = 'bottom')
annotate_figure(Challenge_plot, top = text_grob('DEG Abundance Dependent Concordance between RNA-Seq and Microarray',
                                                size = 16))
ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/Challenge_plot.png')
