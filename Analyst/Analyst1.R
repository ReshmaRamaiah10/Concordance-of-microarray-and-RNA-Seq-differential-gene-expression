# Author: Allison Choy
# BF528 Project 3
# Role: Analyst

# This script is meant to run limma to determine DEGs for a dataset in microarray analysis
# Then it will order the results by padj and write out the files of this along with the top 10 genes
# ranked by p-values
# This will also plot the histogram as well as a volcano plot of the significant genes for each treatment

# Set directory
setwd('/projectnb/bf528/users/tinman_2022/project_3/Analyst')

# Installing (if necessary) and setting up libraries
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.14")
  BiocManager::install("limma")
}
if (!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
  install.packages("ggpubr")
}

# Libraries needed for this script
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggpubr))

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_3_mic_info.csv',
                    as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

# subset the full expression matrix to just those in this comparison
lef.rma.subset <- rma[paste0('X',samples$array_id[(samples$chemical == 'LEFLUNOMIDE')|(samples$chemical == 'Control' & samples$vehicle == 'CORN_OIL_100_%')])]
flu.rma.subset <- rma[paste0('X',samples$array_id[(samples$chemical == 'FLUCONAZOLE')|(samples$chemical == 'Control' & samples$vehicle == 'CORN_OIL_100_%')])]
ifo.rma.subset <- rma[paste0('X',samples$array_id[(samples$chemical == 'IFOSFAMIDE')|(samples$chemical == 'Control' & samples$vehicle == 'SALINE_100_%')])]

affy_map <- read_csv('/project/bf528/project_3/refseq_affy_map.csv')

#############################################################################################
## LEFLUNOMIDE
# construct a design matrix modeling treatment vs control for use by limma
lef_design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle == 'CORN_OIL_100_%'],
    levels=c('Control','LEFLUNOMIDE')
  )
)
colnames(lef_design) <- c('Intercept','LEFLUNOMIDE')

# run limma
lef_fit <- lmFit(lef.rma.subset, lef_design)
lef_fit <- eBayes(lef_fit)
lef_t <- topTable(lef_fit, coef=2, n=nrow(lef.rma.subset), adjust='BH')

# Filtered results according to the padj values (< 0.05)
sorted_lef <- lef_t[order(lef_t$adj.P.Val),]
filtered_lef <- sorted_lef %>% 
  filter(adj.P.Val<0.05) 

# Pulling top 10 DE genes by adjusted p-value
top_10_lef_p <- resort_lef[1:10,]

# Pulling out the top 10 DE genes by p-value for this set
resort_lef <- lef_t[order(lef_t$P.Value),]
top_10_lef <- as_tibble(resort_lef, rownames = 'probeid') %>% 
  left_join(affy_map, c('probeid' = 'PROBEID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  dplyr::select(probeid, SYMBOL, GENENAME, logFC, P.Value, adj.P.Val) %>% 
  head(10)

# Counts the total number of DE genes from their respective analysis
lef_DEG <- filtered_lef %>% summarize(lef_DEGs = n())
# 466 DEGs

#############################################################################################
## FLUCONAZOLE
# construct a design matrix modeling treatment vs control for use by limma
flu_design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle == 'CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(flu_design) <- c('Intercept','FLUCONAZOLE')

# run limma
flu_fit <- lmFit(flu.rma.subset, flu_design)
flu_fit <- eBayes(flu_fit)
flu_t <- topTable(flu_fit, coef=2, n=nrow(flu.rma.subset), adjust='BH')

# Filtered results according to the padj values (< 0.05)
sorted_flu <- flu_t[order(flu_t$adj.P.Val),]
filtered_flu <- sorted_flu %>% 
  filter(adj.P.Val<0.05) 

# Pulling top 10 DE genes by adjusted p-value
top_10_flu_p <- resort_flu[1:10,]

# Pulling out the top 10 DE genes by p-value for this set
resort_flu <- flu_t[order(flu_t$P.Value),]
top_10_flu <- as_tibble(resort_flu, rownames = 'probeid') %>% 
  left_join(affy_map, c('probeid' = 'PROBEID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  #filter(base::unique(SYMBOL)) %>% 
  dplyr::select(probeid, SYMBOL, GENENAME, logFC, P.Value, adj.P.Val) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  distinct(probeid, .keep_all = TRUE) %>% 
  head(10)

# Counts the total number of DE genes from their respective analysis
flu_DEG <-filtered_flu %>% summarize(flu_DEGs = n())
# 1997 DEGs

#############################################################################################
## IFOSFAMIDE
# construct a design matrix modeling treatment vs control for use by limma
ifo_design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle == 'SALINE_100_%'],
    levels=c('Control','IFOSFAMIDE')
  )
)
colnames(ifo_design) <- c('Intercept','IFOSFAMIDE')

# run limma
ifo_fit <- lmFit(ifo.rma.subset, ifo_design)
ifo_fit <- eBayes(ifo_fit)
ifo_t <- topTable(ifo_fit, coef=2, n=nrow(ifo.rma.subset), adjust='BH')

# Fitered results according to the padj values (< 0.05)
sorted_ifo <- ifo_t[order(ifo_t$adj.P.Val),]
filtered_ifo <- sorted_ifo %>% 
  filter(adj.P.Val<0.05)

# Pulling top 10 DE genes by adjusted p-value
top_10_ifo_p <- sorted_ifo[1:10,]

# Pulling out the top 10 DE genes by p-value for this set
resort_ifo <- ifo_t[order(ifo_t$P.Value),]
top_10_ifo <- as_tibble(resort_ifo, rownames = 'probeid') %>% 
  left_join(affy_map, c('probeid' = 'PROBEID')) %>% 
  filter(!is.na(SYMBOL)) %>% 
  dplyr::select(probeid, SYMBOL, GENENAME, logFC, P.Value, adj.P.Val) %>% 
  head(10)

# Counts the total number of DE genes from their respective analysis
ifo_DEG<-filtered_ifo %>% summarize(ifo_DEGs = n())
# 0 DEGs

#############################################################################################
# write out the results to file
write.csv(sorted_lef,'AhR_padj_limma.csv')
write.csv(sorted_flu,'CARPXR_padj_limma.csv')
write.csv(sorted_ifo,'DNADamage_padj_limma.csv')
write.csv(top_10_lef,'top_10_AHR_limma.csv')
write.csv(top_10_flu,'top_10_CARPXR_limma.csv')
write.csv(top_10_ifo,'top_10_DNADamage_limma.csv')

# Combine all DEGs
overall_DEG <- rbind(c(lef_DEG, flu_DEG, ifo_DEG))
write.csv(overall_DEG,'No_of_DEG_limma.csv')

# Combining all top 10 DEGs
top_10_all <- cbind('AhR probeid' = rownames(top_10_lef),
                    #'AhR Log2FC' = top_10_lef$LogFC,
                    'AhR P.Value' = top_10_lef$P.Value,
                    'AhR P.Adj.' = top_10_lef$adj.P.Val,
                    'CAR/PXR probeid' = rownames(top_10_flu),
                    'CAR/PXR P.Value' = top_10_flu$P.Value,
                    'CAR/PXR P.Adj.' = top_10_flu$adj.P.Val,
                    'DNA Damage probeid' = rownames(top_10_ifo),
                    'DNA Damage P.Value' = top_10_ifo$P.Value,
                    'DNA Damage P.Adj.' = top_10_ifo$adj.P.Val
                     )
write.csv(top_10_all,'top_10_diffexp_limma.csv')

#############################################################################################
# Histogram
lef_sig <- sorted_lef[sorted_lef$P.Value<0.05&(abs(sorted_lef$logFC) > log2(1.5)),]
lef_hist <- lef_sig %>% 
  ggplot(aes(x = logFC)) +
  geom_histogram(bins = 50, fill = 'light blue', color = 'black') +
  theme_classic() +
  labs(x = 'log2 Fold Change',
       y = 'Counts for Log2 Fold Change',
       title = 'Leflunomide - MOA: AhR')

flu_sig <- sorted_flu[sorted_flu$P.Value<0.05&(abs(sorted_flu$logFC) > log2(1.5)),]
flu_hist <- flu_sig %>% 
  ggplot(aes(x = logFC)) +
  geom_histogram(bins = 50, fill = 'light blue', color = 'black') +
  theme_classic() +
  labs(x = 'log2 Fold Change',
       y = 'Counts for Log2 Fold Change',
       title = 'Fluconazole - MOA: CAR/PXR')

ifo_sig <- sorted_ifo[sorted_ifo$P.Value<0.05&(abs(sorted_ifo$logFC) > log2(1.5)),]
ifo_hist <- ifo_sig %>% 
  ggplot(aes(x = logFC)) +
  geom_histogram(bins = 50, fill = 'light blue', color = 'black') +
  theme_classic() +
  labs(x = 'log2 Fold Change',
       y = 'Counts for Log2 Fold Change',
       title = 'Ifosfamide - MOA: DNA Damage')

histo <- ggarrange(lef_hist, flu_hist, ifo_hist, nrow = 1, ncol = 3, labels = 'AUTO')
annotate_figure(histo, top = text_grob('Fold Change Values of Significant DE Genes for Microarray Analysis',
                                         size = 16))
ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/toxogroup3_histogram.png')

#############################################################################################
# Scatter Plot
# AhR - Leflunomide

lef_scatter_sig <- as_tibble(sorted_lef, rownames = 'probeid') %>% 
  mutate(top_sig = case_when(probeid %in% top_10_lef$probeid ~ 'Top 10 P.Val',
                            probeid %in% rownames(lef_sig) ~ 'P-Val & Log2FC',
                            TRUE ~ 'NS'))

lef_scatter_sig$top_sig <-factor(lef_scatter_sig$top_sig, levels = c('Top 10 P.Val', 'P-Val & Log2FC', 'NS'))

lef_scatter <- lef_scatter_sig %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(aes(color = top_sig), alpha = 0.8) +
  scale_color_manual(values = c('Top 10 P.Val' = '#F8766D', 'P-Val & Log2FC' = '#619CFF', 'NS' = 'black')) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), color = 'red') +
  labs(x = 'Log2 Fold Change',
       y = '-Log10 of P Values',
       color = 'Significance By:',
       title = 'Leflunomide - MOA: AhR') +
  theme(legend.position = 'bottom')
lef_scatter

#############################################################################################
# Fluconazole

flu_scatter_sig <- as_tibble(sorted_flu, rownames = 'probeid') %>% 
  mutate(top_sig = case_when(probeid %in% top_10_flu$probeid ~ 'Top 10 P.Val',
                            probeid %in% rownames(flu_sig) ~ 'P-Val & Log2FC',
                            TRUE ~ 'NS'))

flu_scatter_sig$top_sig <-factor(flu_scatter_sig$top_sig, levels = c('Top 10 P.Val', 'P-Val & Log2FC', 'NS'))

flu_scatter <- flu_scatter_sig %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(aes(color = top_sig), alpha = 0.8) +
  scale_color_manual(values = c('Top 10 P.Val' = '#F8766D', 'P-Val & Log2FC' = '#619CFF', 'NS' = 'black')) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), color = 'red') +
  labs(x = 'Log2 Fold Change',
       y = '-Log10 of P Values',
       color = 'Significance By:',
       title = 'Fluconazole - MOA: CAR/PXR') +
  theme(legend.position = 'bottom')
flu_scatter

#############################################################################################
# Ifosfamide

ifo_scatter_sig <- as_tibble(sorted_ifo, rownames = 'probeid') %>% 
  mutate(top_sig = case_when(probeid %in% top_10_ifo$probeid ~ 'Top 10 P.Val',
                            probeid %in% rownames(ifo_sig) ~ 'P-Val & Log2FC',
                            TRUE ~ 'NS'))

ifo_scatter_sig$top_sig <-factor(ifo_scatter_sig$top_sig, levels = c('Top 10 P.Val', 'P-Val & Log2FC', 'NS'))

ifo_scatter <- ifo_scatter_sig %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(aes(color = top_sig), alpha = 0.8) +
  scale_color_manual(values = c('Top 10 P.Val' = '#F8766D', 'P-Val & Log2FC' = '#619CFF', 'NS' = 'black')) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), color = 'red') +
  labs(x = 'Log2 Fold Change',
       y = '-Log10 of P Values',
       color = 'Significance By:',
       title = 'Ifosfamide - MOA: DNA Damage') +
  theme(legend.position = 'bottom')
ifo_scatter

#############################################################################################
# Arranging Plots and saving them
volcano <- ggarrange(lef_scatter, flu_scatter, ifo_scatter, nrow = 1, ncol = 3, labels = 'AUTO',
          common.legend = TRUE, legend = 'bottom')
annotate_figure(volcano, top = text_grob('Fold Change vs Nominal P-Values of Significant DE Genes',
                                         size = 16))
ggsave('/projectnb/bf528/users/tinman_2022/project_3/Analyst/lef_flu_ifo_scatter.png')
