#Programmer: Reshma Ramaiah
# This script does a DeSeq2 analysis

# To install packages switch to TRUE: 
if (FALSE){
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

# load libraries
library(tidyverse)
library(DESeq2)
library(stringr)

setwd("/projectnb2/bf528/users/tinman_2022/project_3/Programmer/de_seq_results")

# read control and treatment counts
control_counts <- read_csv("/project/bf528/project_3/samples/control_counts.csv")
control_counts <- control_counts %>% select(Geneid, SRR1178050, SRR1178061, SRR1178063, SRR1178004, SRR1178006, SRR1178013)
treatment_counts <- read_csv("/projectnb2/bf528/users/tinman_2022/project_3/Programmer/featureCounts_output/featureCounts_combined.csv")
# combine counts files
combined_counts <- inner_join(control_counts, treatment_counts, by="Geneid")
combined_counts2 <- subset(combined_counts,rowSums(combined_counts==0)==0)
combined_counts2 <- column_to_rownames(combined_counts2, var="Geneid")
# reading metadata
meta <- read_csv("/project/bf528/project_3/toxgroups/toxgroup_3_rna_info.csv")

control_list <- meta$Run[meta$mode_of_action=="Control"] 
exp_list <- meta$Run[meta$mode_of_action!="Control"]
full_list <- meta$Run
# mode of actions
meta$mode_of_action <- factor(meta$mode_of_action)
modes_of_action <- unique(meta$mode_of_action)
modes_of_action
exp_modes_of_action <- modes_of_action[modes_of_action!="Control"]
exp_modes_of_action

# Initialize list to store deseq results:
deseq_list <- list()  # To store deseq objects
deseq_tib_list  <- list()  # To store results converted to tibble

# Run DESeq for each group
for (exp_group in exp_modes_of_action){
  print("Experimental group:")
  print(exp_group)
  # Get the experimental group sample names:
  temp_exp_sample_list <- meta$Run[meta$mode_of_action==exp_group]
  #Grab list of applicable control samples (based on vehicle value):
  exp_vehicle <- meta$vehicle[meta$Run==temp_exp_sample_list][1]
  print("exp_vehicle:")
  print(exp_vehicle)
  control_list2 <- meta$Run[meta$vehicle==exp_vehicle & meta$mode_of_action=="Control"]
  print("Control list:")
  print(control_list2)
  
  # Subset dataframe:
  temp_tib <- combined_counts2[,c(temp_exp_sample_list, control_list2)]
  #Subset meta data:
  temp_meta <- meta[meta$mode_of_action==exp_group | (meta$mode_of_action=="Control" & meta$vehicle==exp_vehicle),]

  # Check order matches:
  does_it_match <- length(names(temp_tib)) == sum(names(temp_tib)==temp_meta$Run)
  if(does_it_match){
    # create the DESeq object
    dds1 <- DESeqDataSetFromMatrix(
      countData = temp_tib,
      colData = temp_meta,
      design= ~ mode_of_action
    )
    
    # relevel mode_of_action as factor
    dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')
    
    # run DESeq
    dds1 <- DESeq(dds1)
    # Adjust for our actual samples:
    res1 <- results(dds1, contrast=c('mode_of_action', exp_group,'Control'))
    res1 <- lfcShrink(dds1, coef=2)
    
    # Order by adjusted p-val
    res2 <- as_tibble(res1, rownames=NA)
    res2 <- res2 %>% arrange(padj)
    
    # Append results to lists:
    deseq_list[[exp_group]] <- res1
    deseq_tib_list[[exp_group]] <- res2
    
    # Strip out weird characters from exp groups:
    exp_group_str <- str_replace_all(exp_group, "[[:punct:]]", "")
    
    # write out DE results
    write.csv(res2, paste(exp_group_str, 'deseq_results1.csv', sep="_"))
    
    # write out matrix of normalized counts
    write.csv(counts(dds1,normalized=TRUE),
              paste(exp_group_str, 'deseq_norm_counts1.csv', sep="_"))
    
  } else {
    print("There is an issue with sample order prior to deseq analysis")
  }
}

# Initialize list to store reformatted tibbles:
deseq_tib_list2 <- list()

for (exp_group in exp_modes_of_action){
  print(exp_group)
  # Convert rownames to own column:
  deseq_tib_list2[[exp_group]] <- rownames_to_column(deseq_tib_list[[exp_group]], var="Geneid")
  
  # Change column names:
  colnames(deseq_tib_list2[[exp_group]]) <- paste(exp_group, colnames(deseq_tib_list2[[exp_group]]), sep="_")
} 

# Create CSV with top 10 results:
top_ten_tib <- bind_cols(deseq_tib_list2$AhR[c(1:10), c(1, 6)],
                         deseq_tib_list2$`CAR/PXR`[c(1:10), c(1, 6)], 
                         deseq_tib_list2$DNA_Damage[c(1:10), c(1, 6)])
write.csv(top_ten_tib, "top_ten_diff_exp.csv", row.names=FALSE)


# Initialize list:
num_sig <- list()
deseq_tib_list_sig_only <- list()

for (exp_group in exp_modes_of_action){
  deseq_tib_list_sig_only[[exp_group]] <- deseq_tib_list[[exp_group]][deseq_tib_list[[exp_group]]$padj < 0.05,]
  num_sig[[exp_group]] <- length(deseq_tib_list_sig_only[[exp_group]]$padj)
} 

# Number of significant genes for each:
num_sig

# Write to CSV:
write.csv(data.frame(num_sig),"num_sig_genes_padj.csv", row.names=FALSE)

# Create histograms for significant genes
for (exp_group in exp_modes_of_action){
  # Strip out weird characters from exp groups:
  exp_group_str <- str_replace_all(exp_group, "[[:punct:]]", "")
  temp_filename <- paste(exp_group_str, "siggene_hist.jpg", sep="_")
  jpeg(temp_filename)
  hist(deseq_tib_list_sig_only[[exp_group]]$log2FoldChange, 
       breaks=30, xlim=c(-6, 6),
       main=exp_group, xlab="log2 Fold Change")
  dev.off()
} 

# Create scatter plots of fold change vs nominal p-value for significant genes
for (exp_group in exp_modes_of_action){
  # Strip out weird characters from exp groups:
  exp_group_str <- str_replace_all(exp_group, "[[:punct:]]", "")
  temp_filename <- paste(exp_group_str, "siggene_scatter.jpg", sep="_")
  jpeg(temp_filename)
  plot(deseq_tib_list_sig_only[[exp_group]]$log2FoldChange, 
       -log(deseq_tib_list_sig_only[[exp_group]]$pvalue),
       main=exp_group, xlab="log2 Fold Change", ylab="-log(nominal p-value)")
  dev.off()
} 

