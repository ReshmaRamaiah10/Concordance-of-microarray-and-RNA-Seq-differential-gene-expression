# Programmer: Reshma Ramaiah
# This script combines featureCounts output files into one csv file
# Additionally, it creates boxplots of distribution of counts

setwd("/projectnb/bf528/users/tinman_2022/project_3/Programmer/featureCounts_output/")

#Install janitor: 
if (FALSE){
  install.packages("janitor", repos = "http://cran.us.r-project.org")
}

# Import libraries:
library(tidyverse)
library(janitor)

# Get list of applicable file names:
files <- list.files(path="/projectnb/bf528/users/tinman_2022/project_3/Programmer/featureCounts_output",pattern="*.txt$")
num_files <- length(files)

# For every file:
for (i in 1:num_files){
  my_file <- files[i]
  temp_ds <- read.table(my_file, header=TRUE)
  # Grab 1st (Geneid) and last (counts)
  temp_ds <- temp_ds[,c(1, ncol(temp_ds))]
  temp_tib <- as_tibble(temp_ds)
  # Start combined tibble if first file, otherwise merge with existing combined tibble:
  if (i==1){
    comb_tib <- temp_tib
  }else{
    comb_tib <- full_join(comb_tib, temp_tib, by="Geneid")
  }
}

# Check column names:
names(comb_tib)
# Shorten column names:
names(comb_tib) <- sub("X.projectnb2.bf528.users.tinman_2022.project_3.star_results.", "", names(comb_tib))
names(comb_tib) <- sub("Aligned.sortedByCoord.out.bam", "", names(comb_tib))

# Write csv:
write.csv(comb_tib, file="featureCounts_combined.csv", row.names=FALSE)

# First boxplot of distribution of gene counts
jpeg("boxplot1.jpg")
boxplot(comb_tib[2:10], horizontal=TRUE, ylim=c(0,1000), las=1,
        pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green"),
        main="Distribution of Gene Counts", xlab="Gene Counts")
legend("topright", legend=c("CAR/PXR", "AhR", "DNA_Damage"), 
       fill=c("green", "red", "blue"), cex=0.75)
dev.off()

# Add 0.0001 for a logarithmic plot
comb_tib_shifted <- comb_tib
comb_tib_shifted[,c(2:10)] <- comb_tib[,c(2:10)] + 0.0001

# Second boxplot of distribution of gene counts in log scale
jpeg("boxplot2.jpg")
boxplot(comb_tib_shifted[2:10], horizontal=TRUE, las=1, log="x",
        pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green"),
        main="Distribution of Gene Counts", xlab="Gene Counts (log)")
legend("topright", legend=c("CAR/PXR", "AhR", "DNA_Damage"), 
       fill=c("green", "red", "blue"), cex=0.75)
dev.off()