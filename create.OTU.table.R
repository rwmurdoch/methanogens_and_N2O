## This is a simple script for combining outputs from qiime2 into a single OTU table which includes taxonomy and sequences
## The inputs are generated from qiime tools exports of sequences, feature table (via biom convert), and taxonomy

library(dplyr)
library(Biostrings)

#read in the feature table (this is already sorted high to low by combined abundance)
feature.table <- read.csv("exports/otu_table.txt", skip = 1, comment.char = "", sep = "\t")
colnames(feature.table)[1] <- "OTU_ID"

#read in the taxonomy
taxonomy <- read.csv("exports/taxonomy.tsv",sep="\t")
colnames(taxonomy) [1] <- "OTU_ID"

#read in the sequence fasta
q2.fasta <- Biostrings::readDNAStringSet("exports/dna-sequences.fasta")
OTU_ID <- names(q2.fasta)
sequence <- paste(q2.fasta)
q2.seq.tab <- data.frame(OTU_ID,sequence)

#join it all together
OTU_table <- dplyr::left_join(feature.table, taxonomy, by = "OTU_ID")
OTU_table <- dplyr::left_join(OTU_table,q2.seq.tab,by="OTU_ID")

#write it
write.csv(OTU_table,"exports/combined_table.csv")
