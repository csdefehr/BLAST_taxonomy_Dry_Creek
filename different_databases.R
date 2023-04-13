setwd("/home/chadwick/Documents/C3VIP/taxonomyBLAST")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
library(rBLAST)
#library(RSQLite)
# Make new db from cultured fasta file
#makeblastdb("TrimmedSagebrushPhyllospherePhylogeny.fasta", dbtype = "nucl")
bl1 <- blast(db = "./TrimmedSagebrushPhyllospherePhylogeny.fasta", type = "blastn")

# Make db from searching ITS_RefSeq_Fungi with term Leaf results
#makeblastdb("FungiOrganismAND.fasta", dbtype = "nucl")
bl2 <- blast(db = "PRJNA177353_AND_Leaf.fasta", type = "blastn")

# Make ITS_RefSeq_Fungi db
bl3 = blast(db = "./ITS_eukaryote_sequences/ITS_eukaryote_sequences", type = "blastn")

#Create list of sequences for querying
seq <- readDNAStringSet("dna-sequences.fasta")

#Blast sequences from ASV-sequence file against new database
test1 <- predict(bl1, seq[100])
test2 <- predict(bl2, seq[100])
test3 <- predict(bl3, seq[100])
ncbi_fulltax <- read.csv("DC_taxonomy_ncbi_full.csv")





