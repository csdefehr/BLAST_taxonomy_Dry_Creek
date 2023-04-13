setwd("/home/chadwick/Documents/C3VIP/taxonomyBLAST")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
library(rBLAST)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("taxize")
library(taxize)
#install.packages("taxonomizr")
library(taxonomizr)
#prepareDatabase(getAccessions=FALSE)
#Read in community matrix
communitymatrix <- read.csv("Data_DCsp_ASV.csv", header = TRUE, row.names = 1, 
                            check.names=FALSE)

#Function for finding top(n) of community matrix
abundance.by.sample <- function(matrix, number_of_ASV) {
  matrix <-as.matrix(matrix)  
  ranked_asv_df <- data.frame(matrix(ncol=number_of_ASV,nrow = 0))
  
  for (i in row.names(matrix)){
    output = names(sort(matrix[i,], decreasing = TRUE))[1:number_of_ASV]
    ranked_asv_df <- rbind.data.frame(ranked_asv_df, output)
  }
  row.names(ranked_asv_df) <- row.names(matrix) 
  colnames(ranked_asv_df) <- c(1:number_of_ASV)  
  return (ranked_asv_df)
}

#top_ten_ASV <- abundance.by.sample( communitymatrix, 10)

#BLAST:
#retrieve reference database
#download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences.tar.gz",
#                "ITS_eukaryote_sequences.tar.gz", mode='wb')

#Create BLAST database from download
bl = blast(db = "./ITS_eukaryote_sequences/ITS_eukaryote_sequences", type = "blastn")

#Create list of sequences for querying
seq <- readDNAStringSet("dna-sequences.fasta")

# Function to BLAST a given column of ranked ASV's
blast.ranked.asv <- function( ranked_asv_df, asv_ranking, initial_index_of_results,
                              final_index_of_results){
  blast_results <- data.frame(matrix(ncol=12,nrow = 0))
  
  for (i in row.names(ranked_asv_df)){
    raw_predict <- (predict(bl, seq[ranked_asv_df[i,asv_ranking],]))
    raw_predict <- raw_predict[order(-raw_predict$Bits),]
    blast_results <- 
    rbind.data.frame(blast_results, 
                     raw_predict[initial_index_of_results:final_index_of_results,])
  }
  return(blast_results)
  }

#ASV_blast <- blast.ranked.asv(top_ten_ASV, 1, 1, 1)

#Taxonomy:
#Function to find NCBI taxonomy of BLAST results
taxonomy.of.blast <- function( blast_results){
  taxonomy_by_sample <- data.frame(matrix(ncol=3, nrow=0))
  
  for (i in blast_results[,2]){
    first_taxa = getTaxonomy(genbank2uid(i), 'nameNode.sqlite')[7]   
    taxonomy_by_sample <- rbind.data.frame(taxonomy_by_sample, first_taxa)
  }
  return(taxonomy_by_sample)
  } 

#taxonomy_by_sample <- taxonomy.of.blast(ASV_blast)

#Pipe them together
final.function <- function(matrix, number_of_ASV, asv_ranking, initial_index, final_index){
test <- abundance.by.sample(matrix, number_of_ASV) %>% 
  blast.ranked.asv(asv_ranking, initial_index, final_index) %>% 
  taxonomy.of.blast()
  return(test)
}


#final_product <- data.frame(matrix(ncol=0,nrow = 193))
#final_product <- map_dfc(1:5, ~final.function(communitymatrix, 5,1, .x, .x))

taxa.by.asv <- function(matrix, asv_rank, n_of_blast_results){
  final_product <- map_dfc(1:n_of_blast_results,
                           ~final.function(matrix, n_of_blast_results,asv_rank, .x, .x))
  final_product <- as.data.frame(final_product)
  row.names(final_product) <- row.names(matrix)
  colnames(final_product) <- 1:n_of_blast_results
  return(final_product)
  }

test <- taxa.by.asv(communitymatrix, 1, 1) 
test2 <- taxa.by.asv(communitymatrix, 2, 1)


test4 <- data.frame(matrix(ncol=0,nrow = 193))

test2 <- for(i in 1:5){
  test3 <- taxa.by.asv(communitymatrix, i, 1)
  test4 <- cbind.data.frame(test4, test3[,1] )
  return(test4)
}

test <- map_dfc(1:10, ~blast.ranked.asv(top_ten_ASV, .x,1,1) %>% taxonomy.of.blast())
top_ten_species <- test
top_ten_species <- rownames_to_column(top_ten_species, "SampleID")

