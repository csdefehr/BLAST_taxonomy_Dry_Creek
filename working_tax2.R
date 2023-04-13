tax_try <- c()
for (a in 1:length(seq@ranges@NAMES)){
  print(paste("Currently running ASV:", seq@ranges@NAMES[a], sep = " "))
  print(seq[a])
  tax1 <- predict(bl2, seq[a]) #search database for sequence
  print(tax1[1:10,])
  tax1 <- subset(tax1, tax1$Bits == max(tax1$Bits)) #select only the entries with the highest number of bits
  taxnames <- c("superkingdom","phylum","class","order","family", "genus","species"  )
  tax7 <- data.frame(getTaxonomy(genbank2uid(tax1$SubjectID), 'nameNode.sqlite'), check.names = FALSE)
  colnames(tax7) <- taxnames
  print(tax7)
  if (length(unique(tax7[,7])) == 1) {
      tax7 <- tax7[1,]
    } else if (length(unique(tax7[,6])) == 1) {
      tax7 <- c(tax7[1,1:6], NA)
    } else if (length(unique(tax7[,6])) == 1) {
      tax7 <- c(tax7[1,1:6], NA)
    } else if (length(unique(tax7[,5])) == 1) {
      tax7 <- c(tax7[1,1:5], NA, NA)
    } else if (length(unique(tax7[,4])) == 1) {
      tax7 <- c(tax7[1,1:4], NA, NA, NA)
    } else if (length(unique(tax7[,3])) == 1) {
      tax7 <- c(tax7[1,1:3], NA, NA, NA, NA)
    } else if (length(unique(tax7[,2])) == 1) {
      tax7 <- c(tax7[1,1:2], NA, NA, NA, NA, NA)
    } else if (length(unique(tax7[,1])) == 1) {
      tax7 <- c(tax7[1,1], NA, NA, NA, NA, NA, NA)
    }
    names(tax7) <- taxnames
    tax_try <- rbind(tax_try,tax7)
    colnames(tax_try) <- taxnames
}


