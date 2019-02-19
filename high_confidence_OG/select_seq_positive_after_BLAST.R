library(Biostrings)



##### for EDS1 #####
# read blast results
blast_results <- read.table("AtEDS1_BLAST_matches.txt", sep = "\t")

# select only seq matching EDS1 as a significant hit
blast_results <- blast_results[blast_results$V2 %in% c("AT3G48090.1", "AT3G48080.1"),]
sign_matches <- as.vector(unique(blast_results$V1))

# extract sequences for significant matches
seq_all <- readAAStringSet("EDS1_all_orthologs.fa")
length(seq_all) # original dataset 85 seq, all are sign EDS1 matches length(sign_matches) -> 85 too
seq_sign <- seq_all[seq_all@ranges@NAMES %in% sign_matches]
length (seq_sign) # ok
writeXStringSet(seq_sign, "EDS1_sign_BLAST.fa")
rm(blast_results, sign_matches, seq_all, seq_sign)


##### for PAD4 #####
# read blast results
blast_results <- read.table("AtPAD4_BLAST_matches.txt", sep = "\t")

# select only seq matching PAD4 as a significant hit
blast_results <- blast_results[blast_results$V2 %in% c("AT3G52430.1"),]
sign_matches <- as.vector(unique(blast_results$V1))

# extract sequences for significant matches
seq_all <- readAAStringSet("PAD4_all_orthologs.fa")
length(seq_all) # original dataset 105 seq, all are sign PAD4 matches length(sign_matches) -> 105 too
seq_sign <- seq_all[seq_all@ranges@NAMES %in% sign_matches]
length (seq_sign) # ok
writeXStringSet(seq_sign, "PAD4_sign_BLAST.fa")
rm(blast_results, sign_matches, seq_all, seq_sign)


##### for SAG101 #####
# read blast results
blast_results <- read.table("AtSAG101_BLAST_matches.txt", sep = "\t")

# select only seq matching SAG101 as a significant hit
blast_results <- blast_results[blast_results$V2 %in% c("AT5G14930.2"),]
sign_matches <- as.vector(unique(blast_results$V1))

# extract sequences for significant matches
seq_all <- readAAStringSet("SAG101_all_orthologs.fa")
length(seq_all) # original dataset 155 seq, sign SAG101 matches length(sign_matches) -> 102, OG contaminated
seq_sign <- seq_all[seq_all@ranges@NAMES %in% sign_matches]
length (seq_sign) # ok
writeXStringSet(seq_sign, "SAG101_sign_BLAST.fa")
