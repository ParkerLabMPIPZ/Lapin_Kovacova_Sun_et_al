### extract and align sequences aligned to a query in blast
library(Biostrings)
library(msa)

# files with blast output
files_with_output <- c("EP_AtEDS1_BLAST_matches.txt",
                       "EP_AtPAD4_BLAST_matches.txt",
                       "EP_AtSAG101_BLAST_matches.txt")

# read complete seq database
db <- readAAStringSet("proteome_Phytozome_selected_species.fa", format = "fasta")

# read blast output
blast_output <- data.frame()
for (i in 1:length(files_with_output)){
  temp <- read.delim(files_with_output[i], skip = 5, header = FALSE)
  temp <- temp[-nrow(temp),] # last row is blast comment
  blast_output <- rbind(blast_output,temp)
  rm(temp)
}


# give columns names for easier navigation
colnames(blast_output) <- c("query",
                            "subject",
                            "e_value",
                            "start_in_query",
                            "end_in_query",
                            "start_in_subject",
                            "end_in_subject")

# calculate match size
blast_output$match_size <- blast_output$end_in_subject - blast_output$start_in_subject

# plot domain size distribution
pdf("EP_domain_size.pdf")
hist(blast_output$match_size,
     main = NA,
     xlab = "length of Arabidopsis EP domain matches in plants",
     ylab = "# of sequences")
dev.off()

# extract seq data for a given domain
selected_seq <- AAStringSet()
for (i in 1:nrow(blast_output)){
temp <- subseq(db[which(db@ranges@NAMES %in% blast_output$subject[i])],
               start = blast_output$start_in_subject[i],
               end = blast_output$end_in_subject[i])
selected_seq <- c(selected_seq, temp)
rm(temp)
}

# remove duplicated seq
selected_seq <- selected_seq[!duplicated(selected_seq)]

# align seq
aligned_seq <- msa(selected_seq, method = "Muscle")

msaPrettyPrint(x=aligned_seq,output="tex",file = "alignedSeq_Muscle.tex",alFile = "alignedSeq.fasta",
               showLogo = "top",shadingColors="blues",
               showConsensus = "none",
               consensusThreshold = 50,
               askForOverwrite=FALSE)

