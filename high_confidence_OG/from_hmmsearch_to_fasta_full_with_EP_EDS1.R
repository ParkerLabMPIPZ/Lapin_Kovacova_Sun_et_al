library(Biostrings)
INPUT_FASTA_FILE<-"EDS1_sign_BLAST.fa"
RESULT_FILE_FROM_HMM<-"./HMMoutput/per_domain_hits_EP_EDS1.txt"
OUTPUT_FILE<-"EDS1_sequences_with_EP.fa"
OUTPUT_ALL_MATCH_INFORMATION <- "./HMMoutput/EDS1_EP_characteristics_all_matches.txt"
OUTPUT_HISTOGRAM <- "length_distribution_EDS1_EP.pdf"
MINIMAL_DOMAIN_LENGTH <- 50 # check histogram to determine the appropriate size
DOMAIN_NAME <- "EP"

################################ the rest of the code does not need adjustment unless STO format is not 1.0 #########
###specify input fasta file
input_fasta<-readAAStringSet(INPUT_FASTA_FILE)
#get nr of seq in the input file
n.seq<-length(input_fasta) #it is however not absolutely critical
print(n.seq)


### select significant hits and their positions
HMM_result <- read.table(RESULT_FILE_FROM_HMM, header = FALSE, skip = 3)
# keep only sign hits sequence E<=0.00001
domain_match <- HMM_result[HMM_result$V7<=0.00001,c("V1","V18","V19")]
colnames(domain_match) <- c("Seq_Name", "Domain_start","Domain_end")

domain_match<-domain_match[complete.cases(domain_match),]

# determine matched domain size
domain_match$domain_size <- domain_match$Domain_end - domain_match$Domain_start

# since each protein can have several domains, it is useful to keep them all but under different names
domain_match$Seq_name_individ_matches <- paste(domain_match$Seq_Name,
                                               domain_match$Domain_start,
                                               domain_match$Domain_end, sep="_")

write.table(domain_match, OUTPUT_ALL_MATCH_INFORMATION,
            quote= FALSE, row.names = FALSE, sep = "\t")

#plot histogram of length distribution
pdf(OUTPUT_HISTOGRAM)
hist(domain_match$domain_size,
     main = paste0(DOMAIN_NAME, "-domain length"),
     xlab = "domain length",
     ylab = "counts",
     col = "#bdbdbd")
dev.off()

########write selected seq into a file, check the df domain_match to see appr. size of the domain
# keep only entries with length above threshold
domain_match <- domain_match[domain_match$domain_size >= MINIMAL_DOMAIN_LENGTH,]

# select seq
selected_seq <- input_fasta[input_fasta@ranges@NAMES %in% unique(domain_match$Seq_Name)]
length(selected_seq) # all input seq contain the domain? compare with print(n.seq). If the same, all seq contain domain

#write extracted seq into a fasta file
writeXStringSet(selected_seq,
                file=OUTPUT_FILE)
