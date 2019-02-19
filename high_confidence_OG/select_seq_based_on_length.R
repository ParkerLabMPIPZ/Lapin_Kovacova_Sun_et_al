library(Biostrings)

#### remove too long and too short EDS1, PAD4, SAG101 sequences ####

# read sequences
EDS1 <- readAAStringSet("EDS1_sequences_with_EP.fa")
PAD4 <- readAAStringSet("PAD4_sequences_with_EP.fa")
SAG101 <- readAAStringSet("SAG101_sequences_with_EP.fa")

# check sequence length distribution
pdf("selected_seq_length_distribution.pdf")
par(mfrow=c(3,1))
hist(EDS1@ranges@width,
     xlim = c(0,2000),
     main = "putative EDS1 sequence length",
     col = "red",
     xlab = "AA seq length",
     ylab = "nr. seq")

hist(PAD4@ranges@width,
     xlim = c(0,2000),
     main = "putative PAD4 sequence length",
     col = "blue",
     xlab = "AA seq length",
     ylab = "nr. seq")

hist(SAG101@ranges@width,
     xlim = c(0,2000),
     main = "putative SAG101 sequence length",
     col = "green",
     xlab = "AA seq length",
     ylab = "nr. seq")
dev.off()

# set limits for seq
min_length <- 400
max_length <- 1200

# select seq with length between limits
EDS1 <- EDS1[EDS1@ranges@width >= min_length & EDS1@ranges@width <= max_length]
PAD4 <- PAD4[PAD4@ranges@width >= min_length & PAD4@ranges@width <= max_length]
SAG101 <- SAG101[SAG101@ranges@width >= min_length & SAG101@ranges@width <= max_length]

# check sequence length distribution after removing too long and too short seq
pdf("selected_seq_length_distribution_400_1200aa.pdf")
par(mfrow=c(3,1))
hist(EDS1@ranges@width,
     xlim = c(0,2000),
     main = "putative EDS1 sequence length",
     col = "red",
     xlab = "AA seq length",
     ylab = "nr. seq")

hist(PAD4@ranges@width,
     xlim = c(0,2000),
     main = "putative PAD4 sequence length",
     col = "blue",
     xlab = "AA seq length",
     ylab = "nr. seq")

hist(SAG101@ranges@width,
     xlim = c(0,2000),
     main = "putative SAG101 sequence length",
     col = "green",
     xlab = "AA seq length",
     ylab = "nr. seq")
dev.off()

# write seq in files
writeXStringSet(EDS1, "EDS1_high_confidence_400_1200.fa")
writeXStringSet(PAD4, "PAD4_high_confidence_400_1200.fa")
writeXStringSet(SAG101, "SAG101_high_confidence_400_1200.fa")

# write names of selected seq into a separate file
write.table(EDS1@ranges@NAMES, "names_EDS1_high_confidence_400_1200.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PAD4@ranges@NAMES, "names_PAD4_high_confidence_400_1200.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(SAG101@ranges@NAMES, "names_SAG101_high_confidence_400_1200.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
