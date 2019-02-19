# Pipeline to filter OrthoMCL orthogroups for EDS1, PAD4 and SAG101
the result is a set of high confidence OGs added to the database with presense/absence profiles

### make a database for BLAST
makeblastdb -in TAIR10_pep_20110103_representative_gene_model_updated.fa -dbtype prot

### run BLASTP
blastp -db TAIR10_pep_20110103_representative_gene_model_updated.fa -query EDS1_all_orthologs.fa -evalue 0.00001 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out AtEDS1_BLAST_matches.txt

blastp -db TAIR10_pep_20110103_representative_gene_model_updated.fa -query PAD4_all_orthologs.fa -evalue 0.00001 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out AtPAD4_BLAST_matches.txt

blastp -db TAIR10_pep_20110103_representative_gene_model_updated.fa -query SAG101_all_orthologs.fa -evalue 0.00001 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out AtSAG101_BLAST_matches.txt

blastp -db TAIR10_pep_20110103_representative_gene_model_updated.fa -query PlantImmunity100829_all_orthologs.fa -evalue 0.00001 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out AtADR1_BLAST_matches.txt

blastp -db TAIR10_pep_20110103_representative_gene_model_updated.fa -query PlantImmunity100583_all_orthologs.fa -evalue 0.00001 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out AtNRG1_BLAST_matches.txt

### select sequences matching seq of interest in TAIR10
run "select_seq_positive_after_BLAST.R"

### remove illegal sign in PAD4_sign_BLAST.fa -> "! " in lines 411 and 445

### test BLAST-sign hits for presence of EP domain
hmmsearch -o ./HMMoutput/hmmsearch_EP_EDS1.txt -A ./HMMoutput/alignment_EP_EDS1.txt --domtblout ./HMMoutput/per_domain_hits_EP_EDS1.txt --incE 0.00001 ./HMMs/EP_domain.hmm EDS1_sign_BLAST.fa

hmmsearch -o ./HMMoutput/hmmsearch_EP_PAD4.txt -A ./HMMoutput/alignment_EP_PAD4.txt --domtblout ./HMMoutput/per_domain_hits_EP_PAD4.txt --incE 0.00001 ./HMMs/EP_domain.hmm PAD4_sign_BLAST.fa

hmmsearch -o ./HMMoutput/hmmsearch_EP_SAG101.txt -A ./HMMoutput/alignment_EP_SAG101.txt --domtblout ./HMMoutput/per_domain_hits_EP_SAG101.txt --incE 0.00001 ./HMMs/EP_domain.hmm SAG101_sign_BLAST.fa

### select sequences with EP domain >50 aa
run "from_hmmsearch_to_fasta_full_with_EP_EDS1.R"

run "from_hmmsearch_to_fasta_full_with_EP_PAD4.R"

run "from_hmmsearch_to_fasta_full_with_EP_SAG101.R"

### select sequences with length 400-1200 aa
run "select_seq_based_on_length.R"
