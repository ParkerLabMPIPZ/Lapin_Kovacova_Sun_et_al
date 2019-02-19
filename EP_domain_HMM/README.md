# Pipeline to build EP domain HMM
### run blastp against nonredundant protein RefSeq on BLAST website evalue=0.01

>AtEDS1_EP_domain Q9XF23 385-623
LKKLAWIEDEYKPKCQAHKNGYYDSFKVSNEENDFK
ANVKRAELAGVFDEVLGLLKKCQLPDEFEGDIDWIKLATRYRRLVEPLDIANYHRHLKNE
DTGPYMKRGRPTRYIYAQRGYEHHILKPNGMIAEDVFWNKVNGLNLGLQLEEIQETLKNS
GSECGSCFWAEVEELKGKPYEEVEVRVKTLEGMLREWITAGEVDEKEIFLEGSTFRKWWI
TLPKNHKSHSPLRDYMMDEITDT

>AtPAD4_EP_domain Q9S745 300-541
SAELANELASVLPARLEIQWYKDRCDASEEQLGYYDFFKRYSLKRDFKVNMSRIRLAKFWD
TVIKMVETNELPFDFHLGKKWIYASQFYQLLAEPLDIANFYKNRDIKTGGHYLEGNRPKR
YEVIDKWQKGVKVPEECVRSRYASTTQDTCFWAKLEQAKEWLDEARKESSDPQRRSLLRE
KIVPFESYANTLVTKKEVSLDVKAKNSSYSVWEANLKEFKCKMGYENEIEMVVDESDAME
T

>AthSAG101_EP_domain Q4F883 291-537
DMMFKKLNDM
KISMAYIEWYKKKCKEVKIGYYDRFKTQLAFPSKEFDINIKNHHKSELNRFWKSVVEEVE
RRPQSDASILKRRFLFSGNNYRRMIEPLDIAEYYLEGRKEYRTTGRSHHYVMLEKWFGME
SILIEKERCKKRDLSDLLTFDSCFWAEVEDSLIVINQLNTTVGMRDDVREVLTRKLVEFE
GYVWEIITKREVSPEIFLEESSFMKWWKEYKKIKGFNSSYLTEFMNTRKYESYGKSQ

### against refseq89:

EDS1: one monocot and no gymnosperms, the database is incomplete

PAD4: only eudicots, the database is incomplete

SAG101: only eudicots

### since plant diversity is not well represented in refseq89, search for EP-containing proteins will be performed with custom Phytozome db
#### proteomes are obtained via BioMart for 32 species on 07.08.2018
1.	Physcomitrella patens
2.	Sphagnum fallax
3.	Marchantia polymorpha v3.1
4.	Selaginella moellendorffii
5.	Ananas comosus
6.	Musa acuminata
7.	Amborella trichopoda v1.0
8.	Oryza sativa
9.	Brachypodium distachyon
10.	Zea mays
11.	Sorghum bicolor
12.	Aquilegia coerulea
13.	Amaranthus hypochondriacus
14.	Kalanchoe laxiflora
15.	Daucus carota
16.	Solanum tuberosum
17.	Solanum lycopersicum
18.	Eucalyptus grandis 
19.	Vitis vinifera Genoscope
20.	Linum usitatissimum
21.	Populus trichocarpa
22.	Citrus clementina
23.	Theobroma cacao
24.	Carica papaya
25.	Arabidopsis thaliana TAIR10
26.	Capsella rubella
27.	Cucumis sativus
28.	Fragaria vesca
29.	Glycine max
30.	Volvox carteri
31.	Micromonas sp. RCC299
32.	Chlamydomonas reinhardtii


### make a database for BLAST
makeblastdb -in proteome_Phytozome_selected_species.fa -dbtype prot
#### with parse_seqid option warnings are produced, without parsing speciesID is kept -> sufficient for further analysis

### run BLASTP
blastp -db proteome_Phytozome_selected_species.fa -query AtEDS1_EP.fa -evalue 0.01 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out EP_AtEDS1_BLAST_matches.txt

blastp -db proteome_Phytozome_selected_species.fa -query AtPAD4_EP.fa -evalue 0.01 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out EP_AtPAD4_BLAST_matches.txt

blastp -db proteome_Phytozome_selected_species.fa -query AtSAG101_EP.fa -evalue 0.01 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out EP_AtSAG101_BLAST_matches.txt

### extract and align (MUSCLE) sequences matched with EP domains - run R script "extract_BLAST_aligned_seq_EP.R"


### build HMM
hmmbuild EP_domain.hmm alignedSeq.fasta
