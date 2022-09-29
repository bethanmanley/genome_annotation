########################################################  
# Remove protein coding genes from repeats modeled by EDTA
########################################################  
library(data.table)
library(tidyverse)

#InterPro scan the EDTA repeat library
./interproscan.sh -i Rirr.curated_primary.no_mt.unscrubbed.fa.mod.EDTA.TElib.fa -t n -goterms -iprlookup -pa -d /TE/
  
#Convert InterProScan output as csv, then:
repeat_ipr <- read.csv("Rirr.curated_primary.no_mt.unscrubbed.fa.mod.EDTA.TElib.csv", header=F, stringsAsFactors=FALSE,  fill=TRUE)

#Split TE column to separate orfs
repeat_ipr <- repeat_ipr %>% 
  extract(V1, into = c("TE", "ORF"), "(.*)_([^_]+)$")

#Column V6 contains protein domain names
IPRs <- data.frame(unique(repeat_ipr$V6))
# To filter out :
Cys/Met metabolism
Ion transport protein
Transferase(Phosphotransferase) 
Actin
BTB
Potassium Channel
Myosin
Glycogen
G1 TO S PHASE TRANSITION 1
PROTEIN CBG23806
RING/U-box
Ricin B-like lectins
ELECTRON TRANSPORT OXIDOREDUCTASE
PALMITOYLTRANSFERASE
Argonaute linker 2 domain
Frizzled cysteine-rich domain
HOMOCYSTEINE/CYSTEINE SYNTHASE
F-box-like
MAC/Perforin domain
TLDc domain
HAMP domain
Acetyl-CoA synthetase-like
CheY-like
THROMBOXANE-A SYNTHASE
p450
TREHALOSE-PHOSPHATASE
FAD binding domain
FERM domain profile
K+ potassium transporter
MyTH4
IQ motif profile
PROTEIN TAR1
Type I antifreeze protein signature
OTU domain profile
TPR
Ubiquitin-like
Multiheme cytochromes
REJ domain profile
Crinkler effector protein N-terminal domain
Collagenase
F-box
Tetratricopeptide
INVOLUCRIN
BPG-independent PGAM N-terminus (iPGM_N)
HMGL-like
Beta-lactamase
TIP49 P-loop domain
nexin
TPR repeat profile
Ribonuclease III
2,3-Bisphosphoglycerate-independent phosphoglycerate mutase, substrate-binding domain
NERD domain profile
KELCH
sel1
N-6 Adenine-specific DNA methylases signature.
S-adenosyl-L-methionine-dependent methyltransferases
RNA METHYLASE-RELATED
Caspase
F-actin
von Willebrand factor, type A domain
AUTOPHAGY-RELATED 2, ISOFORM A
Protein kinases ATP-binding region signature.
tRNA
fungal_TF_MHR
ATP-DEPENDENT RNA HELICASE
MYND
P450
SGNH
NAD(P)H
Aldolase
homocitrate
Vps51/Vps67
SERINE/THREONINE-PROTEIN
amino acid

toremove <- dplyr::filter(repeat_ipr, grepl("Cys/Met metabolism|Ion transport protein|Transferase(Phosphotransferase) |Actin|BTB|Potassium Channel|Myosin|Glycogen|G1 TO S PHASE TRANSITION 1|PROTEIN CBG23806|RING/U-box|Ricin B-like lectins|ELECTRON TRANSPORT OXIDOREDUCTASE|PALMITOYLTRANSFERASE|Argonaute linker 2 domain|Frizzled cysteine-rich domain|HOMOCYSTEINE/CYSTEINE SYNTHASE|F-box-like|MAC/Perforin domain|TLDc domain|HAMP domain|Acetyl-CoA synthetase-like|CheY-like|THROMBOXANE-A SYNTHASE|p450|TREHALOSE-PHOSPHATASE|FAD binding domain|FERM domain profile|K+ potassium transporter|MyTH4|IQ motif profile|PROTEIN TAR1|Type I antifreeze protein signature|OTU domain profile|TPR|Ubiquitin-like|Multiheme cytochromes|REJ domain profile|Crinkler|collagenase|f-box|Tetratricopeptide|Tetratricopeptide-like|involucrin|iPGM_N|HMGL-like|Beta-lactamase|TIP49|nexin|TPR repeat profile|Ribonuclease III|phosphoglycerate|NERD|KELCH|sel1|DNA methylases signature|S-adenosyl-L-methionine-dependent methyltransferases|RNA METHYLASE-RELATED|Caspase|F-actin|von Willebrand factor, type A domain|AUTOPHAGY-RELATED 2, ISOFORM A|Protein kinases ATP-binding region signature.|tRNA|fungal_TF_MHR|ATP-DEPENDENT RNA HELICASE|MYND|P450|SGNH|NAD(P)H|Aldolase|homocitrate|Vps51/Vps67|SERINE/THREONINE-PROTEIN|amino acid", V6, ignore.case=TRUE))
toremove_list <- data.frame(unique(toremove$TE))

consensus_OG <- read.table("Rirr.curated_primary.no_mt.unscrubbed.fa.mod.EDTA.TElib.fa", header = F, sep="\t", quote="", comment.char = "")
consensus_OG <- as.data.frame(matrix(consensus_OG$V1, ncol = 2, byrow = TRUE))
consensus_OG$V1 <- gsub(">", "", consensus_OG$V1)
length(unique(consensus_OG$V1)) #1786 consensus sequences
names(toremove_list)[1] <- "V1"
antijoin <- anti_join(consensus_OG, toremove_list, by='V1')

#Reformat into fasta, and rename TE IDs with RepeatMasker format
antijoin$V1 <- paste(">", antijoin$V1, sep="")
antijoin$V1 <- gsub("DNA/DTA", "DNA/hAT", antijoin$V1)
antijoin$V1 <- gsub("DNA/DTC", "DNA/CMC-EnSpm", antijoin$V1)
antijoin$V1 <- gsub("DNA/DTH", "DNA/PIF-Harbinger", antijoin$V1)
antijoin$V1 <- gsub("DNA/DTM", "DNA/MULE-MuDR", antijoin$V1)
antijoin$V1 <- gsub("DNA/DTT", "DNA/TcMar-Mariner", antijoin$V1)
antijoin$V1 <- gsub("TIR/Tc1_Mariner", "DNA/PIF-Harbinger", antijoin$V1)
antijoin$V1 <- gsub("polinton", "DNA/Maverick", antijoin$V1)
antijoin$V1 <- gsub("LTR/unknown", "LTR", antijoin$V1)
antijoin$V1 <- gsub("DIRS", "LTR/DIRS", antijoin$V1)
antijoin$V1 <- gsub("LINE/unknown", "LINE", antijoin$V1)
antijoin$V1 <- gsub("DNA/Helitron", "RC/Helitron", antijoin$V1)

unique(antijoin$V1)

antijoin2 <- antijoin %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "TE",
    values_to = "seq")
antijoin2$TE <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_curatedrepeatlibrary_firstpass.fasta", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)


########################################################  
#Identify LINEs from Repeat Modeler output
########################################################  

#Remove line breaks from consensus seq fasta
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' consensi.fa.classified > consensi.fa.classified.fa

Rir_RM <- read.table("consensi.fa.classified.fa", header = F, sep="\t", quote="", comment.char = "")
Rir_RM <- as.data.frame(matrix(Rir_RM$V1, ncol = 2, byrow = TRUE))
length(unique(Rir_RM$V1)) #1700 families
#Simplify fasta headers 
Rir_RM$V1 <- sub(" .*", "", Rir_RM$V1)
#Isolate LINEs
Rir_LINEs <- dplyr::filter(Rir_RM, grepl("LINE", V1, ignore.case=TRUE))
#Reformat into fasta
Rir_LINEs <- Rir_LINEs %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "TE",
    values_to = "seq")
Rir_LINEs$TE <- NULL
write.table(Rir_LINEs, paste(out_dir,"consensi.fa.classified.LINEs.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

#Blast 'consensi.fa.classified.LINEs.fa' against 'Manley.EDTA.TElib.mod.fa' To find unclassified repeat families that are LINEs
~/ncbi-blast-2.13.0+/bin/makeblastdb -in ~/genomics/rirman22/TE/Rirr.curated_primary.no_mt.unscrubbed.fa.mod.EDTA.TElib_curated1.fa -title RirTEdb -dbtype nucl -out  ~/genomics/rirman22/TE/RirTE_ntdb
~/ncbi-blast-2.13.0+/bin/blastn -task blastn -query ~/genomics/rirman22/TE/Rmodeler/RM_1513.WedJun292001092022/consensi.fa.classified.LINEs.fa -db ~/genomics/rirman22/TE/RirTE_ntdb -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 300 -max_hsps 1 -evalue 1e-10 -out ~/genomics/rirman22/TE/Rmodeler/Rir_LINEs_blastn.out
#Classified by blastn of RepeatModeler LINEs to EDTA consensus sequences. 
TE_00000012#Unknown > rnd-1_family-139#LINE/I # But has no predicted protein domain. Leave uknown
TE_00000509#Unknown > rnd-3_family-305#LINE   # But has no predicted protein domain. Leave uknown
TE_00000149#Unknown > rnd-1_family-661#LINE/Proto1  # RNAseH domain, this is a LINE
TE_00000338#Unknown > rnd-1_family-661#LINE/Proto1  # But has no predicted protein domain. Leave uknown
TE_00000956#Unknown > rnd-5_family-1722#LINE/L2     # DNAseI, endonuclease domains, this is a LINE
TE_00000385#Unknown > rnd-4_family-612#LINE         # But has no predicted protein domain. Leave uknown

# Based on this and on checking protein-coding domains, I rename the following 'unknowns'
TE_00000149#Unknown > TE_00000149#LINE
TE_00000956#Unknown > TE_00000956#LINE

#TEs with a DUF 659 domain are DNA/hAT: I re-classified the following:
antijoin$V1 <- gsub("TE_00000149#Unknown", "TE_00000149#LINE", antijoin$V1)
antijoin$V1 <- gsub("TE_00000956#Unknown", "TE_00000956#LINE", antijoin$V1)
antijoin$V1 <- gsub("TE_00000126#Unknown", "TE_00000126#DNA/hAT", antijoin$V1)
antijoin$V1 <- gsub("TE_00001358#DNA/MULE-MuDR", "TE_00001358#DNA/hAT", antijoin$V1)
antijoin$V1 <- gsub("TE_00001486#DNA/CMC-EnSpm", "TE_00001486#DNA/MULE-MuDR", antijoin$V1)
antijoin$V1 <- gsub("TE_00001595#DNA/hAT", "TE_00001595#DNA/MULE-MuDR", antijoin$V1)
antijoin$V1 <- gsub("TE_00001365#DNA/hAT", "TE_00001365#DNA/MULE-MuDR", antijoin$V1)
antijoin$V1 <- gsub("TE_00001707#DNA/Helitron", "TE_00001707#Unknown", antijoin$V1)

antijoin2 <- antijoin %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "TE",
    values_to = "seq")
antijoin2$TE <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_curatedrepeatlibrary.fasta", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

# Then I run RepeatMsker using this curated TE library.

