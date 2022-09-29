
setwd("/genes/functional_annot/Illumina")
setwd("/genes/functional_annot/Illumina_ONT")

#gffread -E Rhizophagus_irregularis_DAOM197198.gff3 -T -o Rhizophagus_irregularis_DAOM197198.gtf

#Functional annotation
GO <- read.table("Rhizophagus_irregularis_DAOM197198.annotations.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
GO <- GO %>% mutate_all(na_if,"")

################################################################################          
# Curate gff3
################################################################################          
gene.gff3 <- read.table("Rhizophagus_irregularis_DAOM197198.gff3", header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")

# Transposon IPR domains to remove. Check folowing files for descriptions
# TEIPR <- unique(read.table("TE_IPRdomains.tsv", sep="\t"))
# Genes to keep (manually checked and selected)
# tokeep <- unique(read.table("tokeep.tsv", sep="\t"))

# Filter transposon-related InterPro domains
toremove <- dplyr::filter(GO, grepl("IPR041078|IPR004242|IPR027806|IPR041577|IPR043128|IPR018289|IPR045249|IPR002492|IPR045259|IPR025525|IPR040521|IPR041457|IPR000477|IPR038717|IPR000477|IPR041373|IPR041577|IPR004875|IPR031052|IPR008906|IPR018473|IPR036691|IPR036236|IPR013087|IPR012337|IPR006758|IPR011010|IPR045455", InterPro, ignore.case=TRUE))

# ***ILLUMINA*** Keep genes like telomerase (for ex. g12233 has a reverse transcriptase domain but isn't a TE)
toremove_final <- dplyr::filter(toremove, !grepl("g12233|g7167|g1782|g3543|g12233|g7951|g6777|g16786|g24896|g29429|g12721|g16420|g10014|g12233|g9634|g18092", GeneID, ignore.case=TRUE))
toremove_final <- data.frame(toremove_final[, c("GeneID")])
names(toremove_final)[1] <- "yup"

# ***ILLUMINA+ONT*** Keep genes like telomerase (for ex. g12233 has a reverse transcriptase domain but isn't a TE)
toremove_final <- dplyr::filter(toremove, !grepl("g12233|g7167|g1782|g3543|g12233|g7951|g6777|g16786|g24896|g29429|g12721|g16420|g10014|g12233|g9634|g18092", GeneID, ignore.case=TRUE))
toremove_final <- data.frame(yup = c(toremove_final[,"GeneID"], toremove_final[,"TranscriptID"]))

#split gff3 and isolate geneid
gene.gff3_split1 <- data.frame(str_split_fixed(gene.gff3$V9, ";", 3))
gene.gff3_split1$X1 <- gsub("\\ID=","", gene.gff3_split1$X1)
gene.gff3_split1$X1 <- gsub("\\..*","",gene.gff3_split1$X1)
gene.gff3_plusids <- cbind(gene.gff3, gene.gff3_split1$X1)
names(gene.gff3_plusids)[10] <- "yup"


#Remove 'toremove_final' from our gff3 based on partial string match
antijoin <- anti_join(gene.gff3_plusids, toremove_final, by='yup')
antijoin$yup <- NULL

#Check original number of genes
OG_genes <- gene.gff3_plusids[gene.gff3_plusids$V3 == "gene", ] 
#Illumina = 31666
#Illumina+ONT = 31654
Filtered_genes <- antijoin[antijoin$V3 == "gene", ]
#Illumina = 30230
#Illumina+ONT = 30209

write.table(antijoin, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_Illumina_curated.gff3", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

write.table(antijoin, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_Illumina+ONT_curated.gff3", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)


################################################################################          
# Curate gtf
################################################################################ 
gene.gtf <- read.table("Rhizophagus_irregularis_DAOM197198.gtf", header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")

#split gtf and isolate geneid ***Illumina***
gene.gtf_split1 <- data.frame(str_split_fixed(gene.gtf$V9, ";", 2))
gene.gtf_split1$X1 <- gsub("transcript_id \"","", gene.gtf_split1$X1)
gene.gtf_split1$X1 <- gsub("\\..*","",gene.gtf_split1$X1)
gene.gtf_plusids <- cbind(gene.gtf, gene.gtf_split1$X1)
names(gene.gtf_plusids)[10] <- "yup"

#split gtf and isolate geneid ***Illumina+ONT***
gene.gtf_split1 <- data.frame(str_split_fixed(gene.gtf$V9, ";", 2))
gene.gtf_split1$X2 <- gsub("gene_id \"","", gene.gtf_split1$X2)
gene.gtf_split1$X2 <- gsub("[^A-Za-z0-9]", "", gene.gtf_split1$X2)
gene.gtf_plusids <- cbind(gene.gtf, gene.gtf_split1$X2)
names(gene.gtf_plusids)[10] <- "yup"

#Remove 'toremove_final' from our gtf based on partial string match
antijoin2 <- anti_join(gene.gtf_plusids, toremove_final, by='yup')
antijoin2$yup <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_Illumina_curated.gtf", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_Illumina+ONT_curated.gtf", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)


################################################################################          
# Curate cds-transcripts
################################################################################ 
perl -pe '/^>/ ? print "\n" : chomp' Rhizophagus_irregularis_DAOM197198.cds-transcripts.fa > Rhizophagus_irregularis_DAOM197198.cds-transcripts_l.fa
perl -pe '/^>/ ? print "\n" : chomp' Rhizophagus_irregularis_DAOM197198.mrna-transcripts.fa > Rhizophagus_irregularis_DAOM197198.mrna-transcripts_2.fa
perl -pe '/^>/ ? print "\n" : chomp' Rhizophagus_irregularis_DAOM197198.proteins.fa > Rhizophagus_irregularis_DAOM197198.proteins_2.fa

setwd("/genes/functional_annot/Illumina")
setwd("/genes/functional_annot/Illumina_ONT")

gene.cds <- read.csv("Rhizophagus_irregularis_DAOM197198.cds-transcripts.fa", header=F)
gene.cds2 <- as.data.frame(matrix(gene.cds$V1, ncol = 2, byrow = TRUE))
#Split fasta header
gene.cds2_split1 <- data.frame(str_split_fixed(gene.cds2$V1, " ", 2))
gene.cds2 <- cbind(gene.cds2, gene.cds2_split1$X2)
names(gene.cds2)[3] <- "yup"

antijoin <- anti_join(gene.cds2, toremove_final, by='yup')
antijoin$yup <- NULL

#Reformat to fasta 
antijoin2 <- antijoin %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "ID",
    values_to = "seq")
antijoin2$ID <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.cds-transcripts_Illumina_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.cds-transcripts_Illumina+ONT_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

################################################################################          
# Curate mRNA-transcripts
################################################################################ 

setwd("/genes/functional_annot/Illumina")
setwd("/genes/functional_annot/Illumina_ONT")

gene.cds <- read.csv("Rhizophagus_irregularis_DAOM197198.mrna-transcripts.fa", header=F)
gene.cds2 <- as.data.frame(matrix(gene.cds$V1, ncol = 2, byrow = TRUE))
#Split fasta header
gene.cds2_split1 <- data.frame(str_split_fixed(gene.cds2$V1, " ", 2))
gene.cds2 <- cbind(gene.cds2, gene.cds2_split1$X2)
names(gene.cds2)[3] <- "yup"

antijoin <- anti_join(gene.cds2, toremove_final, by='yup')
antijoin$yup <- NULL

#Reformat to fasta 
antijoin2 <- antijoin %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "ID",
    values_to = "seq")
antijoin2$ID <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.mrna-transcripts_Illumina_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.mrna-transcripts_Illumina+ONT_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

################################################################################          
# Curate proteins
################################################################################ 

setwd("/genes/functional_annot/Illumina")
setwd("/genes/functional_annot/Illumina_ONT")

gene.cds <- read.csv("Rhizophagus_irregularis_DAOM197198.proteins.fa", header=F)
gene.cds2 <- as.data.frame(matrix(gene.cds$V1, ncol = 2, byrow = TRUE))
#Split fasta header
gene.cds2_split1 <- data.frame(str_split_fixed(gene.cds2$V1, " ", 2))
gene.cds2 <- cbind(gene.cds2, gene.cds2_split1$X2)
names(gene.cds2)[3] <- "yup"

antijoin <- anti_join(gene.cds2, toremove_final, by='yup')
antijoin$yup <- NULL

#Reformat to fasta 
antijoin2 <- antijoin %>% 
  pivot_longer(
    cols = c(V1, V2), 
    names_to = "ID",
    values_to = "seq")
antijoin2$ID <- NULL

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.proteins_Illumina_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

write.table(antijoin2, paste(out_dir,"Rhizophagus_irregularis_DAOM197198.proteins_Illumina+ONT_curated.fa", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)


################################################################################          
# Compare functional annotations (Illumina-based VS Illumina+ONT-based)
################################################################################

#Functional annotation
setwd("/genes/functional_annot/Illumina")
GO_Illumina <- read.table("Rhizophagus_irregularis_DAOM197198.annotations.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
GO_Illumina <- GO_Illumina %>% mutate_all(na_if,"")
length(unique(GO_Illumina$GeneID))
setwd("/genes/functional_annot/Illumina_ONT")
GO_Illumina_ONT <- read.table("Rhizophagus_irregularis_DAOM197198.annotations.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
GO_Illumina_ONT <- GO_Illumina_ONT %>% mutate_all(na_if,"")
length(unique(GO_Illumina_ONT$GeneID))

#Number of genes WITHOUT Product name 
GO_Illumina %>% 
  filter(Product != "hypothetical protein") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
31666-31111
555
GO_Illumina_ONT %>% 
  filter(Product != "hypothetical protein") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
31654-31089
565

#Number of genes with PFAM
GO_Illumina %>% 
  filter(PFAM != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
11616
GO_Illumina_ONT %>% 
  filter(PFAM != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
11660

#Number of genes with InterPro
GO_Illumina %>% 
  filter(InterPro != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
15614
GO_Illumina_ONT %>% 
  filter(InterPro != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
15668

#Number of genes with GO.Terms
GO_Illumina %>% 
  filter(GO.Terms != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
10563
GO_Illumina_ONT %>% 
  filter(GO.Terms != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
10664

#Number of genes with Secreted
GO_Illumina %>% 
  filter(Secreted != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
1181
GO_Illumina_ONT %>% 
  filter(Secreted != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
1192

#Number of genes with CAZyme
GO_Illumina %>% 
  filter(CAZyme != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
175
GO_Illumina_ONT %>% 
  filter(CAZyme != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
175

#Number of genes with antiSMASH
#Number of genes with CAZyme
GO_Illumina %>% 
  filter(antiSMASH != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
60
GO_Illumina_ONT %>% 
  filter(antiSMASH != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
60

#N of genes with any annotation
GO_Illumina %>% 
  filter(PFAM != "NA"| InterPro != "NA" | Secreted != "NA" | CAZyme != "NA" | antiSMASH != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
16330

GO_Illumina_ONT %>% 
  filter(PFAM != "NA"| InterPro != "NA" | Secreted != "NA" | CAZyme != "NA" | antiSMASH != "NA") %>%
  group_by(GeneID) %>%
  summarise(Unique_Elements = n_distinct(GeneID))
16401


# create a df (these numbers exclude isoforms)
N <- c("10563", "11616", "15614", "175", "312", "1181", "60", "16330", "10664", "11660", "15668", "175", "312", "1192", "60", "16410")
Category <- c("GO terms", "PFAM", "InterPro", "CAZYmes", "BUSCO_fungi", "Secretome", "antiSMASH_clusters", "Total_annotated", "GO terms", "PFAM", "InterPro", "CAZYmes", "BUSCO_fungi", "Secretome", "antiSMASH_clusters", "Total_annotated")
Annot <- c("Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT")

df <- data.frame(N, Category, Annot)
df$N <- as.numeric(df$N)
df$Annot <- factor(df$Annot, levels = c("Illumina", "Illumina+ONT"))

ggplot(df, aes(x = Category, y = N, group = Annot)) +
  geom_col(aes(fill = N), width = 0.7) +
  ylim(0,19000) +
  facet_grid(~Annot) +
  theme_classic() +
  geom_text(aes(label=N),position="stack",vjust=-1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of functional annotations") +
  theme(legend.position="none")

################################################################################          
# Compare BUSCOs
################################################################################

#Illumina
C:96.9%[S:89.0%,D:7.9%],F:1.7%,M:1.4%,n:290

281	Complete BUSCOs (C)
258	Complete and single-copy BUSCOs (S)
23	Complete and duplicated BUSCOs (D)
5	Fragmented BUSCOs (F)
4	Missing BUSCOs (M)
290	Total BUSCO groups searched

#Illumina+ONT
C:96.9%[S:89.3%,D:7.6%],F:1.4%,M:1.7%,n:290

281	Complete BUSCOs (C)
259	Complete and single-copy BUSCOs (S)
22	Complete and duplicated BUSCOs (D)
4	Fragmented BUSCOs (F)
5	Missing BUSCOs (M)
290	Total BUSCO groups searched

# create a df 
N <- c("96.9", "89.0", "7.9", "1.7", "1.4", "96.9", "89.3", "7.6", "1.7", "1.4")
Category <- c("Complete", "Complete single copy", "Complete duplicated", "Fragmented","Missing", "Complete", "Complete single copy", "Complete duplicated", "Fragmented", "Missing")
Annot <- c("Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT", "Illumina+ONT")

df <- data.frame(N, Category, Annot)
df$N <- as.numeric(df$N)
df$Annot <- factor(df$Annot, levels = c("Illumina", "Illumina+ONT"))

ggplot(df, aes(x = Category, y = N, group = Annot)) +
  geom_col(aes(fill = N), width = 0.7) +
  ylim(0,100) +
  facet_grid(~Annot) +
  theme_classic() +
  geom_text(aes(label=N),position="stack",vjust=-1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of functional annotations") +
  theme(legend.position="none")






