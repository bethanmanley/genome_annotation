library('ggplot2')
library(stringr)
library("GenomicRanges")
library("GenomicFeatures")
library(rtracklayer)
library(chromPlot)
library(dplyr)
library(tidytext)
library(textrank)
library(ggrepel)



setwd("/Users/alexandradallaire/genomics/Manley_2022/transcript_assembly/genes")
gene.gff <- as.data.frame(import.gff("Rhizophagus_irregularis_DAOM197198_curated.gtf"))
x_genes <- gene.gff[gene.gff$type == "transcript", ]
x_genes$length = x_genes$end - x_genes$start

GO <- read.table("Rhizophagus_irregularis_DAOM197198.annotations.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
GO <- GO %>% mutate_all(na_if,"")

setwd("/Users/alexandradallaire/genomics/predsym/hourglass/Rirr_PS")
PS <- read.csv(file = 'PS_filt.csv', head=TRUE, sep=";")
names(PS)[1] <- "transcript_id"

PS_j <- left_join(PS, x_genes, by=c("transcript_id"))

#Add Gene annotation
names(GO)[2] <- "transcript_id"
GO$transcript_id <- gsub("\\-T",".t", GO$transcript_id)

PS_j_info <- left_join(PS_j, GO, by=c("transcript_id"))
PS_j_info <- PS_j_info[!is.na(PS_j_info$seqnames),]


head(PS_j_info)
colnames(PS_j_info)
PS_j_info$gDNA <- NULL
PS_j_info$mRNA <- NULL
PS_j_info$CDS.transcript <- NULL
PS_j_info$Translation <- NULL

################################################################################          
#### INTERPRO DOMAIN
################################################################################          
#Extract words from InterPro column
Data <- PS_j_info %>% 
  group_by(PS) %>% 
  mutate(InterPro = paste0(InterPro, collapse = ";")) %>%
  select(PS, InterPro) %>% 
  distinct()

#Tokenise the InterPro words
myTokens <- Data %>% 
  unnest_tokens(word, InterPro) %>% 
  anti_join(stop_words, by = "word")

#Keep only words that start with 'ipr'
myTokens <- myTokens %>%
  filter(str_detect(word, "^ipr"))

#Count occurence of each IPR domain, grouped by PS
counts <- myTokens %>% 
  group_by(PS, word) %>% 
  summarize(count=n())

#Keep top 5 for each PS group
top5 <- counts %>%
  arrange(desc(count)) %>%
  group_by(PS) %>% 
  slice(1:5)

#Plot 1
# top 3
top5 <- read.table("top5.PS.IPRs_curated.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
top3 <- top5 %>%
  arrange(desc(count)) %>%
  group_by(PS) %>% 
  slice(1:3)

ggplot(top3, aes(x = PS, y = log(count), fill = word, label = word)) +
  geom_bar(position="stack", stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), angle = 30) +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip()

ggplot(top3 %>% group_by(PS) %>% arrange(count, by_group = TRUE),
       aes(fill = word, y = log(count), x = PS, group = PS, label=id)) + 
  geom_col(position = "stack") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), angle = 30) +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip() +
  xlab("Phylostratum") +
  ylab("log(count)")
