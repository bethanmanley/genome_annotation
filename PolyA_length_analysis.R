library(dplyr)
library(stringr)
library(ggplot2)
require(scales)
library(tidyr)

#All polyA tails of all transcripts detected (not clustered by 'gene' yet)
file <- read.table("polyAsite_analysis.out", header = F, sep="\t", stringsAsFactors=FALSE, quote="")

#Keep only rows that start with stretch of white space (trimmedSeq)
tails <- file %>% 
  filter(grepl('          ', V1))

#Calculate number of As in each tail
tails$number.of.A <- str_count(tails$V1, "A")
tails$number.of.char <- nchar(tails$V1)-10
tails$perc.A <- (tails$number.of.A)/(tails$number.of.char)*100

#Keep only rows that have >80%A (this number is based on PASA pipeline)
tails_filt <- tails %>% 
  filter(perc.A >= 80)

nrow(tails_filt)
nrow(tails)
1155347/1262165*100
#91.5% transcripts retained after filtering for %A

#Plot histogram of polyA tail length!
ggplot(tails_filt, aes(x=number.of.char)) + 
  geom_density(color = 4, fill = 4, alpha = 0.25) +
  theme_classic() +
  ggtitle("poly(A) tail length distribution") +
  xlab("poly(A) tail length") +
  ylab("Density")

#Log10 scale
ggplot(tails_filt, aes(x=number.of.char)) + 
  geom_density(color = 4, fill = 4, alpha = 0.25) +
  theme_classic() +
  ggtitle("poly(A) tail length distribution") +
  xlab("log10(poly(A) tail length)") +
  ylab("Density") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b")  +
  geom_vline(aes(xintercept = 42.447), linetype = "dashed") 

tails_filt$lognchar = log10(tails_filt$number.of.char)
mean(tails_filt$lognchar)
#mean = 1.556
mean(tails_filt$number.of.char)
#42.4447

# Remove rows with the following strings: this evens out the number of rows per transcript
file_filt1 <- file %>% 
  filter(!grepl('OK polyA site candidate.', V1))
file_filt2 <- file_filt1 %>% 
  filter(!grepl('*polyA site', V1))
#Remove first row
file_filt2 <- file_filt2 %>% slice(-1)
file_filt3 <- file_filt2 %>% 
  filter(!grepl('IGNORE', V1))
file_filt4 <- file_filt3 %>% 
  filter(!grepl('gneomic_seq_part', V1))
#Shape into df
df <- as.data.frame(matrix(file_filt4$V1, ncol = 8, byrow = TRUE))
# Remove string from V3
df$V3 <- gsub("PercentA_trimmedSeq: ", "", df$V3)
#Again keep only rows with V3 >=80
df_filt <- df %>% 
  filter(V3 >= 80)

# MEME ANALYSIS FOR POLYA SIGNAL
file2 <- read.table("Rhizophagus_irregularis_DAOM197198_pasa.polyAsites.fasta", header = F, sep="\t", stringsAsFactors=FALSE, quote="")
file2df <- as.data.frame(matrix(file2$V1, ncol = 2, byrow = TRUE))
# file2df -> summarizes all mapped polyA sites supported by the transcripts
#Split at upper/lowercase. Upper=3'end, last bp is the bd to which first A of the polyA tail is added. Lower=bit of genome after the gene
file2df$UPPER <- sub("([A-Z]+)[a-z].*", "\\1", file2df$V2)
file2df$LOWER <- str_extract_all(file2df$V2, '[[:lower:]]+')
file2df$V2 <- NULL
file2df$LOWER <- NULL

file2df$V1 <- gsub(":.*","", file2df$V1)

file2df_fasta <- file2df %>% 
  pivot_longer(
    cols = c(V1, UPPER), 
    names_to = "TE",
    values_to = "seq")
file2df_fasta$TE <- NULL
write.table(file2df_fasta, paste(out_dir,"for_polyAsignal_search_clusteredtranscripts.fa",sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

meme for_polyAsignal_search_clusteredtranscripts.fa -dna -oc polyAmeme_clustered -nostatus -time 14400 -mod zoops -nmotifs 18 -minw 6 -maxw 6 -objfun classic -markov_order 0

#Hierarchically search for derivatives of the canonical polyA signal
library(data.table)
AATAAA <- file2df[file2df$UPPER %like% "AATAAA", ] #
AATAAA <- file2df %>% dplyr::filter(stringr::str_detect(UPPER, "AATAAA", negate=T))

AUUAAA <- AATAAA[AATAAA$UPPER %like% "ATTAAA", ] #28045
AUUAAA <- AATAAA %>% dplyr::filter(stringr::str_detect(UPPER, "ATTAAA", negate=T))

AAUAUA <- AUUAAA[AUUAAA$UPPER %like% "AATATA", ] #11489
AAUAUA <- AUUAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AAtAtA", negate=T))

TATAAA <- AAUAUA[AAUAUA$UPPER %like% "TATAAA", ] #4581
TATAAA <- AAUAUA %>% dplyr::filter(stringr::str_detect(UPPER, "TATAAA", negate=T))

AATAAT <- TATAAA[TATAAA$UPPER %like% "AATAAT", ] #1876
AATAAT <- TATAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AATAAT", negate=T))

AGTAAA <- TATAAA[TATAAA$UPPER %like% "AGTAAA", ] #517
AGTAAA <- TATAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AGTAAA", negate=T))

AACAAA <- AGTAAA[AGTAAA$UPPER %like% "AACAAA", ] #431
AACAAA <- AGTAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AACAAA", negate=T))

AAGAAA <- AACAAA[AACAAA$UPPER %like% "AAGAAA", ] #481
AAGAAA <- AACAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AAGAAA", negate=T))

AATACA <- AAGAAA[AAGAAA$UPPER %like% "AATACA", ] #431
AATACA <- AAGAAA %>% dplyr::filter(stringr::str_detect(UPPER, "AATACA", negate=T))

CATAAA <- AATACA[AATACA$UPPER %like% "CATAAA", ] #212
CATAAA <- AATACA %>% dplyr::filter(stringr::str_detect(UPPER, "CATAAA", negate=T))

GATAAA <- CATAAA[CATAAA$UPPER %like% "GATAAA", ] #212
GATAAA <- CATAAA %>% dplyr::filter(stringr::str_detect(UPPER, "GATAAA", negate=T))

TATGAA <- GATAAA[GATAAA$UPPER %like% "TATGAA", ] #336
TATGAA <- GATAAA %>% dplyr::filter(stringr::str_detect(UPPER, "TATGAA", negate=T))

AATAAC <- TATGAA[TATGAA$UPPER %like% "AATAAC", ] #196
AATAAC <- TATGAA %>% dplyr::filter(stringr::str_detect(UPPER, "AATAAC", negate=T))

AATAGA <- AATAAC[AATAAC$UPPER %like% "AATAGA", ] #266
AATAGA <- AATAAC %>% dplyr::filter(stringr::str_detect(UPPER, "AATAGA", negate=T))

AATAAG <- AATAGA[AATAGA$UPPER %like% "AATAAG", ] #117
AATAAG <- AATAGA %>% dplyr::filter(stringr::str_detect(UPPER, "AATAAG", negate=T))

ACTAAA <- AATAAG[AATAAG$UPPER %like% "ACTAAA", ] #169
ACTAAA <- AATAAG %>% dplyr::filter(stringr::str_detect(UPPER, "ACTAAA", negate=T))

GATGAA <- ACTAAA[ACTAAA$UPPER %like% "GATGAA", ] #169
GATGAA <- ACTAAA %>% dplyr::filter(stringr::str_detect(UPPER, "GATGAA", negate=T))

CATGAA <- GATGAA[GATGAA$UPPER %like% "CATGAA", ] #169
CATGAA <- GATGAA %>% dplyr::filter(stringr::str_detect(UPPER, "CATGAA", negate=T))



