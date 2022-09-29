library(dplyr)
library(stringr)
library(ggplot2)
require(scales)
library(tidyr)

gene.gff <- as.data.frame(import.gff("Rhizophagus_irregularis_DAOM197198_Illumina+ONT_curated.gff3"))
gene <- gene.gff[gene.gff$type == "gene", ] #31652
gene$length = gene$end - gene$start
mRNA <- gene.gff[gene.gff$type == "mRNA", ] #32706
mRNA$length = mRNA$end - mRNA$start

#5'UTR
fiveprime <- gene.gff[gene.gff$type == "five_prime_UTR", ] #7604
fiveprime$length = fiveprime$end - fiveprime$start
mean(fiveprime$length) #116.0079 nt

#Plot histogram of 5'UTR tail length
#Log10 scale
ggplot(fiveprime, aes(x=length)) + 
  geom_density(color = 4, fill = 4, alpha = 0.25) +
  theme_classic() +
  ggtitle("5'UTR length distribution") +
  xlab("log10(5'UTR length)") +
  ylab("Density") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b")  +
  geom_vline(aes(xintercept = 116.0079), linetype = "dashed")

#3'UTR
threeprime <- gene.gff[gene.gff$type == "three_prime_UTR", ] #9674
threeprime$length = threeprime$end - threeprime$start
mean(threeprime$length) #249.5724 nt

ggplot(threeprime, aes(x=length)) + 
  geom_density(color = 4, fill = 4, alpha = 0.25) +
  theme_classic() +
  ggtitle("3'UTR length distribution") +
  xlab("log10(3'UTR length)") +
  ylab("Density") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b")  +
  geom_vline(aes(xintercept = 249.5724), linetype = "dashed")

