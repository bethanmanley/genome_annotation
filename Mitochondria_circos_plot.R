library(circlize)
library(dplyr)
library("GenomicRanges")
library("GenomicFeatures")
library(rtracklayer)

# Make chromosome file
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' Rhizophagus_irregularis_DAOM197198_mtDNA.fasta
contig <- c("mitochondria")
chromStart <- c("0")
chromEnd <- c("70793")
gieStain <- c("gneg")
chr_file <- data.frame(contig, chromStart, chromEnd, gieStain)

#Mito gene annotation
gff <- data.frame(import.gff("Rhizophagus_irregularis_DAOM197198_mtDNA.gff"))
colnames(gff) <- c("Seqname", "Start", "End", "width",  "Strand", "Source", "Type", "Score", "Phase", "ID") 
gff$Seqname <- "mitochondria"
gff$value <- as.numeric("1")
genes.file <- data.frame(gff$Seqname, gff$Start, gff$End, gff$value, gff$ID, gff$Type)
colnames(genes.file) <- c("Seqname", "Start", "End", "value", "ID", "Type")
#Simplify IDs
genes.file$ID <- sub(" .*", "", genes.file$ID)
#Isolate tRNA and CDS rows for plotting
target <- c("CDS", "tRNA", "rRNA")
genes.file.plot <- filter(genes.file, Type %in% target)


sapply(genes.file, class)

#Initialize plot
circos.clear()
circos.par("start.degree" = 90, "track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(chr_file)
#circos.genomicDensity(genes.file, col = c("#0000FF80"), track.height = 0.1)
circos.genomicLabels(genes.file.plot, labels.column=5, side = "inside",
                     col = as.numeric(factor(genes.file.plot[[6]])), line_col = as.numeric(factor(genes.file.plot[[6]])))


