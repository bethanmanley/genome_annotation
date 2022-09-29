
#Make TSV file for g::profiler with format:
#"GO TERM" , “GENE PRODUCT ID” and “GO NAME”  correspond to the annotation, ID and description columns.
# Example:
#GO:0003677  A0A015I1C0 DNA binding
annot <- read.table("all.annotations.txt", header = F, sep="\t", stringsAsFactors=FALSE, quote="")
annot <- annot %>% mutate_all(na_if,"")

unique(annot$V2)
#Keep go_function, go_component and go_process rows
GO <- dplyr::filter(annot, grepl("go_function|go_component|go_process", V2, ignore.case=F))

#Split GO NAME column into 3 parts 
GO.split <- data.frame(str_split_fixed(GO$V3, "\\|", 3))
#Add string to GO TERMs
string <- "GO:"
GO.split$X4 <- paste(string, GO.split$X2, sep="")

#Make df
GO_df <- cbind(GO.split$X4, GO$V1, GO.split$X1)

write.table(GO_df, paste(out_dir,"Rhizophagus_irregularis_DAOM197198_GOterms.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

#Go to: https://biit.cs.ut.ee/gmt-helper/
#Convert *.tsv to .gmt

#Then go to http://geneontology.org/docs/download-ontology/ and download the latest go-basic.obo file. 
#Upload the GMT file and the obo file to the tab “Reannotate GMT using OBO”. 
#The output is a fully annotated gmt file (Rhizophagus_irregularis_DAOM197198_GOterms_final.gmt)
#There was a problem in g:profiler's GMT processing, which they fixed upon request.
#The fixed/functional .gmt file is then: Rhizophagus_irregularis_DAOM197198_GOterms_final_fixed.gmt 

#This gmt file was validated and can be used with the following token:
gp__xfGY_dQeI_yx4

#Check which GO-terms are enriched at different phyloranks
setwd("/Users/alexandradallaire/genomics/Manley_2022/geneage/Rirr_usen/")
PS <- read_tsv("1432141_phylostrata_assignation_filtered.tsv")

out_dir <- ("/Users/alexandradallaire/genomics/Manley_2022/geneage/Rirr_usen/")
PS9 <- dplyr::filter(PS, grepl("9", PS, ignore.case=F))
write.table(PS9, paste(out_dir,"PS9.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS8 <- dplyr::filter(PS, grepl("8", PS, ignore.case=F))
write.table(PS8, paste(out_dir,"PS8.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS7 <- dplyr::filter(PS, grepl("7", PS, ignore.case=F))
write.table(PS7, paste(out_dir,"PS7.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS6 <- dplyr::filter(PS, grepl("6", PS, ignore.case=F))
write.table(PS6, paste(out_dir,"PS6.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS5 <- dplyr::filter(PS, grepl("5", PS, ignore.case=F))
write.table(PS5, paste(out_dir,"PS5.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS4 <- dplyr::filter(PS, grepl("4", PS, ignore.case=F))
write.table(PS4, paste(out_dir,"PS4.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS3 <- dplyr::filter(PS, grepl("3", PS, ignore.case=F))
write.table(PS3, paste(out_dir,"PS3.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS2 <- dplyr::filter(PS, grepl("2", PS, ignore.case=F))
write.table(PS2, paste(out_dir,"PS2.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS1 <- dplyr::filter(PS, grepl("1", PS, ignore.case=F))
write.table(PS1, paste(out_dir,"PS1.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

#MAKE FINAL TABLE FOR PLOT
enriched_GO <- read.table("Enriched_GOterms.tsv", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
enriched_GO$padj <- gsub("×10", "E", as.character(enriched_GO$padj))
enriched_GO$padj <- as.numeric(enriched_GO$padj)
enriched_GO <- select(enriched_GO, 1:4)
write.table(enriched_GO, paste(out_dir,"Enriched_GOterms_final.tsv", sep="/"), col.names=T, quote=F, sep="\t", row.names=F)









#Check which GO-terms are enriched at different phyloranks, FOR HIGH CONFIDENCE GENES
setwd("/Users/alexandradallaire/genomics/Manley_2022/geneage/final_HDF/")
PS <- read_tsv("1432141_high_confidence_gene_ages_filtered.tsv")

out_dir <- ("/Users/alexandradallaire/genomics/Manley_2022/geneage/final_HDF/")

PS8 <- dplyr::filter(PS, grepl("8", PS, ignore.case=F))
write.table(PS8, paste(out_dir,"PS8.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS7 <- dplyr::filter(PS, grepl("7", PS, ignore.case=F))
write.table(PS7, paste(out_dir,"PS7.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS6 <- dplyr::filter(PS, grepl("6", PS, ignore.case=F))
write.table(PS6, paste(out_dir,"PS6.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS5 <- dplyr::filter(PS, grepl("5", PS, ignore.case=F))
write.table(PS5, paste(out_dir,"PS5.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

PS4 <- dplyr::filter(PS, grepl("4", PS, ignore.case=F))
write.table(PS4, paste(out_dir,"PS4.tsv", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)


#MAKE FINAL TABLE FOR PLOT
enriched_GO <- read.table("Enriched_GOterms.tsv", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
enriched_GO$padj <- gsub("×10", "E", as.character(enriched_GO$padj))
enriched_GO$padj <- as.numeric(enriched_GO$padj)
enriched_GO <- select(enriched_GO, 1:4)
write.table(enriched_GO, paste(out_dir,"Enriched_GOterms_final.tsv", sep="/"), col.names=T, quote=F, sep="\t", row.names=F)


 