
samtools faidx Rirr.curated_primary.no_mt.unscrubbed.fa
cut -f1,2 Rirr.curated_primary.no_mt.unscrubbed.fa.fai > Rhiir_chrfile

cd genomics/rirman22/intron
Sort the GFF files:
cat PS1_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > PS1_curated_s.gff3
cat PS2_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS2_curated_s.gff3
cat PS3_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS3_curated_s.gff3
cat PS4_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS4_curated_s.gff3
cat PS5_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS5_curated_s.gff3
cat PS6_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS6_curated_s.gff3
cat PS7_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS7_curated_s.gff3
cat PS8_curated.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k4,4n -k5,5rn -k3,3r"}' > PS8_curated_s.gff3

#Split gffs into pos and neg strand
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS1_curated_s.gff3 > PS1_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS1_curated_s.gff3 > PS1_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS2_curated_s.gff3 > PS2_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS2_curated_s.gff3 > PS2_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS3_curated_s.gff3 > PS3_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS3_curated_s.gff3 > PS3_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS4_curated_s.gff3 > PS4_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS4_curated_s.gff3 > PS4_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS5_curated_s.gff3 > PS5_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS5_curated_s.gff3 > PS5_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS6_curated_s.gff3 > PS6_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS6_curated_s.gff3 > PS6_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS7_curated_s.gff3 > PS7_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS7_curated_s.gff3 > PS7_curated_s_pos.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "-") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS8_curated_s.gff3 > PS8_curated_s_neg.gff3
awk 'BEGIN{OFS="\t"} {if ($7 == "+") print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' PS8_curated_s.gff3 > PS8_curated_s_pos.gff3

#Get intergenic regions
bedtools complement -i PS1_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS1_intergenic_sorted_neg.bed
bedtools complement -i PS1_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS1_intergenic_sorted_pos.bed

bedtools complement -i PS2_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS2_intergenic_sorted_neg.bed
bedtools complement -i PS2_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS2_intergenic_sorted_pos.bed

bedtools complement -i <(cat PS3_curated_s_neg.gff3 | sort -k1,1V -k4,4n -k5,5rn -k3,3r) -g ~/genomics/rirman22/Rhiir_chrfile > PS3_intergenic_sorted_neg.bed
bedtools complement -i PS3_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS3_intergenic_sorted_pos.bed

bedtools complement -i PS4_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS4_intergenic_sorted_neg.bed
bedtools complement -i PS4_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS4_intergenic_sorted_pos.bed

bedtools complement -i PS5_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS5_intergenic_sorted_neg.bed
bedtools complement -i PS5_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS5_intergenic_sorted_pos.bed

bedtools complement -i PS6_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS6_intergenic_sorted_neg.bed
bedtools complement -i PS6_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS6_intergenic_sorted_pos.bed

bedtools complement -i PS7_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS7_intergenic_sorted_neg.bed
bedtools complement -i PS7_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS7_intergenic_sorted_pos.bed

bedtools complement -i PS8_curated_s_neg.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS8_intergenic_sorted_neg.bed
bedtools complement -i PS8_curated_s_pos.gff3 -g ~/genomics/rirman22/Rhiir_chrfile > PS8_intergenic_sorted_pos.bed




Next, intron is the complement of intergenic and exonic regions. Extract exonic coordinates in BED format:
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS1_curated_s_neg.gff3 > PS1_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS1_curated_s_pos.gff3 > PS1_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS2_curated_s_neg.gff3 > PS2_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS2_curated_s_pos.gff3 > PS2_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS3_curated_s_neg.gff3 | sort -k1,1 -k2,2n > PS3_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS3_curated_s_pos.gff3 | sort -k1,1 -k2,2n > PS3_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS4_curated_s_neg.gff3 > PS4_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS4_curated_s_pos.gff3 > PS4_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS5_curated_s_neg.gff3 > PS5_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS5_curated_s_pos.gff3 > PS5_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS6_curated_s_neg.gff3 > PS6_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS6_curated_s_pos.gff3 > PS6_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS7_curated_s_neg.gff3 > PS7_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS7_curated_s_pos.gff3 > PS7_exon_sorted_pos.bed

awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS8_curated_s_neg.gff3 > PS8_exon_sorted_neg.bed
awk 'BEGIN{OFS="\t"} {if ($3 == "exon") print $1, $4-1, $5}' PS8_curated_s_pos.gff3 > PS8_exon_sorted_pos.bed


# Use BEDtools complement to get the introns:
bedtools complement -i <(cat PS1_exon_sorted_neg.bed PS1_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS1_intron_sorted_neg.bed
bedtools complement -i <(cat PS1_exon_sorted_pos.bed PS1_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS1_intron_sorted_pos.bed

bedtools complement -i <(cat PS2_exon_sorted_neg.bed PS2_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS2_intron_sorted_neg.bed
bedtools complement -i <(cat PS2_exon_sorted_pos.bed PS2_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS2_intron_sorted_pos.bed

bedtools complement -i <(cat PS3_exon_sorted_neg.bed PS3_intergenic_sorted_neg.bed | sort -V -k1,1 -k2,2) -g ~/genomics/rirman22/Rhiir_chrfile > PS3_intron_sorted_neg.bed
bedtools complement -i <(cat PS3_exon_sorted_pos.bed PS3_intergenic_sorted_pos.bed | sort -V -k1,1 -k2,2) -g ~/genomics/rirman22/Rhiir_chrfile > PS3_intron_sorted_pos.bed

bedtools complement -i <(cat PS4_exon_sorted_neg.bed PS4_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS4_intron_sorted_neg.bed
bedtools complement -i <(cat PS4_exon_sorted_pos.bed PS4_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS4_intron_sorted_pos.bed

bedtools complement -i <(cat PS5_exon_sorted_neg.bed PS5_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS5_intron_sorted_neg.bed
bedtools complement -i <(cat PS5_exon_sorted_pos.bed PS5_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS5_intron_sorted_pos.bed

bedtools complement -i <(cat PS6_exon_sorted_neg.bed PS6_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS6_intron_sorted_neg.bed
bedtools complement -i <(cat PS6_exon_sorted_pos.bed PS6_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS6_intron_sorted_pos.bed

bedtools complement -i <(cat PS7_exon_sorted_neg.bed PS7_intergenic_sorted_neg.bed | sort -V -k1,1 -k2,2) -g ~/genomics/rirman22/Rhiir_chrfile > PS7_intron_sorted_neg.bed
bedtools complement -i <(cat PS7_exon_sorted_pos.bed PS7_intergenic_sorted_pos.bed | sort -V -k1,1 -k2,2) -g ~/genomics/rirman22/Rhiir_chrfile > PS7_intron_sorted_pos.bed

bedtools complement -i <(cat PS8_exon_sorted_neg.bed PS8_intergenic_sorted_neg.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS8_intron_sorted_neg.bed
bedtools complement -i <(cat PS8_exon_sorted_pos.bed PS8_intergenic_sorted_pos.bed | sort -k1,1 -k2,2n) -g ~/genomics/rirman22/Rhiir_chrfile > PS8_intron_sorted_pos.bed


bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS1_intron_sorted_neg.bed > PS1_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS1_intron_sorted_pos.bed > PS1_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS2_intron_sorted_neg.bed > PS2_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS2_intron_sorted_pos.bed > PS2_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS3_intron_sorted_neg.bed > PS3_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS3_intron_sorted_pos.bed > PS3_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS4_intron_sorted_neg.bed > PS4_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS4_intron_sorted_pos.bed > PS4_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS5_intron_sorted_neg.bed > PS5_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS5_intron_sorted_pos.bed > PS5_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS6_intron_sorted_neg.bed > PS6_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS6_intron_sorted_pos.bed > PS6_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS7_intron_sorted_neg.bed > PS7_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS7_intron_sorted_pos.bed > PS7_intron_sorted_pos.fa

bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS8_intron_sorted_neg.bed > PS8_intron_sorted_neg.fa
bedtools getfasta -fi ~/genomics/rirman22/Rirr.curated_primary.no_mt.unscrubbed.fa -bed PS8_intron_sorted_pos.bed > PS8_intron_sorted_pos.fa

#Reverse complement the negative strand
revseq PS1_intron_sorted_neg.fa -reverse -complement -outseq PS1_intron_sorted_neg_revcomp.fa
revseq PS2_intron_sorted_neg.fa -reverse -complement -outseq PS2_intron_sorted_neg_revcomp.fa
revseq PS3_intron_sorted_neg.fa -reverse -complement -outseq PS3_intron_sorted_neg_revcomp.fa
revseq PS4_intron_sorted_neg.fa -reverse -complement -outseq PS4_intron_sorted_neg_revcomp.fa
revseq PS5_intron_sorted_neg.fa -reverse -complement -outseq PS5_intron_sorted_neg_revcomp.fa
revseq PS6_intron_sorted_neg.fa -reverse -complement -outseq PS6_intron_sorted_neg_revcomp.fa
revseq PS7_intron_sorted_neg.fa -reverse -complement -outseq PS7_intron_sorted_neg_revcomp.fa
revseq PS8_intron_sorted_neg.fa -reverse -complement -outseq PS8_intron_sorted_neg_revcomp.fa


cat PS1_intron_sorted_neg_revcomp.fa PS1_intron_sorted_pos.fa > PS1_intron.fa
cat PS2_intron_sorted_neg_revcomp.fa PS2_intron_sorted_pos.fa > PS2_intron.fa
cat PS3_intron_sorted_neg_revcomp.fa PS3_intron_sorted_pos.fa > PS3_intron.fa
cat PS4_intron_sorted_neg_revcomp.fa PS4_intron_sorted_pos.fa > PS4_intron.fa
cat PS5_intron_sorted_neg_revcomp.fa PS5_intron_sorted_pos.fa > PS5_intron.fa
cat PS6_intron_sorted_neg_revcomp.fa PS6_intron_sorted_pos.fa > PS6_intron.fa
cat PS7_intron_sorted_neg_revcomp.fa PS7_intron_sorted_pos.fa > PS7_intron.fa
cat PS8_intron_sorted_neg_revcomp.fa PS8_intron_sorted_pos.fa > PS8_intron.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS1_intron.fa > PS1_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS2_intron.fa > PS2_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS3_intron.fa > PS3_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS4_intron.fa > PS4_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS5_intron.fa > PS5_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS6_intron.fa > PS6_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS7_intron.fa > PS7_intron_m.fa

awk '!/^>/ { printf "%s", $0; n = "\n" }
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' PS8_intron.fa > PS8_intron_m.fa

#Extract first 12 nt at the beginning of introns
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS1_intron_m.fa > PS1_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS2_intron_m.fa > PS2_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS3_intron_m.fa > PS3_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS4_intron_m.fa > PS4_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS5_intron_m.fa > PS5_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS6_intron_m.fa > PS6_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS7_intron_m.fa > PS7_intron_first12.fa
awk -v OFS="\n" '{getline seq} {print $0, substr(seq,1,12)}' PS8_intron_m.fa > PS8_intron_first12.fa

#Clean
seqtk seq -L 8 PS1_intron_first12.fa > PS1_intron_first12_f.fa
seqtk seq -L 8 PS2_intron_first12.fa > PS2_intron_first12_f.fa
seqtk seq -L 8 PS3_intron_first12.fa > PS3_intron_first12_f.fa
seqtk seq -L 8 PS4_intron_first12.fa > PS4_intron_first12_f.fa
seqtk seq -L 8 PS5_intron_first12.fa > PS5_intron_first12_f.fa
seqtk seq -L 8 PS6_intron_first12.fa > PS6_intron_first12_f.fa
seqtk seq -L 8 PS7_intron_first12.fa > PS7_intron_first12_f.fa
seqtk seq -L 8 PS8_intron_first12.fa > PS8_intron_first12_f.fa

meme PS1_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS2_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS3_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS4_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS5_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS6_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS7_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS8_intron_first12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0


#Last 12 nt
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS1_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS1_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS2_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS2_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS3_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS3_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS4_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS4_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS5_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS5_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS6_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS6_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS7_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS7_intron_last12.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' PS8_intron_m.fa |\
awk -F '\t' '{x=12;L=length($2);printf("%s\n%s\n",$1,(L<=x?$2:substr($2,1+L-x,x)));}' > PS8_intron_last12.fa

seqtk seq -L 8 PS1_intron_last12.fa > PS1_intron_last12_f.fa
seqtk seq -L 8 PS2_intron_last12.fa > PS2_intron_last12_f.fa
seqtk seq -L 8 PS3_intron_last12.fa > PS3_intron_last12_f.fa
seqtk seq -L 8 PS4_intron_last12.fa > PS4_intron_last12_f.fa
seqtk seq -L 8 PS5_intron_last12.fa > PS5_intron_last12_f.fa
seqtk seq -L 8 PS6_intron_last12.fa > PS6_intron_last12_f.fa
seqtk seq -L 8 PS7_intron_last12.fa > PS7_intron_last12_f.fa
seqtk seq -L 8 PS8_intron_last12.fa > PS8_intron_last12_f.fa

meme PS1_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS2_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS3_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS4_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS5_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS6_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 8 -minw 4 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS7_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
meme PS8_intron_last12_f.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 4 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0
