
#EggnogMapper 1.2.7 with Diamond in more-sensitive mode (corresponds to diamond.2.0.7's very-sensitive mode). Compared to default, which is diamond.2.0.7's sensitive mode)
screen -S egg
conda activate eggnog-mapper
emapper.py -m diamond --sensmode more-sensitive -i Rhizophagus_irregularis_DAOM197198.proteins.fa -o eggnog_vs_annot


#Secondary metabolism and transmembrane domain annot antiSMASH version 6.0.1
https://fungismash.secondarymetabolites.org/
antiSMASH 6.0: improving cluster detection and comparison capabilities
Kai Blin, Simon Shaw, Alexander M Kloosterman, Zach Charlop-Powers, Gilles P van Weezel, Marnix H Medema, & Tilmann Weber
Nucleic Acids Research (2021) doi: 10.1093/nar/gkab335.

Submit genome sequence (Rirr.curated_primary.no_mt.unscrubbed.fa).
Download .gbk file
