#!/bin/bash

#SBATCH -J funup6
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH --mem=240G
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -o fun.funupdate6.out
#SBATCH -e fun.funupdate6.err
#SBATCH --mail-user=ad984@cam.ac.uk

#Sub-sample 500k Nanopore reads
#seqtk sample -2 -s100 G5_Riirr_cDNA_trimmed_bcode1.fastq 500000 > G5_Riirr_cDNA_trimmed_bcode1.sub1.fastq

funannotate update -i Rhizophagus_irregularis_DAOM197198.gbk \
--fasta Rirr.curated_primary.no_mt.unscrubbed.fa \
--gff Rhizophagus_irregularis_DAOM197198_curated.gff3 \
--nanopore_cdna G5_Riirr_cDNA_trimmed_bcode1.sub1.fastq \
--species "Rhizophagus irregularis" \
--strain "DAOM197198" \
--jaccard_clip \
--max_intronlen 5000 \
--cpus 20 \
--memory 100G \
-o update6/
