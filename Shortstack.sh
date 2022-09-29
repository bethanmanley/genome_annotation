#!/bin/bash

#SBATCH -J shortstack_man22
#SBATCH --mail-type=ALL
#SBATCH -N 1
#SBATCH -n 8 # Number of cores = cpu
#SBATCH -t 0-24:00
#SBATCH -p IACT
#SBATCH -o shortstack_man22_09-02-22.out
#SBATCH -e shortstack_man22_09-02-22.err
#SBATCH --mail-user=ad984@cam.ac.uk

#Prepped small RNA fastq files from Dallaire et al 2021: column-purifed and oxidised samples
#cat col-2_rRNAdep_dcarotadep.fq col-3_rRNAdep_dcarotadep.fq c-4-ox_rRNAdep_dcarotadep.fq 48-4-ox_rRNAdep_dcarotadep.fq > Rhiir_sRNA_mix3.fq

#Count number of reads in small RNA fastq file (the one curated for running Shortstack in Dallaire et al 2021). i.e. rRNA and carrot filtered out
#echo $(cat Rhiir_sRNA_mix3.fq|wc -l)/4|bc
#70,956,710 reads

#conda activate shortstack

#Run for Manley_2022 assembly. -mmap option : random (random for multi-mappers)
ShortStack --readfile Rhiir_sRNA_mix3.fq  --genomefile Rirr.curated_primary.no_mt.unscrubbed.fa --dicermin 20 --dicermax 27 --foldsize 300 --pad 200 --mincov 10.0rpm --strand_cutoff 0.8 --mmap r

#Make chromosome file for plotting
#samtools faidx Rirr.curated_primary.no_mt.unscrubbed.fa
#cut -f1,2 Rirr.curated_primary.no_mt.unscrubbed.fa.fai > rirman22.chromosome.file
