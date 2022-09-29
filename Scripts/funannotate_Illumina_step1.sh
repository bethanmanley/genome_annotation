#!/bin/bash

#SBATCH -J step1
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH -o fun.step1.out
#SBATCH -e fun.step1.err
#SBATCH --mail-user=ad984@cam.ac.uk

# Funannotate step1 - Look up Interpro domains

# First remove asterisk from protein protein sequences

#InterProScan-5.55-88.0
#with options '-goterms -iprlookup'
./interproscan.sh -i augustus.hints.aa.fa -f XML -goterms -iprlookup -pa -d funct_annot_r
