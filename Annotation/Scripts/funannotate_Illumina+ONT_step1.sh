#!/bin/bash

#SBATCH -J step1
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH -n 16
#SBATCH -o fun.postupdate.step1.out
#SBATCH -e fun.postupdate.step1.err
#SBATCH --mail-user=ad984@cam.ac.uk

# Funannotate step1 - Look up Interpro domains

# Remove asterisks from protein protein sequences if needed.

#InterProScan-5.55-88.0
#with options '-goterms -iprlookup'
./interproscan.sh -cpu 16 -i Rhizophagus_irregularis_DAOM197198.proteins.fa -f XML -goterms -iprlookup -pa -d funct_annot_postupdate/
