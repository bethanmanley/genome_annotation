#!/bin/bash

#SBATCH -J rhiir_22
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH -o edta_man2022_31-01-22.out
#SBATCH -e edta_man2022_31-01-22.err
#SBATCH --mail-user=ad984@cam.ac.uk

perl ~/EDTA/EDTA.pl --genome Rirr.curated_primary.no_mt.unscrubbed.fa --species others --step all --overwrite 0 --sensitive 1 --anno 1 --evaluate 0 --threads 42
