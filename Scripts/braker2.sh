#!/bin/bash
#SBATCH -J braker2
#SBATCH -n 3
#SBATCH -o Guppy_5_Riirr/shasta0.5/braker2/braker.o
#SBATCH -e Guppy_5_Riirr/shasta0.5/braker2/braker.e
#SBATCH -p 1804
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bm502@cam.ac.uk

conda activate braker2

SPECIES=$1
WD=$2
GENOME=$3
BAM=$4
CORES=$5


mkdir -p $WD

miniconda3/bin/braker.pl \
--workingdir=$WD \
--genome=$GENOME \
--GENEMARK_PATH=/mnt/home1/miska/bm502/Apps/gmes_linux_64 \
--fungus \
--softmasking \
--gff3 \
--bam=$BAM \
--cores=$CORES
