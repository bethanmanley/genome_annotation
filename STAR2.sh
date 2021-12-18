#SBATCH -J star2
#SBATCH -n 5
#SBATCH -o ONT_Transcriptome_G5Riirr/Illumina_RNA_STAR/softmasked/star.o
#SBATCH -e ONT_Transcriptome_G5Riirr/Illumina_RNA_STAR/softmasked/star.e
#SBATCH -p IACT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bm502@cam.ac.uk

conda activate braker2

WD=$1
GENOME_DIR=$2
GENOME=$3

STAR=/mnt/home1/miska/bm502/Apps/STAR-2.7.6a/bin/Linux_x86_64_static/STAR

# Index genome 

$STAR \
--genomeDir $GENOME_DIR \
--genomeFastaFiles $GENOME \
--runMode genomeGenerate \
--genomeSAindexNbases 12 


# Run STAR

for R in $WD/*_1.fq;
do
    prefix=$(basename $R _1.fq)
    $STAR --genomeDir $GENOME_DIR --readFilesIn "$R" "${R/_1.fq/_2.fq}" --outFilterMultimapNmax 20 --outFileNamePrefix $WD/$prefix ;
done 


