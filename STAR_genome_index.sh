#SBATCH -J star_index
#SBATCH -n 5
#SBATCH -o ONT_Transcriptome_G5Riirr/Illumina_RNA_STAR/softmasked/star.o
#SBATCH -e ONT_Transcriptome_G5Riirr/Illumina_RNA_STAR/softmasked/star.e
#SBATCH -p IACT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bm502@cam.ac.uk

GENOME_DIR=$1
GENOME=$2

Apps/STAR-2.7.6a/bin/Linux_x86_64_static/STAR \
--genomeDir $GENOME_DIR
--genomeFastaFiles $GENOME
--runMode genomeGenerate \
--genomeSAindexNbases 12 

