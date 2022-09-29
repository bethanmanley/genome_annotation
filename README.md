# genome_annotation
Series of scripts associated with the genome annotation of Rhizophagus irregularis. 

# Annotation using BRAKER2 for Rhizophagus irregularis
# Using BRAKER2 from conda. Create in this order

conda create --name braker2 augustus braker2


## Download GeneMark + GeneMark key
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_5Rxuj/gmes_linux_64.tar.gz

wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_f_2SU/gm_key_64.gz
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
