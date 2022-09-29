#!/bin/bash

#SBATCH -J step3
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH -o fun.step3.out
#SBATCH -e fun.step3.err
#SBATCH --mail-user=ad984@cam.ac.uk

#Setup busco database for fungi_odb9 (default was dikarya)
#cd ~/funannotate_db
#funannotate setup -d ~/funannotate_db -b fungi

#install signalp 5.0. Download from https://services.healthtech.dtu.dk/service.php?SignalP-5.0
#untar on local machine and transfer the signalp bin file (33Mb) to ~/anaconda3/envs/funannotate/bin, and transfer the signalp lb files into anaconda3/envs/funannotate/lib
#check that it's installed: funannotate check --show-versions # should see 'signalp: 5.0b' in list of external dependencies

funannotate annotate --gff braker.gff3 --fasta Rirr.curated_primary.no_mt.unscrubbed.fa.classified.hs.masked.fa -s "Rhizophagus irregularis" -o funct_annot_r/ \
--sbt Genbanktemplate.sbt --eggnog eggnog_vs_annot.emapper.annotations --antismash Rirr.curated_primary.no_mt.unscrubbed.gbk \
--iprscan augustus.hints.aa.fa.xml  \
--isolate DAOM197198 --busco_db fungi --cpus 2
