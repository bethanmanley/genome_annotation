#!/bin/bash

#SBATCH -J rircur2
#SBATCH --mail-type=ALL
#SBATCH -p 1804
#SBATCH -n 20
#SBATCH -o curated2_repeat_masking.out
#SBATCH -e curated2_repeat_masking.err
#SBATCH --mail-user=ad984@cam.ac.uk

#conda activate EDTA
RepeatMasker -pa 2 -s -a -inv -dir ./curated2/RepMask_Rir_curated2 -no_is -norna -xsmall -nolow -div 40 -lib Rirr.curated_primary.no_mt.unscrubbed.fa.mod.EDTA.TElib_curated2.fa -cutoff 225 Rirr.curated_primary.no_mt.unscrubbed.fa

cd ~/anaconda3/envs/EDTA/share/RepeatMasker
calcDivergenceFromAlign.pl -s RepMask_Rir_curated2.divsum Rirr.curated_primary.no_mt.unscrubbed.fa.align
cd RepMask_Rir_curated2
tail -n 72 RepMask_Rir_curated2.divsum > RepMask_Rir_curated2.Kimura_distance

cd ~/anaconda3/envs/EDTA/share/RepeatMasker/util
./createRepeatLandscape_mod1.pl -div RepMask_Rir_curated2/RepMask_Rir_curated2.divsum -g 146773001  > RepMask_Rir_curated2_RMRL.html
./createRepeatLandscape_mod2.pl -div RepMask_Rir_curated2/RepMask_Rir_curated2.divsum -g 146773001  > RepMask_Rir_curated2_classonly_RMRL.html
