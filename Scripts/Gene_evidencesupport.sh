
################################################################################################
# Check how well the short-reads support the short-read-guided BRAKER2-gene annotation
################################################################################################
./selectSupportedSubsets.py augustus.hints.gtf genemark_hintsfile.gff --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT
./generateReport.py augustus.hints.gtf genemark_hintsfile.gff BRAKER_prediction_report.pdf

Genes fully supported by external evidence: 9848 (31.1%)
Genes partially supported by external evidence: 15232 (48.1%)
Genes unsupported by any external evidence: 16435 (51.9%)
***Careful: all single-exon/intronless genes are *un*supported (the hints comes from introns)

Gene count: 31667 Single-exon genes: 8841 Multi-exon genes: 22826
Introns per gene: 2.98
Introns per multi-exon gene: 4.14

################################################################################################
# Check how well the long-reads support the short-read-guided BRAKER2-gene annotation
################################################################################################
#Generate long-read-hints (bam2hints is part of augustus Trinity)
source ~/anaconda3/etc/profile.d/conda.sh
conda activate funannotate

#Run funannotate predict to generate long-read hints file
funannotate predict -i Rirr.curated_primary.no_mt.unscrubbed.fa -o predict \
      --species "Rhizophagus irregularis" --strain DAOM197198 \
      --transcript_evidence long-reads.mapped.fasta \
      --rna_bam nano_cDNA.coordSorted.bam \
      --cpus 12

#Let run for ~45min, keep only hints.BAM.gff (intron evidence)
#Use hints file
./selectSupportedSubsets.py augustus.hints.gtf hints.BAM.gff --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT
./generateReport.py augustus.hints.gtf hints.BAM.gff long-read-hints_short-read-annot_prediction_report.pdf

Genes fully supported by external evidence: 8294 (26.19%)
Genes partially supported by external evidence: 12976 (40.98%)
Genes unsupported by any external evidence: 18691 (59.02%)


################################################################################################
# Check how well the long-reads support the PASA-long-read-guided gene annotation
################################################################################################
#Convert gff3 to gtf. Run on own machine [../transcript_assembly/genes/update6]
gffread -E Rhizophagus_irregularis_DAOM197198.gff3 -T -o Rhizophagus_irregularis_DAOM197198.gtf
#Run on own machine [../transcript_assembly/genes]
./selectSupportedSubsets.py update6/Rhizophagus_irregularis_DAOM197198.gtf hints.BAM.gff --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT

################################################################################################
# Check how well the short-reads support the PASA-long-read-guided gene annotation
################################################################################################
#Run on own machine [../transcript_assembly/genes]
./selectSupportedSubsets.py update6/Rhizophagus_irregularis_DAOM197198.gtf genemark_hintsfile.gff --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT
