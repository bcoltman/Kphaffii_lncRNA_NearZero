#!/usr/bin/env bash
#### Carry out initial lncRNA prediction using feelnc
#### inputs are: 1) reference genome fasta 2) reference genome GTF
#### 3) stringtie GTF 4) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-f = Reference genome sequence "
        echo "-G = Reference genome GTF"
        echo "-g = Stringtie GTF"
	echo "-p = num processors"
        echo "-o = Output directory"
	echo "-c = coding fasta"
	echo "-n = non-coding fasta"
        exit 2
fi
while getopts f:G:g:p:o:c:n: option
  do
    case "${option}"
      in
      f) REF_SEQ=${OPTARG};;
      G) REF_GTF=${OPTARG};;
      g) STR_GTF=${OPTARG};;
      p) THREADS=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
      c) CODING_FA=${OPTARG};;
      n) NONCODING_FA=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

echo $OUT_DIR




#FEELnc_filter.pl \
#-i  transcriptome_assembly/non_protein_coding_stringtie.gtf \
#-a reference_genome/cbs7435_codingtranscripts.gtf \
#-o lncrna_annotation/FEELnc_test/file.log \
#--monoex=1 \
#--size=200 \
#-p 30 \
#> lncrna_annotation/FEELnc_test/candidate_lncRNA.gtf


# filter transcripts overlapping with sense protein coding exons. In our GTF, no protein_coding flag in 
# the GTF file so do not set -b flag for iotype filtering 
#-b transcript_biotype=protein_coding,pseudogene 
# the gffread in the previous step keeps all transcripts that do not have CDS features, 
# this is only known protein coding transcripts
# Keeping monoex to deal with misc_RNA biotype - monoex=1 keeps all 
FEELnc_filter.pl \
-i $STR_GTF \
-a $REF_GTF \
-o $OUT_DIR/file.log \
--monoex=1 \
--size=200 \
-p $THREADS \
> $OUT_DIR/candidate_lncRNA.gtf

# create a gtf for known protein coding transcripts - already done
#awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' \
#$REF_GTF | \
#grep 'protein_coding' \
#> $OUT_DIR/known_mrna.gtf

#determine protein coding potential
FEELnc_codpot.pl \
-i $OUT_DIR/candidate_lncRNA.gtf \
-m shuffle \
-a $REF_GTF \
-g $REF_SEQ \
--verbosity=0 \
--keeptmp \
--outdir $OUT_DIR/feelnc_codpot_out/ 

FEELnc_codpot.pl \
-i $OUT_DIR/candidate_lncRNA.gtf \
-l $NONCODING_FA \
-a $CODING_FA \
-g $REF_SEQ \
--verbosity=0 \
--keeptmp \
--outdir $OUT_DIR/feelnc_codpot_out_integrated/ 
