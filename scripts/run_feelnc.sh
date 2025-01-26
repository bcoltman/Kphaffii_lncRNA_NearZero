#!/usr/bin/env bash
#### Carry out initial lncRNA prediction using feelnc


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
