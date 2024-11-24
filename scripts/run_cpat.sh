#!/usr/bin/env bash
#### Determine the protein-coding potential of putative lncRNA transcripts using CPAT
#### This script classifies lncRNA sequences using CPAT 
#### Inputs:
####   -f: FASTA file containing putative lncRNA transcripts
####   -o: Output directory for classification results
####   -c: codinga genes fasta
####   -n: non-coding genes fasta



if (($# == 0)); then
        echo "Usage:"
        echo "-f = lncRNA fasta sequence"
        echo "-o = Output directory"
	echo "-c = coding test file"
	echo "-n = noncoding test file"
        exit 2
fi
while getopts f:o:c:n: option
  do
    case "${option}"
      in
      f) LNCRNA_SEQ=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
      c) CODING_FILE=${OPTARG};;
      n) NON_CODING_FILE=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

#get mouse CPAT data if needed
#if [ ! -f $OUT_DIR/Mouse_logitModel.RData ]; then
#wget https://ayera.dl.sourceforge.net/project/rna-cpat/v1.2.2/prebuilt_model/Mouse_logitModel.RData \
#-P $OUT_DIR
#fi

#if [ ! -f $OUT_DIR/Mouse_Hexamer.tsv ]; then
#wget https://ayera.dl.sourceforge.net/project/rna-cpat/v1.2.2/prebuilt_model/Mouse_Hexamer.tsv \
#-P $OUT_DIR
#fi

if [ ! -f $OUT_DIR/combined_hexamers.tsv ]; then
make_hexamer_tab.py \
-c $CODING_FILE \
-n $NON_CODING_FILE \
> $OUT_DIR/combined_hexamers.tsv

make_logitModel.py \
-x $OUT_DIR/combined_hexamers.tsv \
-c $CODING_FILE \
-n $NON_CODING_FILE \
-o $OUT_DIR/combined \
--log-file $OUT_DIR/make_logitModel_run_info.log
fi

# run CPAT
cpat.py \
--antisense --top-orf=5 \
-g $LNCRNA_SEQ \
-o $OUT_DIR/CPAT.analysis \
-x $OUT_DIR/combined_hexamers.tsv \
-d $OUT_DIR/combined.logit.RData \
--log-file $OUT_DIR/CPAT_run_info.log

# END
