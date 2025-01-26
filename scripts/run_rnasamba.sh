#!/usr/bin/env bash
#### Determine the protein-coding potential of putative lncRNA transcripts using RNAsamba
#### This script classifies lncRNA sequences using RNAsamba with various training models to assess coding potential.
#### Inputs:
####   -f: FASTA file containing putative lncRNA transcripts
####   -o: Output directory for classification results
#### Models used:
####   - Pre-trained, shuffled, and combined (CPPred+Kp) models for ensemble classification



if (($# == 0)); then
        echo "Usage:"
        echo "-f = lncRNA fasta sequence"
        echo "-o = Output directory"
#	echo "-l = lncrna annotation folder"
        exit 2
fi
while getopts f:o: option #l: option
  do
    case "${option}"
      in
      f) LNCRNA_SEQ=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
#      l) LNC_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi


rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_pretrained.tsv \
$LNCRNA_SEQ \
$OUT_DIR/partial_length_weights.hdf5


rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_combined.tsv \
$LNCRNA_SEQ \
$OUT_DIR/combined.hdf5


rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_shuffled.tsv \
$LNCRNA_SEQ \
$OUT_DIR/pichia_shuffle_nc.hdf5



rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_ensemble.tsv \
$LNCRNA_SEQ \
$OUT_DIR/partial_length_weights.hdf5 $OUT_DIR/combined.hdf5

echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_pretrained.tsv | wc -l) "Noncoding predicted by Pre-trained model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_shuffled.tsv | wc -l) "Noncoding predicted by shuffled model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_combined.tsv | wc -l) "Noncoding predicted by CPPred+Kp model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_ensemble.tsv | wc -l) "Noncoding predicted by Ensemble (CPPred+Kp, & Pre-trained) model"


# END
