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

#get S. cerevisiae CPPred data if needed
#if [ ! -f $OUT_DIR/S.cerevisiae_small_coding_RNA.fa ]; then
#echo "wgetting S.cerevisiae small coding RNA" && \
#wget http://rnabinding.com/CPPred/S.cerevisiae-testing/S.cerevisiae_small_coding_RNA.fa -P $OUT_DIR
#fi

#if [ ! -f $OUT_DIR/S.cerevisiae_small_ncrna.fa ]; then
#echo "wgetting S.cerevisiae small non coding RNA" && \
#wget http://rnabinding.com/CPPred/S.cerevisiae-testing/S.cerevisiae_small_ncrna.fa -P $OUT_DIR 
#fi

#if [ ! -f $OUT_DIR/S.cerevisiae_coding_RNA.fa ]; then
#echo "wgetting S.cerevisiae coding RNA" && \
#wget http://rnabinding.com/CPPred/S.cerevisiae-testing/S.cerevisiae_coding_RNA.fa -P $OUT_DIR
#fi

#if [ ! -f $OUT_DIR/S.cerevisiae_ncrna.fa ]; then
#echo "wgetting S.cerevisiae non coding RNA" && \
#wget http://rnabinding.com/CPPred/S.cerevisiae-testing/S.cerevisiae_ncrna.fa -P $OUT_DIR
#fi

# F_TMP=$(ls $FEEL_DIR/tmp/ | cut -d "_" -f1 | uniq | sort -n | tail -1 )

# echo "Using $FEEL_DIR/tmp/${F_TMP}_"



rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_pretrained.tsv \
$LNCRNA_SEQ \
$OUT_DIR/partial_length_weights.hdf5

# run RNAsamba
#rnasamba train \
#-v 2 -s 3 \
#$OUT_DIR/combined.hdf5 \
#$LNC_DIR/combined_coding_rna_test.fa \
#$LNC_DIR/combined_non_coding_rna_test.fa

rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_combined.tsv \
$LNCRNA_SEQ \
$OUT_DIR/combined.hdf5

#rnasamba train \
#-v 2 -s 3 \
#$OUT_DIR/pichia_shuffle_nc.hdf5 \
#$LNC_DIR/FEELnc/tmp/${F_TMP}_candidate_lncRNA.gtf.coding_orf.fa \
#$LNC_DIR/FEELnc/tmp/${F_TMP}_candidate_lncRNA.gtf.noncoding_orf.fa

rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_shuffled.tsv \
$LNCRNA_SEQ \
$OUT_DIR/pichia_shuffle_nc.hdf5

# run RNAsamba
#rnasamba train \
#-v 2 -s 3 \
#$OUT_DIR/s_cerevisiae.hdf5 \
#$OUT_DIR/S.cerevisiae_coding_RNA.fa \
#$OUT_DIR/S.cerevisiae_ncrna.fa 

#rnasamba train \
#-v 2 -s 3 \
#$OUT_DIR/s_cerevisiae_small.hdf5 \
#$OUT_DIR/S.cerevisiae_small_coding_RNA.fa \
#$OUT_DIR/S.cerevisiae_small_ncrna.fa 

rnasamba classify \
-v 1 \
$OUT_DIR/RNAsamba_classification_ensemble.tsv \
$LNCRNA_SEQ \
$OUT_DIR/partial_length_weights.hdf5 $OUT_DIR/combined.hdf5

echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_pretrained.tsv | wc -l) "Noncoding predicted by Pre-trained model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_shuffled.tsv | wc -l) "Noncoding predicted by shuffled model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_combined.tsv | wc -l) "Noncoding predicted by CPPred+Kp model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $OUT_DIR/RNAsamba_classification_ensemble.tsv | wc -l) "Noncoding predicted by Ensemble (CPPred+Kp, & Pre-trained) model"

#$OUT_DIR/pichia_shuffle_nc.hdf5 $OUT_DIR/partial_length_weights.hdf5 $OUT_DIR/combined.hdf5

#rnasamba train \
#-v 2 -s 3 \
#$OUT_DIR/pichia_shuffle_nc.hdf5 \
#$LNC_DIR/FEELnc/tmp/${F_TMP}_candidate_lncRNA.gtf.coding_orf.fa \
#$LNC_DIR/FEELnc/tmp/${F_TMP}_candidate_lncRNA.gtf.noncoding_orf.fa
## $FEEL_DIR/tmp/*_candidate_lncRNA.gtf.*_rna.fa \ 
## $FEEL_DIR/tmp/1483320_candidate_lncRNA.gtf.noncoding_rna.fa 
## if FEElnc is ran a couple of times, 
## then get multiple JOB ids, maybe just output ushuffle in future use

##another training can then create ensemble model 


#rnasamba classify \
#-v 1 \
#$OUT_DIR/RNAsamba_classification.tsv \
#$LNCRNA_SEQ \
#$OUT_DIR/s_cerevisiae.hdf5 $OUT_DIR/s_cerevisiae_small.hdf5 $OUT_DIR/pichia_shuffle_nc.hdf5

# END
