#!/usr/bin/env bash
#### lncRNA and coding RNA classification pipeline
#### This script tests the different classifiers on known coding and non-coding RNA sequences, 
#### evaluates various prediction models, and reports misclassification rates for each method.
#### Inputs:
####   -g: File with known coding RNA sequences in FASTA format
####   -n: File with known non-coding RNA sequences in FASTA format
####   -l: Directory for lncRNA annotations and output files


if (($# == 0)); then
        echo "Usage:"
        echo "-g = Known coding fa"
	echo "-n = Known non-coding fa"
	echo "-l = lncRNA annotation directory"
        exit 2
fi
while getopts g:n:l: option
  do
    case "${option}"
      in
      g) CODING_FA_TEST=${OPTARG};;
      n) NONCODING_FA_TEST=${OPTARG};;
      l) LNC_DIR=${OPTARG};;
    esac
done


mkdir -p $LNC_DIR/test_classifiers/coding \
	$LNC_DIR/test_classifiers/non_coding \
	$LNC_DIR/CPAT $LNC_DIR/RNAsamba \
	$LNC_DIR/CPC2


####################
### CPAT


make_hexamer_tab.py \
-c $LNC_DIR/data_fasta/Integrated_coding_RNA_test.fa \
-n $LNC_DIR/data_fasta/Integrated_ncRNA_test.fa \
> $LNC_DIR/CPAT/combined_hexamers.tsv

make_logitModel.py \
-x $LNC_DIR/CPAT/combined_hexamers.tsv \
-c $LNC_DIR/data_fasta/Integrated_coding_RNA_test.fa \
-n $LNC_DIR/data_fasta/Integrated_ncRNA_test.fa \
-o $LNC_DIR/CPAT/combined \
--log-file $LNC_DIR/CPAT/make_logitModel_run_info.log


Rscript ./scripts/10Fold_CrossValidation.r "${LNC_DIR}/CPAT/combined.feature.xls" "${LNC_DIR}/CPAT/10Fold_CV"

# run CPAT
cpat.py \
--antisense --top-orf=5 \
-g $CODING_FA_TEST \
-o $LNC_DIR/test_classifiers/coding/CPAT.analysis \
-x $LNC_DIR/CPAT/combined_hexamers.tsv \
-d $LNC_DIR/CPAT/combined.logit.RData \
--log-file $LNC_DIR/test_classifiers/coding/CPAT_run_info.log

# run CPAT
cpat.py \
--antisense --top-orf=5 \
-g $NONCODING_FA_TEST \
-o $LNC_DIR/test_classifiers/non_coding/CPAT.analysis \
-x $LNC_DIR/CPAT/combined_hexamers.tsv \
-d $LNC_DIR/CPAT/combined.logit.RData \
--log-file $LNC_DIR/test_classifiers/non_coding/CPAT_run_info.log

################
## CPC2
bin/CPC2_standalone-1.0.1/bin/CPC2.py \
  -i $NONCODING_FA_TEST -o $LNC_DIR/test_classifiers/non_coding/CPC2.analysis

bin/CPC2_standalone-1.0.1/bin/CPC2.py \
  -i $CODING_FA_TEST -o $LNC_DIR/test_classifiers/coding/CPC2.analysis

#################
## RNAsamba


curl https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/partial_length_weights.hdf5 -o $LNC_DIR/RNAsamba/partial_length_weights.hdf5

## Classify with pre-trained model

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_coding_pretrained.tsv \
$CODING_FA_TEST \
$LNC_DIR/RNAsamba/partial_length_weights.hdf5

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_non_coding_pretrained.tsv \
$NONCODING_FA_TEST \
$LNC_DIR/RNAsamba/partial_length_weights.hdf5

## Train new model with integrated data from CPPpred and then classify
cp data/integrated.hdf5 $LNC_DIR/RNAsamba/integrated.hdf5

#rnasamba train \
#-v 2 -s 3 \
#$LNC_DIR/RNAsamba/integrated.hdf5 \
#$LNC_DIR/data_fasta/Integrated_coding_RNA_test.fa \
#$LNC_DIR/data_fasta/Integrated_ncRNA_test.fa

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_coding_integrated.tsv \
$CODING_FA_TEST \
$LNC_DIR/RNAsamba/integrated.hdf5

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_non_coding_integrated.tsv \
$NONCODING_FA_TEST \
$LNC_DIR/RNAsamba/integrated.hdf5


## Train new model with shuffled FEELnc data and then classify

shuffle_coding=$(find "$LNC_DIR/FEELnc/feelnc_codpot_out/tmp" -type f -name '*_candidate_lncRNA.gtf.coding_orf.fa' | sort -n | tail -n 1)
shuffle_noncoding=$(find "$LNC_DIR/FEELnc/feelnc_codpot_out/tmp" -type f -name '*_candidate_lncRNA.gtf.noncoding_orf.fa' | sort -n | tail -n 1)

echo "Using:"
echo $shuffle_coding
echo $shuffle_noncoding

cp data/pichia_shuffle_nc.hdf5 $LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5
#rnasamba train \
#-v 2 -s 3 \
#$LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5 \
#$shuffle_coding \
#$shuffle_noncoding


rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_coding_shuffled.tsv \
$CODING_FA_TEST \
$LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_non_coding_shuffled.tsv \
$NONCODING_FA_TEST \
$LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_non_coding_ensemble.tsv \
$NONCODING_FA_TEST \
$LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5 $LNC_DIR/RNAsamba/partial_length_weights.hdf5 $LNC_DIR/RNAsamba/integrated.hdf5

rnasamba classify \
-v 1 \
$LNC_DIR/test_classifiers/RNAsamba_coding_ensemble.tsv \
$CODING_FA_TEST \
$LNC_DIR/RNAsamba/pichia_shuffle_nc.hdf5 $LNC_DIR/RNAsamba/partial_length_weights.hdf5 $LNC_DIR/RNAsamba/integrated.hdf5



#### Train merged integrated and KP data

cp data/combined.hdf5 $LNC_DIR/RNAsamba/combined.hdf5
#rnasamba train \
#-v 2 -s 3 \
#$LNC_DIR/RNAsamba/combined.hdf5 \
#$LNC_DIR/data_fasta/combined_coding_rna_test.fa \
#$LNC_DIR/data_fasta/combined_non_coding_rna_test.fa


### Count missed

n_coding=$(sed -n '/^>/p' $CODING_FA_TEST | wc -l)
n_noncoding=$(sed -n '/^>/p' $NONCODING_FA_TEST | wc -l)
CPAT_coding=$(awk -F'\t' '$11 < 0.4 { print $11 }'  $LNC_DIR/test_classifiers/coding/CPAT.analysis.ORF_prob.best.tsv | wc -l)
CPAT_noncoding=$(awk -F'\t' '$11 > 0.4 { print $11 }'  $LNC_DIR/test_classifiers/non_coding/CPAT.analysis.ORF_prob.best.tsv | wc -l)
CPC2_coding=$(awk -F'\t' '$8 == "noncoding" {print $8}' $LNC_DIR/test_classifiers/coding/CPC2.analysis.txt | wc -l)
CPC2_noncoding=$(awk -F'\t' '$8 == "coding" {print $8}' $LNC_DIR/test_classifiers/non_coding/CPC2.analysis.txt | wc -l)


echo "RNAsamba non-coding"
echo  $(awk '{ if ($3 == "coding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_non_coding_pretrained.tsv | wc -l) "Out of" $n_noncoding "Non-coding incorrectly predicted by pre-trained model"
echo  $(awk '{ if ($3 == "coding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_non_coding_integrated.tsv | wc -l) "Out of" $n_noncoding "Non-coding incorrectly predicted by integrated (CPP) model"
echo  $(awk '{ if ($3 == "coding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_non_coding_shuffled.tsv | wc -l) "Out of" $n_noncoding "Non-coding incorrectly predicted by shuffled model"
echo  $(awk '{ if ($3 == "coding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_non_coding_ensemble.tsv | wc -l) "Out of" $n_noncoding "Non-coding incorrectly predicted by ensemble model"

echo "RNAsamba coding"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_coding_pretrained.tsv | wc -l) "Out of" $n_coding "Coding incorrectly predicted by pre-trained model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_coding_integrated.tsv | wc -l) "Out of" $n_coding "Coding incorrectly predicted by integrated (CPP) model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_coding_shuffled.tsv | wc -l) "Out of" $n_coding "Coding incorrectly predicted by shuffled model"
echo  $(awk '{ if ($4 == "noncoding") print $0}' $LNC_DIR/test_classifiers/RNAsamba_coding_ensemble.tsv | wc -l) "Out of" $n_coding "Coding incorrectly predicted by ensemble model"

echo "CPAT"

echo  $CPAT_coding " Out of " $n_coding "Coding incorrectly predicted by CPAT"
echo  $CPAT_noncoding " Out of " $n_noncoding "Non-coding incorrectly predicted by CPAT"

echo "CPC2"
echo  $CPC2_coding " Out of " $n_coding "Coding incorrectly predicted by CPC2"
echo  $CPC2_noncoding " Out of " $n_noncoding "Non-coding incorrectly predicted by CPC2"

# END
