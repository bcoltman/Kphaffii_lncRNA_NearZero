#!/usr/bin/env bash
#### HMMER search of TransDecoder protein sequences against Pfam database
#### This adapted script performs a HMMER search using specified protein sequences and e-value cutoff.
#### Inputs:
####   1) Number of threads
####   2) E-value cutoff
####   3) Input protein sequences
####   4) Output directory for results


if (($# == 0)); then
        echo "Usage:"
        echo "-t = threads"
        echo "-e = E-Value threshold"
        echo "-p = protein protein sequences"
        echo "-o = Output directory"
        exit 2
fi
while getopts t:e:p:o: option
  do
    case "${option}"
      in
      t) THREADS=${OPTARG};;
      e) EVALUE=${OPTARG};;
      p) PROTEIN=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

if [ ! -d $OUT_DIR ]; then
mkdir -p $OUT_DIR
fi

if [  -f $OUT_DIR/Pfam-A.hmm.h3i ]; then
rm $OUT_DIR/Pfam-A.hmm.h3i
fi


if [ ! -f $OUT_DIR/Pfam-A.hmm ]; then
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
-P $OUT_DIR
gunzip $OUT_DIR/Pfam-A.hmm.gz
fi


hmmpress $OUT_DIR/Pfam-A.hmm
hmmscan \
--cpu $THREADS \
-E $EVALUE \
--domtblout $OUT_DIR/pfam.domtblout \
$OUT_DIR/Pfam-A.hmm \
$PROTEIN

# make a list of ids for filtering from PFAM table
awk '{print $4}' $OUT_DIR/pfam.domtblout | \
grep -E 'PP7435|MSTRG' | uniq | sed 's/.p.*//g' > \
$OUT_DIR/pfam_domain_transcripts
