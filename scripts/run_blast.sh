#!/usr/bin/env bash


if (($# == 0)); then
        echo "Usage:"
        echo "-s = SWISSPROT BLAST output dir"
        echo "-e = E-Value threshold"
        echo "-p = protein sequence"
        echo "-t = Processor numbers"
        exit 2
fi
while getopts s:e:p:s:t: option
  do
    case "${option}"
      in
      e) EVALUE=${OPTARG};;
      p) PROTEIN=${OPTARG};;
      s) SWISSPROT_OUTDIR=${OPTARG};;
      t) THREADS=${OPTARG};;
    esac
done

# Swiss prot
if [ ! -d $SWISSPROT_OUTDIR ]; then
mkdir -p $SWISSPROT_OUTDIR
fi

if [ ! -f $SWISSPROT_OUTDIR/uniprot_sprot.fasta ]; then
echo "Wgetting current uniprot release" && \
wget -q ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
-P $SWISSPROT_OUTDIR 
wget -q ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt \
-P $SWISSPROT_OUTDIR 
gunzip $SWISSPROT_OUTDIR/uniprot_sprot.fasta.gz
fi

if [ ! -f $SWISSPROT_OUTDIR/uniprot_sprot.fasta.phr ]; then
echo "Making blast DB with current uniprot release" && \
makeblastdb -in $SWISSPROT_OUTDIR/uniprot_sprot.fasta -dbtype prot
fi

echo "Blast P" 
blastp \
-query $PROTEIN \
-db $SWISSPROT_OUTDIR/uniprot_sprot.fasta  \
-max_target_seqs 1 \
-outfmt 6 \
-evalue $EVALUE \
-num_threads $THREADS \
> $SWISSPROT_OUTDIR/blastp.outfmt6

echo "----------------------"
echo "BLAST steps complete"
echo "----------------------"
echo "protein BLAST output written to $SWISSPROT_OUTDIR/blastp.outfmt6"
