#!/usr/bin/env bash
#### blast against RFAM,SWISSPROT and MIRBASE with feelnc transcripts and
#### transdecoder sequences
#### inputs are: 1-3)output directories 4) evalue cutoff 5) transcript sequences
#### 6) protein sequences 6) threads
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-s = SWISSPROT BLAST output dir"
        echo "-r = RNA central dir"
        echo "-e = E-Value threshold"
        echo "-n = FEELNc lncRNA nucleotide sequence"
        echo "-p = protein sequence"
        echo "-t = Processor numbers"
        exit 2
fi
while getopts r:s:e:p:n:s:t: option
  do
    case "${option}"
      in
      r) RNACENTRAL_DIR=${OPTARG};;
      e) EVALUE=${OPTARG};;
      p) PROTEIN=${OPTARG};;
      n) NUCLEOTIDE=${OPTARG};;
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

echo "Search RNAcentral.py "

#python scripts/search_RNAcentral_fromfile.py $RNACENTRAL_DIR 1e-5 90
python scripts/search_RNAcentral.py $NUCLEOTIDE $RNACENTRAL_DIR $EVALUE 90

echo "----------------------"
echo "BLAST steps complete"
echo "----------------------"
echo "protein BLAST output written to $SWISSPROT_OUTDIR/blastp.outfmt6"
echo "RNA central output written to $RNACENTRAL_DIR"
