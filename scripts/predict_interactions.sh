#!/usr/bin/env bash
#### Triplex Forming Oligonucleotide (TFO) Prediction and Interaction Analysis
#### This script predicts TFO interactions between lncRNAs and target genes.
#### Inputs:
####   -l: Path to lncRNA FASTA file
####   -f: Path to reference genome FASTA file
####   -g: Combined GTF file with transcript annotations
####   -s: StringTie GTF file with transcript annotations
####   -r: Reference genome directory
####   -o: Output directory for results
####   -k: GTF file with lncRNA annotations
####   -c: File with chromosome lengths
####   -p: GTF file with coding gene annotations



if (($# == 0)); then
        echo "Usage:"
        echo "-l = lncRNA fasta"
        echo "-f = fasta file genome"
        echo "-g = Combined GTF file"
        echo "-s = stringtie GTF"
        echo "-r = Reference genome dir"
        echo "-o = path for output"
		echo "-k = lncRNA GTF"
		echo "-c = chromosome lengths"
		echo "-p = coding genes GTF"
        exit 2
fi
while getopts l:f:g:s:r:o:k:c:p: option
  do
    case "${option}"
      in
      l) LNC_FA=${OPTARG};;
      f) G_FA=${OPTARG};;
      g) GTF=${OPTARG};;
      s) STRGTF=${OPTARG};;
      r) REFDIR=${OPTARG};;
      o) OUTDIR=${OPTARG};;
      k) LNC_GTF=${OPTARG};;
      c) CHROM_SIZES=${OPTARG};;
      p) CODING=${OPTARG};;
    esac
done

#!/usr/bin/env bash

# Ensure output directories exist
mkdir -p $OUTDIR/triplexes
mkdir -p $OUTDIR/nearby

# Generate BED file for reference genes
gffread --keep-genes -E $GTF | \
	awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1, $4-1, $5, $9, $6, $7}}' | \
	sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > $OUTDIR/nearby/CBS7435_genes.bed

# Generate BED file for lncRNA genes
gffread --keep-genes -E $LNC_GTF | \
	awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1, $4-1, $5, $9, $6, $7}}' | \
	sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > $OUTDIR/nearby/lncRNA_genes.bed

# Find nearby genes within 1000 bp window around lncRNA genes
bedtools window -w 1000 \
	-a $OUTDIR/nearby/lncRNA_genes.bed \
	-b $OUTDIR/nearby/CBS7435_genes.bed \
	> $OUTDIR/nearby/nearby.bed
awk -F'\t' '$4 != $10 { print $4, $10 }' $OUTDIR/nearby/nearby.bed > $OUTDIR/nearby/nearby_genes.txt

# Expand CBS7435 genes by 1000 bp on both sides and extract sequences
bedtools slop -i $OUTDIR/nearby/CBS7435_genes.bed \
	-g reference_genome/$CHROM_SIZES \
	-b 1000 \
	> $OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed
bedtools getfasta -name \
	-fi reference_genome/$G_FA \
	-bed $OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed \
	> $OUTDIR/triplexes/CBS7435_genes_1000UpDown.fa

# Predict Triplex Forming Oligos (TFOs) with Triplexator
triplexator -l 20 -e 20 -c 2 -g 20 -mf \
	-fr on -mrl 7 -mrp 3 -rm 2 -p 32 \
	-ds $OUTDIR/triplexes/CBS7435_genes_1000UpDown.fa \
	-ss $LNC_FA \
	-od $OUTDIR/triplexes \
	-o triplexator_predicted_interactions_merged.txt

# Map transcripts to gene IDs
awk -F '[\t";]' '{if ($3=="transcript") {print $13, $10}}' $STRGTF > $REFDIR/transcript_genes.txt

# Generate list of interacting gene pairs
awk 'NR==FNR {ar[$1]=$2; next} ($1 in ar) && $1 != $2 { $1=ar[$1]; print $1, $4 }' \
	$REFDIR/transcript_genes.txt $OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt | \
	sort -u > $OUTDIR/triplexes/interacting_pairs.txt

# Filter for coding gene interactions
awk 'NR==FNR {filter[$1]; next} ($2 in filter)' \
	$CODING \
	$OUTDIR/triplexes/interacting_pairs.txt \
	> $OUTDIR/triplexes/interacting_pairs_coding.txt

# Subset and format interacting genes for .bed file
awk 'BEGIN { OFS="\t" } FNR==NR {a[$0]++; next} { for (i in a) { split(i, one); if (one[4] == $4) print one[1], one[2]+$5, one[2]+$6, $4 } }' \
	$OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed \
	$OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt \
	> $OUTDIR/triplexes/triplexator_interactions.bed

# Sort and merge BED file to identify TFO interaction sites
sort -u -k1,1 -k2,2n \
	$OUTDIR/triplexes/triplexator_interactions.bed \
	> $OUTDIR/triplexes/triplexator_interactions_sort.bed

mergeBed -i $OUTDIR/triplexes/triplexator_interactions_sort.bed \
	-c 4 \
	-o distinct,count \
	> $OUTDIR/triplexes/triplexator_interactions_sort_merge.bed

# Split merged entries into separate lines if there are multiple interactions
awk '{
  OFS="\t"; if ($4 !~ ",") {print $0} 
  else {split($4, chr, ","); for (i = 1; i <= length(chr); i++) {print $1, $2, $3, chr[i], $5}}
}' $OUTDIR/triplexes/triplexator_interactions_sort_merge.bed \
	> $OUTDIR/triplexes/triplexator_interactions_sort_merge_split.bed



