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
        echo "-g = GFF3 file"
        echo "-s = stringtie GTF"
        echo "-o = path for output"
	echo "-k = lncRNA GTF"
	echo "-c = chromosome lengths"
        exit 2
fi
while getopts l:f:g:s:o:k:c: option
  do
    case "${option}"
      in
      l) LNC_FA=${OPTARG};;
      f) G_FA=${OPTARG};;
      g) GFF3=${OPTARG};;
      s) STRGTF=${OPTARG};;
      o) OUTDIR=${OPTARG};;
      k) LNC_GTF=${OPTARG};;
      c) CHROM_SIZES=${OPTARG};;
    esac
done

#!/usr/bin/env bash

# Ensure output directories exist
mkdir -p $OUTDIR/triplexes
mkdir -p $OUTDIR/nearby
mkdir -p $OUTDIR/bed



#############################################
# Step A: Prepare Gene Annotation BED files  #
#############################################

#####
#i. Reference (pre-lncRNA annotation)
#####

#### Genes

# Generate BED file for reference genes
awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1, $4-1, $5-1, $9, $6, $7}}' "$GFF3" | \
	sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > "$OUTDIR/bed/cbs7435_all_genes.bed"

# Generate BED file for reference genes with CDS
gffread --keep-genes -E -C "$GFF3" | \
	awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1, $4-1, $5-1, $9, $6, $7}}' | \
	sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' >  "$OUTDIR/bed/cbs7435_coding_genes.bed"

#### Exons


# Generate BED file for reference exons
awk -F '[\t;]' 'OFS="\t" {if ($3=="exon") {print $1, $4-1, $5-1, $9, $6, $7}}' "$GFF3" | \
	sed -n 's/\(.*\).*Parent=\(.*\)/\1\2/p' > "$OUTDIR/bed/cbs7435_all_exons.bed"

## Replace transcript names with gene names
awk 'BEGIN { OFS="\t" } NR==FNR { map[$1] = $2; next } { if ($4 in map) $4 = map[$4]; print }' \
     <(gffread -T "$GFF3" | awk -F '[\t";]' 'BEGIN { OFS="\t" } { if ($3=="transcript") print $10, $13 }' -) \
     "$OUTDIR/bed/cbs7435_all_exons.bed" \
     > "$OUTDIR/bed/cbs7435_all_exons_.bed"

mv "$OUTDIR/bed/cbs7435_all_exons_.bed" "$OUTDIR/bed/cbs7435_all_exons.bed"


# Generate BED file for reference genes with CDS
gffread --keep-genes -E -C "$GFF3" -o- | \
	awk -F '[\t;]' 'OFS="\t" {if ($3=="exon") {print $1, $4-1, $5-1, $9, $6, $7}}' | \
	sed -n 's/\(.*\).*Parent=\(.*\)/\1\2/p' >  "$OUTDIR/bed/cbs7435_coding_exons.bed"

## Replace transcript names with gene names
awk 'BEGIN { OFS="\t" } NR==FNR { map[$1] = $2; next } { if ($4 in map) $4 = map[$4]; print }' \
     <(gffread -T "$GFF3" | awk -F '[\t";]' 'BEGIN { OFS="\t" } { if ($3=="transcript") print $10, $13 }' -) \
     "$OUTDIR/bed/cbs7435_coding_exons.bed" \
     > "$OUTDIR/bed/cbs7435_coding_exons_.bed"

mv "$OUTDIR/bed/cbs7435_coding_exons_.bed" "$OUTDIR/bed/cbs7435_coding_exons.bed"

#### Introns

# Compute intron regions by subtracting exons from genes
bedtools subtract -a "$OUTDIR/bed/cbs7435_all_genes.bed" \
	-b "$OUTDIR/bed/cbs7435_all_exons.bed" > \
	"$OUTDIR/bed/cbs7435_all_introns.bed"

# Compute intron regions by subtracting exons from genes
bedtools subtract -a "$OUTDIR/bed/cbs7435_coding_genes.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_exons.bed" > \
	"$OUTDIR/bed/cbs7435_coding_introns.bed"



# Generate promoter regions using TSS positions from the GTF file
awk -F '[\t;]' 'OFS="\t" {
	# For positive strand, TSS is at the gene start
	# For negative strand, TSS is at the gene end
	start = ($6 == "+") ? ($2) : ($3);
	print $1, start, start, $4, $5, $6;
        }' "$OUTDIR/bed/cbs7435_all_genes.bed" > "$OUTDIR/bed/cbs7435_all_genes_tss.bed"

# Extend each TSS to create promoter regions - upstream relative to strand
bedtools slop -i "$OUTDIR/bed/cbs7435_all_genes_tss.bed" \
    -g "$CHROM_SIZES" \
    -l 1500 -r 250 -s > "$OUTDIR/bed/cbs7435_all_genes_promoter.bed"

# Extend each TSS to create promoter regions - upstream relative to strand
bedtools slop -i "$OUTDIR/bed/cbs7435_all_genes.bed" \
    -g "$CHROM_SIZES" \
    -l 1500 -r 250 -s > "$OUTDIR/bed/cbs7435_all_genes_and_promoter.bed"

# Generate promoter regions using TSS positions from the GTF file
awk -F '[\t;]' 'OFS="\t" {
	# For positive strand, TSS is at the gene start (adjusted to 0-based)
	# For negative strand, TSS is at the gene end (adjusted to 0-based)
	start = ($6 == "+") ? ($2) : ($3);
	print $1, start, start, $4, $5, $6;
        }'  "$OUTDIR/bed/cbs7435_coding_genes.bed" > "$OUTDIR/bed/cbs7435_coding_genes_tss.bed"

# Extend each TSS to create promoter regions - upstream relative to strand
bedtools slop -i "$OUTDIR/bed/cbs7435_coding_genes_tss.bed" \
    -g "$CHROM_SIZES" \
    -l 1000 -r 250 -s > "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed"

# Extend each TSS to create promoter regions - upstream relative to strand
bedtools slop -i "$OUTDIR/bed/cbs7435_coding_genes.bed" \
    -g "$CHROM_SIZES" \
    -l 1000 -r 250 -s > "$OUTDIR/bed/cbs7435_coding_genes_and_promoter.bed"

# Generate BED file for lncRNA genes
gffread --keep-genes -E $LNC_GTF | \
	awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1, $4-1, $5, $9, $6, $7}}' | \
	sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > "$OUTDIR/nearby/lncrna_genes.bed"


# Find nearby genes within a 2000 bp window around lncRNA genes
bedtools window -w 2000 \
	-a "$OUTDIR/nearby/lncrna_genes.bed" \
	-b "$OUTDIR/bed/cbs7435_all_genes.bed" \
	> "$OUTDIR/nearby/nearby_all.bed"

awk -F'\t' '$4 != $10 { print $4, $10 }' "$OUTDIR/nearby/nearby_all.bed" \
	> "$OUTDIR/nearby/nearby_genes_all.txt"

# Find nearby genes within a 2000 bp window around lncRNA genes
bedtools window -w 2000 \
	-a "$OUTDIR/nearby/lncrna_genes.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_genes.bed" \
	> "$OUTDIR/nearby/nearby_coding.bed"

awk -F'\t' '$4 != $10 { print $4, $10 }' "$OUTDIR/nearby/nearby_coding.bed" \
	> "$OUTDIR/nearby/nearby_genes_coding.txt"

#############################################
# Step B: Run Triplexator and format output
#############################################

#triplexator -l 20 -e 20 -c 2 -g 20 -mf \
#    -fr on -mrl 7 -mrp 3 -rm 2 -p 32 \
#    -ds "$G_FA" \
#    -ss "$LNC_FA" \
#    -od "$OUTDIR/triplexes" \
#    -o triplexator_result.txt

sort -k1,1 -k2,2n \
	$OUTDIR/triplexes/triplexator_result.txt \
	> $OUTDIR/triplexes/triplexator_result_sort.txt

## Replace transcript names of lncRNA with gene names
awk 'BEGIN { OFS="\t" }
     NR==FNR { map[$1] = $2; next }
     { if ($1 in map) $1 = map[$1]; print }' \
     <(awk -F '[\t";]' 'BEGIN { OFS="\t" } { if ($3=="transcript") print $13, $10 }' "$STRGTF") \
     "$OUTDIR/triplexes/triplexator_result_sort.txt" \
     > "$OUTDIR/triplexes/triplexator_result.txt"

rm $OUTDIR/triplexes/triplexator_result_sort.txt

#############################################
# Step C: TTS Classification and Genome Partitioning
#############################################

# Convert Triplexator TSV output to a BED file for TTS regions
# (Columns: Duplex-ID, TTS start (0-based), TTS end, Sequence-ID, placeholder, strand)
awk 'BEGIN {FS="\t"; OFS="\t"} 
     !/^#/ {print $4, $5, $6, $1, ".", $11}' \
     "$OUTDIR/triplexes/triplexator_result.txt" \
     > "$OUTDIR/triplexes/tts_regions.bed"

# Exclusively assign TTS regions into four classes (priority: promoter > exon > intron > intergenic)
bedtools intersect -a "$OUTDIR/triplexes/tts_regions.bed" \
	 -b "$OUTDIR/bed/cbs7435_all_genes.bed" -wo | sort -k1,1 -k2,2n > "$OUTDIR/triplexes/tts_genes.bed"

echo "Assigning TTS regions exclusively..."

# A. Promoter assignment if TTS overlaps with two promoters then both are reported
bedtools intersect -a "$OUTDIR/triplexes/tts_regions.bed" \
	 -b "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed" -wo \
    > "$OUTDIR/triplexes/tts_coding_promoter_exclusive.bed"

echo "Promoter-exclusive TTS regions are in: $OUTDIR/triplexes/tts_coding_promoter_exclusive.bed"

# B. Exclude promoters from TTS set
bedtools intersect -a "$OUTDIR/triplexes/tts_regions.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed" -v \
	> "$OUTDIR/triplexes/tts_not_promoter.bed"

# C. Exon assignment (from remaining TTS) if TTS overlaps with two exons then both are reported
bedtools intersect -a "$OUTDIR/triplexes/tts_not_promoter.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_exons.bed" -wo \
	> "$OUTDIR/triplexes/tts_exon_exclusive.bed"
echo "Exon-exclusive TTS regions are in: $OUTDIR/triplexes/tts_exon_exclusive.bed"

# D. Exclude exonic entries from the remaining TTS
bedtools intersect -a "$OUTDIR/triplexes/tts_not_promoter.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_exons.bed" -v \
	> "$OUTDIR/triplexes/tts_not_promoter_exon.bed"

# E. Intron assignment (from the remaining TTS) if TTS overlaps with two introns then both are reported
bedtools intersect -a "$OUTDIR/triplexes/tts_not_promoter_exon.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_introns.bed" -wo \
	> "$OUTDIR/triplexes/tts_intron_exclusive.bed"
echo "Intron-exclusive TTS regions are in: $OUTDIR/triplexes/tts_intron_exclusive.bed"

# F. Exclude intronic entries from the remaining TTS
bedtools intersect -a "$OUTDIR/triplexes/tts_not_promoter_exon.bed" \
	-b "$OUTDIR/bed/cbs7435_coding_introns.bed" -v \
	> "$OUTDIR/triplexes/tts_not_promoter_exon_intron.bed"

# G. Intergenic assignment: From the remainder, select those overlapping intergenic regions.
# (First, generate intergenic regions as the complement of assigned regions)
cat <(cut -f1-3 "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed") \
    <(cut -f1-3 "$OUTDIR/bed/cbs7435_coding_exons.bed") \
    <(cut -f1-3 "$OUTDIR/bed/cbs7435_coding_introns.bed") | \
    bedtools sort -i - | bedtools merge -i - > "$OUTDIR/bed/cbs7435_assigned_regions.bed"

bedtools complement -i "$OUTDIR/bed/cbs7435_assigned_regions.bed" \
	-g "$CHROM_SIZES" > "$OUTDIR/bed/cbs7435_intergenic.bed"

bedtools intersect -a "$OUTDIR/triplexes/tts_not_promoter_exon_intron.bed" \
	-b "$OUTDIR/bed/cbs7435_intergenic.bed" -wa -u \
	> "$OUTDIR/triplexes/tts_intergenic_exclusive.bed"



echo "Intergenic-exclusive TTS regions are in: $OUTDIR/triplexes/tts_intergenic_exclusive.bed"

#############################################
# Step D: Generate Promoter-Exclusive Interacting Gene Pairs
#############################################

awk 'NR==FNR { print $4, $10}' "$OUTDIR/triplexes/tts_coding_promoter_exclusive.bed" \
	| sort -u > "$OUTDIR/triplexes/interacting_pairs_coding_promoter.txt"

awk 'NR==FNR { print $4, $10}' "$OUTDIR/triplexes/tts_exon_exclusive.bed" \
	| sort -u > "$OUTDIR/triplexes/interacting_pairs_coding_exon.txt"

awk 'NR==FNR { print $4, $10}' "$OUTDIR/triplexes/tts_intron_exclusive.bed" \
	| sort -u > "$OUTDIR/triplexes/interacting_pairs_coding_intron.txt"

cat "$OUTDIR/triplexes/interacting_pairs_coding_promoter.txt" \
	"$OUTDIR/triplexes/interacting_pairs_coding_exon.txt" \
	"$OUTDIR/triplexes/interacting_pairs_coding_intron.txt" > "$OUTDIR/triplexes/interacting_pairs_coding_all.txt"

cut -f10 "$OUTDIR/triplexes/tts_genes.bed" | sort | uniq -c | sort -k1,1 -nr | awk '{print $1 "\t" $2}' > "$OUTDIR/triplexes/interactions_per_gene.txt"
tail -n +2 "$OUTDIR/triplexes/triplexator_result.txt" | cut -f1 | sort | uniq -c | sort -k1,1 -nr | awk '{print $1 "\t" $2}'> "$OUTDIR/triplexes/interactions_per_lncrna.txt"

#############################################
# Step E: Genome Region Coverage and TTS Count Summary
#############################################

echo "Calculating exclusive genome region coverage..."

# 1. Exclusive Promoter Regions
bedtools sort -i "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed" | \
    bedtools merge -i - > "$OUTDIR/bed/exclusive_promoter.bed"

# 2. Exclusive Exon Regions (remove exon portions overlapping promoters)
bedtools intersect -a "$OUTDIR/bed/cbs7435_coding_exons.bed" \
    -b "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed" -v > "$OUTDIR/bed/exons_no_promoter.bed"
bedtools sort -i "$OUTDIR/bed/exons_no_promoter.bed" | bedtools merge -i - > "$OUTDIR/bed/exclusive_exon.bed"

# 3. Exclusive Intron Regions (remove promoter overlaps from introns, then subtract exclusive exons)
bedtools subtract -a "$OUTDIR/bed/cbs7435_coding_introns.bed" \
    -b "$OUTDIR/bed/cbs7435_coding_genes_promoter.bed" > "$OUTDIR/bed/introns_no_promoter.bed"
bedtools subtract -a "$OUTDIR/bed/introns_no_promoter.bed" -b "$OUTDIR/bed/exclusive_exon.bed" | \
    bedtools merge -i - > "$OUTDIR/bed/exclusive_intron.bed"

# 4. Exclusive Intergenic Regions (complement of assigned regions)
cat <(cut -f1-3 "$OUTDIR/bed/exclusive_promoter.bed") <(cut -f1-3 "$OUTDIR/bed/exclusive_exon.bed") <(cut -f1-3 "$OUTDIR/bed/exclusive_intron.bed") | \
    bedtools sort -i - | bedtools merge -i - > "$OUTDIR/bed/assigned_regions.bed"
bedtools complement -i "$OUTDIR/bed/assigned_regions.bed" -g "$CHROM_SIZES" > "$OUTDIR/bed/exclusive_intergenic.bed"

# Calculate total genome length from the chromosome sizes file.
total_genome_length=$(awk '{sum += $2} END {print sum}' "$CHROM_SIZES")
echo "Total genome length (bp): $total_genome_length"

# Loop over each exclusive region file and calculate its total length and fraction.
for region in "$OUTDIR/bed/exclusive_promoter.bed" "$OUTDIR/bed/exclusive_exon.bed" "$OUTDIR/bed/exclusive_intron.bed" "$OUTDIR/bed/exclusive_intergenic.bed"; do
    region_name=$(basename "$region" .bed)
    region_length=$(awk '{sum += ($3 - $2)} END {print sum}' "$region")
    fraction=$(echo "scale=4; $region_length / $total_genome_length" | bc -l)
    
    echo "-----------------------------------------"
    echo "$region_name: total length (bp) = $region_length"
    echo "Fraction of genome defined as $region_name: $fraction"
done

echo "-----------------------------------------"
echo "Counting TTS regions in each exclusive category..."

rm "$OUTDIR/triplexes/tts_not_promoter_exon_intron.bed" \
	"$OUTDIR/triplexes/tts_not_promoter_exon.bed" \
	"$OUTDIR/triplexes/tts_not_promoter.bed" \
	"$OUTDIR/bed/exclusive_promoter.bed" \
	"$OUTDIR/bed/exclusive_exon.bed" \
	"$OUTDIR/bed/exclusive_intron.bed" \
	"$OUTDIR/bed/exclusive_intergenic.bed"

# Count unique TTSs based on chromosome/start/end
num_tts_promoter=$(cut -f1-3 "$OUTDIR/triplexes/tts_coding_promoter_exclusive.bed" | sort -u | wc -l)
num_tts_exon=$(cut -f1-3 "$OUTDIR/triplexes/tts_exon_exclusive.bed" | sort -u | wc -l)
num_tts_intron=$(cut -f1-3 "$OUTDIR/triplexes/tts_intron_exclusive.bed" | sort -u | wc -l)
num_tts_intergenic=$(cut -f1-3 "$OUTDIR/triplexes/tts_intergenic_exclusive.bed" | sort -u | wc -l)
	

# Count the number of TTS regions in each exclusive category.
#num_tts_promoter=$(wc -l < "$OUTDIR/triplexes/tts_coding_promoter_exclusive.bed")
#num_tts_exon=$(wc -l < "$OUTDIR/triplexes/tts_exon_exclusive.bed")
#num_tts_intron=$(wc -l < "$OUTDIR/triplexes/tts_intron_exclusive.bed")
#num_tts_intergenic=$(wc -l < "$OUTDIR/triplexes/tts_intergenic_exclusive.bed")

num_promoter=$(cut -f10 "$OUTDIR/triplexes/tts_coding_promoter_exclusive.bed" | sort -u | wc -l)
num_exon=$(cut -f10 "$OUTDIR/triplexes/tts_exon_exclusive.bed" | sort -u | wc -l)
num_intron=$(cut -f10 "$OUTDIR/triplexes/tts_intron_exclusive.bed" | sort -u | wc -l)


echo "Number of promoters containing TTS: $num_promoter"
echo "Number of exons containing TTS: $num_exon"
echo "Number of introns containing TTS: $num_intron"
echo "-----------------------------------------"
echo "Number of promoter-exclusive TTS regions: $num_tts_promoter"
echo "Number of exon-exclusive TTS regions: $num_tts_exon"
echo "Number of intron-exclusive TTS regions: $num_tts_intron"
echo "Number of intergenic-exclusive TTS regions: $num_tts_intergenic"




