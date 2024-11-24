#!/usr/bin/env bash
## anotation of lncRNAs via RFAM blast and syntheny with human and mouse

if (($# == 0)); then
        echo "Usage:"
        echo "-s lncRNA sequences"
        echo "-o = Output directory"
        exit 2
fi
while getopts g:o:l:a:r: option
  do
    case "${option}"
      in
      a) STR_GTF=${OPTARG};;
      l) LNCRNA_GTF=${OPTARG};;
      g) REF_GTF=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
      r) GENE_ID_NAME=${OPTARG};;
    esac
done

# k) ENS_LNCRNA=${OPTARG};;
# each lncrna is classified using FEELNc
mkdir $OUT_DIR/classification
FEELnc_classifier.pl \
	-i $LNCRNA_GTF \
	-l $OUT_DIR/classification/FEELnc_1st_classification.log \
	-a $REF_GTF \
	> $OUT_DIR/classification/tmp_first_classification.txt


# Create the header
header=$(head -1 "$OUT_DIR/classification/tmp_first_classification.txt")

# Append the new column name to the header
echo -e "${header}\tpartnerRNA_gene_name" > "$OUT_DIR/classification/first_classification.txt"

awk 'BEGIN { FS = OFS = "\t" } 
    FNR==NR { gene_name_map[$1] = $2; next }
    FNR==1 { next }  # Skip the header line of the second file
    { gene_id = $4; gene_name = (gene_id in gene_name_map) ? gene_name_map[gene_id] : "NA"; print $0, gene_name }' \
	$GENE_ID_NAME $OUT_DIR/classification/tmp_first_classification.txt \
	>> $OUT_DIR/classification/first_classification.txt

rm $OUT_DIR/classification/tmp_first_classification.txt

#awk '{ if ($3 == "exon") print $0; }' lncrna_annotation/all_lncrna.gtf | grep -o "transcript_id [^;]\+;" | awk -v FS=" " '{ arr[substr($2, 2, length($2)-3)] += 1; } END { for (t in arr) { if (arr[t] > 1) { print t;} } }'


sed 's/\s/\t/g' $LNCRNA_GTF | \
	awk '$3 == "transcript"'| \
	sed -n 's/.*transcript_id[[:space:]]*"\([^"]*\)".*/\1/p' \
	> $OUT_DIR/monoexonic_filter/first_lncrna_list.txt

## Filter RNA Central searches to only include those lncRNAs decided as "lncRNAs"

awk 'NR==FNR { filter[$1]; next } FNR==1 || $1 in filter' \
	$OUT_DIR/monoexonic_filter/first_lncrna_list.txt \
	$OUT_DIR/RNACentral/RNACentralAll.tsv \
	> $OUT_DIR/RNACentral/RNACentralAllFiltered.tsv



# all monexonic
mkdir $OUT_DIR/monoexonic_filter

# determine which of the putative lncRNAs have 1 exon
sed 's/\s/\t/g' $LNCRNA_GTF | \
	awk '$3 == "exon"'| \
	datamash -s  -g 12 count 3 | \
	awk '$2 < 2' | \
	tr -d \" | \
	tr -d \; | \
	awk '{print $1}' | \
	sort | \
	uniq \
	> $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt

sed 's/\s/\t/g' $LNCRNA_GTF | \
	awk '$3 == "exon"' | \
	datamash -s  -g 12 count 3 | \
	awk '$2 >= 2' | \
	tr -d \" | \
	tr -d \; | \
	awk '{print $1}' | \
	sort | \
	uniq \
	> $OUT_DIR/monoexonic_filter/multiexonic_lncrnas.txt

# one exon lncRNAs are overlapped with the previous RFAM blast
# determine monoexonic with hits in rfam, as well as synteny with human mouse
# rfam blast hits
#awk '{print $1}' $OUT_DIR/RFAM/lncrna.blastn.outfmt6 | sort > \
#$OUT_DIR/monoexonic_filter/rfam_hit_transcripts.txt

#comm -12 $OUT_DIR/monoexonic_filter/rfam_hit_transcripts.txt $OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt >  \
#$OUT_DIR/monoexonic_filter/rfam_monoexonic_lncrnas.txt


# Filtering monoexonic two times, take transcripts that are antisense to protein_coding genes (i.e. 0bp) or transcripts that are more than 500bp away
# Find monoexonic lncRNAs that are > 150 bp away from protein coding gene

awk -F$'\t' '$1 == 1 && $8>=150 {print}' \
	$OUT_DIR/classification/first_classification.txt | \
	awk -F $'\t' '{print $3}' | \
	uniq | \
	sort \
	> $OUT_DIR/monoexonic_filter/feelnc_150bp_away.txt

# Now, if less than < 150 bp, then either, 1) only accept them if they are overlapping & antisense

#awk -F$'\t' '$1 == 1 && $8==0 {print}' lncrna_annotation/classification/first_classification.txt | grep 'antisense' | awk -F $'\t' '{print $3}' | uniq | sort > $OUT_DIR/monoexonic_filter/feelnc_antisense_overlap.txt

# 2) Only accept them if they are antisense

awk -F$'\t' '$1 == 1 {print}' $OUT_DIR/classification/first_classification.txt | \
	grep 'antisense' | \
	awk -F $'\t' '{print $3}' | \
	uniq | \
	sort > $OUT_DIR/monoexonic_filter/feelnc_antisense.txt


# Then, this code works for both cases in 1/2. Just need to uncomment/comment the desired one
cat \
	$OUT_DIR/monoexonic_filter/feelnc_150bp_away.txt \
	$OUT_DIR/monoexonic_filter/feelnc_antisense.txt | sort | uniq \
	> $OUT_DIR/monoexonic_filter/feelnc_antisense_and_150.txt


comm -12 \
	$OUT_DIR/monoexonic_filter/feelnc_antisense_and_150.txt \
	$OUT_DIR/monoexonic_filter/monoexonic_transcripts.txt > \
	$OUT_DIR/monoexonic_filter/antisense_or_150_monoexonic_lncrnas.txt

cat \
	$OUT_DIR/monoexonic_filter/multiexonic_lncrnas.txt \
	$OUT_DIR/monoexonic_filter/antisense_or_150_monoexonic_lncrnas.txt | \
	sort| uniq > $OUT_DIR/monoexonic_filter/conservative_lncrna_list.txt

#$OUT_DIR/monoexonic_filter/rfam_monoexonic_lncrnas.txt \

grep -wFf \
	$OUT_DIR/monoexonic_filter/conservative_lncrna_list.txt \
	$STR_GTF \
	> $OUT_DIR/conservative_lncrna.gtf


## rerun classifer on final lncrna list
FEELnc_classifier.pl \
	-l $OUT_DIR/classification/FEELnc_conservative_classification.log \
	-i $OUT_DIR/conservative_lncrna.gtf \
	-a $REF_GTF \
	> $OUT_DIR/classification/tmp_conservative_classification.txt


# Create the header
header=$(head -1 "$OUT_DIR/classification/tmp_conservative_classification.txt")

# Append the new column name to the header
echo -e "${header}\tpartnerRNA_gene_name" > "$OUT_DIR/classification/conservative_classification.txt"

awk 'BEGIN { FS = OFS = "\t" } 
    FNR==NR { gene_name_map[$1] = $2; next } 
    FNR==1 { next }  # Skip the header line of the second file
    { gene_id = $4; gene_name = (gene_id in gene_name_map) ? gene_name_map[gene_id] : "NA"; print $0, gene_name }' \
	$GENE_ID_NAME $OUT_DIR/classification/tmp_conservative_classification.txt \
	>> $OUT_DIR/classification/conservative_classification.txt

rm $OUT_DIR/classification/tmp_conservative_classification.txt
