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
awk 'FNR==NR {a[$0]++; next} { for (i in a) { split(i, one); if (one[4] == $4) print one[1], one[2]+$5, one[2]+$6, $4 } }' \
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



# gffread --keep-genes -E $GTF |  awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$9,$6,$7}}' | sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > $OUTDIR/nearby/CBS7435_genes.bed

# ## 50,000 up and down of lncRNAs
# ### This one is on the gene level of the lncRNAs
# gffread --keep-genes -E $LNC_GTF |  awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$9,$6,$7}}' | sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > $OUTDIR/nearby/lncRNA_genes.bed

# # This one ise using the transcripts
# #gffread -E "$LNC_GTF" --bed -o $OUT_DIR/lncRNA.bed
# bedtools window -w 1000 -a $OUTDIR/nearby/lncRNA_genes.bed -b $OUTDIR/nearby/CBS7435_genes.bed > $OUTDIR/nearby/nearby.bed

# awk -F'\t' '$4 != $10 { print $4, $10 }' $OUTDIR/nearby/nearby.bed > $OUTDIR/nearby/nearby_genes.txt

# ## Triplexes

# bedtools slop -i $OUTDIR/nearby/CBS7435_genes.bed -g reference_genome/$CHROM_SIZES -b 1000 > $OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed
# bedtools getfasta -name -fi reference_genome/$G_FA -bed $OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed > $OUTDIR/triplexes/CBS7435_genes_1000UpDown.fa

# ## Predict triplex forming oligos (TFOs) and triplex target sites (TTSs) - merge overlapping
# # Find all triplexes between -ss lncrnas and -ds genes, complying to: 
# #	- at least 20 bps in length "-l 20"
# #	- having at most 20% mismatches and errors "-e 20"
# #	- maximum of 2 consecutive errors in triplex rules "-c 2"
# #	- minimum guanine content in triplex features "-g 20"
# #	- Do not filter for low complexity regions "-fr off"
# #	- filtered for low complexity regions of length >= 7 and period <=3 "-fr on -mrl 7 -mrp 3"
# #	- Merge overlapping features to a cluster "-mf"
# #	- Parallelise on duplexes i.e. genes "-rm 2", using "-p 32" cores


# triplexator -l 20 -e 20 -c 2 -g 20 -mf \
	# -fr on -mrl 7 -mrp 3 \
	# -rm 2 -p 32 \
	# -ds $OUTDIR/triplexes/CBS7435_genes_1000UpDown.fa \
	# -ss $LNC_FA \
	# -od $OUTDIR/triplexes \
	# -o triplexator_predicted_interactions_merged.txt

# awk -F '[\t";]' '{if ($3=="transcript") {print $13,$10}}' $STRGTF > $REFDIR/transcript_genes.txt

# awk 'NR==FNR{ar[$1]=$2;next} ($1 in ar) && $1 != $2 { $1=ar[$1]; print $1, $4 }' $REFDIR/transcript_genes.txt $OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt | sort -u  > $OUTDIR/triplexes/interacting_pairs.txt

# awk 'NR==FNR{filter[$1]; next} ($2 in filter)' $CODING $OUTDIR/triplexes/interacting_pairs.txt > $OUTDIR/triplexes/interacting_pairs_coding.txt

# # Subset genes that interact with lncRNAs and saved as .bed

# awk 'FNR==NR {a[$0]++;next}{ for (i in a) { split(i, one); if (one[4] == $4) print one[1]"\t"one[2]+$5"\t"one[2]+$6"\t"$4} }' $OUTDIR/triplexes/CBS7435_genes_1000UpDown.bed $OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt > $OUTDIR/triplexes/triplexator_interactions.bed
# # Sort bed file
# sort -u -k1,1 -k2,2n  $OUTDIR/triplexes/triplexator_interactions.bed > $OUTDIR/triplexes/triplexator_interactions_sort.bed
# ## Merge the TFO interaction sites, adds column indiciating number of merged sites in each site
# mergeBed -i $OUTDIR/triplexes/triplexator_interactions_sort.bed -c 4 -o distinct,count > $OUTDIR/triplexes/triplexator_interactions_sort_merge.bed
# awk '{OFS="\t"; if($4!~",") {print $0} else {split($4,chr,","); size=length(chr); for(i=1; i<=size;i++) {print $1,$2,$3,chr[i],$5}}}' $OUTDIR/triplexes/triplexator_interactions_sort_merge.bed > $OUTDIR/triplexes/triplexator_interactions_sort_merge_split.bed






# #gffread --keep-genes -E transcriptome_assembly/stringtie.all.transcripts.gtf |  awk -F '[\t;]' 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$9,$6,$7}}' | sed -n 's/\(.*\).*ID=\(.*\)/\1\2/p' > functional_predictions/CBS7435_genes.bed

# #awk -v OFS='\t' {'print $1,$2'} $REFDIR/${G_FA}.fai > $REFDIR/cbs7435.chrom.sizes

# #bedtools slop -i functional_predictions/CBS7435_genes.bed -g reference_genome/cbs7435.chrom.sizes -b 1000 > functional_predictions/CBS7435_genes_1000UpDown.bed
# #bedtools getfasta -name -fi $G_FA -bed functional_predictions/CBS7435_genes_1000UpDown.bed > functional_predictions/CBS7435_genes_1000UpDown.fa


# ## Predict triplex forming oligos (TFOs) and triplex target sites (TTSs) - no merge 

# #triplexator -l 20 -e 20 -c 2 -fr off -rm 1 -p 20 -ds $OUTDIR/CBS7435_genes_1000UpDown.fa -ss $LNC_DIR/lncrna_sequences.fa -od $OUTDIR/triplexes -o triplexator_predicted_interactions.txt




# #triplexator -l 20 -e 20 -c 2 -fr off -rm 1 -p 20 -mf -ds functional_predictions/CBS7435_genes_1000UpDown.fa -ss lncrna_annotation/lncrna_sequences.fa -od functional_predictions/triplexes -o triplexator_predicted_interactions_merged.txt



# #awk -F '[\t";]' '{if ($3=="transcript") {print $10,$13}}' $STRGTF > $REFDIR/transcript_genes.txt
 


# #awk 'NR==FNR{ar[$1]=$2;next}($1 in ar) {$1= ar[$1];print $1,$4}' $REFDIR/transcript_genes.txt $OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt | sort -u  | awk '$1 != $2' > $OUTDIR/triplexes/interacting_pairs.txt


# #awk 'NR==FNR{ar[$1]=$2;next}($1 in ar) {$1= ar[$1];print $1,$4}' reference_genome/transcript_genes.txt functional_predictions/triplexes/triplexator_predicted_interactions_merged.txt | sort -u  | awk '$1 != $2' > functional_predictions/triplexes/interacting_pairs.txt


# ## Extract coordinates for the flanking regions

# #echo "ID Start End UpstreamFlank DownstreamFlank" > $OUTDIR/AllTranscriptsCoords.txt 
# #sed -n '/^>/ s/^>\([^ ]*\).*|\(.*\)-\(.*\)|.*padding:\(.*\)|\(.*\) .*/\1 \2 \3 \4 \5/p' $OUTDIR/1500UpDownAllTranscripts.fa >> $OUTDIR/AllTranscriptsCoords.txt





# #awk '{print $1}' differential_expression/results/R10_vs_SS_0.025.tsv | head -n 501 | tail -n 500 > Top500DEGenes.txt
# #awk '{print $1}' differential_expression/results/R10_vs_SS_0.025.tsv | tail -n 500 > Bottom500DEGenes.txt

# ##################################################################################



# ## Extract sequences of lncRNAs and the flanking regions of all transcripts
# ##gffread -w $LNC_DIR/lncrna_sequences.fa -g $G_FA $LNC_DIR/lncrna.gtf 
# ##gffread -W -w $OUTDIR/1500UpDownAllTranscripts.fa --w-add 1500 -g $G_FA $GTF 



# # when looking at the GFF file, lots of unnecessary infromation is present. Also, the annotation is pretty poor. Look at impelemnting it in gffutils, https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html#gffutils.create_db where the {'gene":ID, "transcript":"ID", }

# # Demosntrates that all mRNA have an ID tag that has .t01 for the trasncript ID
# # cat reference_genome/cbs7435.gff3| awk '{ if ($3 ~ "mRNA")print $0 }' | awk -F'[\t;]' '{ if ($9 ~ ".t01"); else print $0}' | wc -l

# ## Create BED file, need to extract the chromosome ID from the GTF and match it with transcript ID

# #awk 'FNR==NR {a[$0]++;next}{ for (i in a) { split(i, one); if (one[4] == $1) print one[4]"\t"$2-$4+one[5]"\t"$2-$4+one[6]"\t"one[4]":"$2-$4"-"$3+$5} }' $OUTDIR/triplexes/triplexator_predicted_interactions_merged.txt $OUTDIR/AllTranscriptsCoords.txt > $OUTDIR/triplexes/triplexator_exp.txt
# #sed -n '/.*transcript\s/ s/\([^ \t]*\).*transcript_id "\([^ ]\+\)".*/\1 \2 /p' $GTF > $OUTDIR/transcripts_chromosome.txt
# #awk -F "[;\t ]"  'FNR==NR {a[$0]++;next}{ for (i in a) {split(i,one); if ( one[2] == $1) print one[1]"\t"$2"\t"$3"\t"$4}}' $OUTDIR/transcripts_chromosome.txt $OUTDIR/triplexes/triplexator_exp.txt > $OUTDIR/triplexes/triplexator_exp.bed

# #sort -u -k1,1 -k2,2n  <(grep -v TFO functional_predictions/triplexes/triplexator_exp.bed) > functional_predictions/triplexes/triplexator_exp_sort.bed
# #mergeBed -i functional_predictions/triplexes/triplexator_exp_sort.bed -c 4 -o distinct,count > functional_predictions/triplexator_exp_sort_merge.bed

# #awk '{OFS="\t"; if($4!~",") {print $0} else {split($4,chr,","); size=length(chr); for(i=1; i<=size;i++) {print $1,$2,$3,chr[i],$5}}}' functional_predictions/triplexes/triplexator_exp_sort_merge.bed > functional_predictions/triplexes/triplexator_exp_sort_merge_split.bed

