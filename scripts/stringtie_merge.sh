#!/usr/bin/env bash
#### Merge individual transcript assemblies into a single reference
#### Inputs:
####   1) Parent directory containing assembled transcripts
####   2) Reference GTF file
####   3) Merged GTF output file
####   4) Reference genome directory


if (($# == 0)); then
        echo "Usage:"
        echo "-t = assembled transcript directory"
        echo "-g = path to reference annotation"
	echo "-m = merged gtf"
        echo "-r = reference genome directory"
        exit 2
fi
while getopts t:g:m:r: option
  do
    case "${option}"
      in
      t) TRANSCRIPT_DIR=${OPTARG};;
      g) GFF=${OPTARG};;
      m) MGTF=${OPTARG};;
      r) REF_DIR=${OPTARG};;
    esac
done

readlink -f $TRANSCRIPT_DIR/individual_gtfs/*.gtf > $TRANSCRIPT_DIR/mergelist.txt

# introduced the -g which decides the gap_len between transcripts to merge
#stringtie \
#--merge $TRANSCRIPT_DIR/mergelist.txt \
#-o $MGTF \
#-G $GFF \
#-f 0.1  \
#-c 10 

stringtie \
--merge $TRANSCRIPT_DIR/mergelist.txt \
-o $MGTF \
-G $GFF \
-F 0.1  \
-c 10 
# create a file liniking stringtie ID to Pichia Genome DB ID DONE

grep -wFf $REF_DIR/protein.coding.genes.txt $MGTF | \
grep -v exon | awk '{print $10, $NF}' | uniq | tr -d \" | tr -d \; > \
$TRANSCRIPT_DIR/stringtie_PGDB_gene_mapping.txt


# make a version of the stringtie assembly for differential exression containing
# transcripts

perl scripts/mstrg_prep.pl $MGTF > $TRANSCRIPT_DIR/stringtie_original_appended.gtf

grep 'MSTRG.*|PP7435.*|PP7435.*' $TRANSCRIPT_DIR/stringtie_original_appended.gtf | \
grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt

grep -v -F -f $TRANSCRIPT_DIR/removed_overlapped_mstrg_transcripts.txt $TRANSCRIPT_DIR/stringtie_original_appended.gtf > \
$TRANSCRIPT_DIR/stringtie_original.appended.fp.filtered.gtf

# remove transcripts without strand information
awk '$7 != "." {print}' $TRANSCRIPT_DIR/stringtie_original.appended.fp.filtered.gtf > \
$TRANSCRIPT_DIR/stringtie.all.transcripts.gtf

gffcompare \
-o $TRANSCRIPT_DIR/gffc_all \
-r $GFF  $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf

# create a stringtie GTF of protein coding genes only

grep -wFf $REF_DIR/protein.coding.genes.txt $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf > \
$TRANSCRIPT_DIR/stringtie.protein.coding.gtf


#######

#Some transcripts being classed as non-coding, despite their gene_id being in the .txt file. Why do they get through?


########

# make to GTF file to search for lncRNAs
# filter transcripts with 1bp exon overlap with annotated protein coding genes
#create bedtools overlap, only for lncRNA
awk '{if($3=="exon"){print $10}}' $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf | sed 's/"//g;s/;//g' | sort | uniq > $TRANSCRIPT_DIR/stringtie.gene.txt

awk '{if($3=="exon"){print}}' $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf > $TRANSCRIPT_DIR/stringtie.exon.gtf

## If use this one, then remove any transcripts overlapping with exons of coding genes
gffread $GFF -C -T | awk '{if($3=="exon"){print}}' > \
$TRANSCRIPT_DIR/reference.exon.gtf

# If using this on, remove any transcripts overlapping with ANY annotated exons
gffread $GFF -T | awk '{if($3=="exon"){print}}' > \
$TRANSCRIPT_DIR/reference.exon.gtf

bedtools intersect -s -u -a $TRANSCRIPT_DIR/stringtie.exon.gtf \
                         -b $TRANSCRIPT_DIR/reference.exon.gtf > \
                          $TRANSCRIPT_DIR/overlap_exon.gtf

awk '{print $10}' $TRANSCRIPT_DIR/overlap_exon.gtf | sed 's/"//g;s/;//g' | sort | uniq  > \
$TRANSCRIPT_DIR/overlap.exon.gene.id.txt

comm -23 $TRANSCRIPT_DIR/stringtie.gene.txt $TRANSCRIPT_DIR/overlap.exon.gene.id.txt | \
grep -wFf - $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf > $TRANSCRIPT_DIR/stringtie_removed_overlapping.gtf

# remove vhh_plasmid 
grep -v ^vhh_plasmid $TRANSCRIPT_DIR/stringtie_removed_overlapping.gtf  > $TRANSCRIPT_DIR/non_protein_coding_stringtie.gtf

# compare to the CBS7435 reference sequence
gffcompare \
-o $TRANSCRIPT_DIR/gffc_nc \
-r $GFF  $TRANSCRIPT_DIR/non_protein_coding_stringtie.gtf



awk '
BEGIN { FS="\t"; OFS="\t" }
$3 == "transcript" {
    ref_gene_id = ""; gene_name = "";
    split($9, fields, ";");
    for (i in fields) {
        if (match(fields[i], /ref_gene_id "([^"]+)"/, arr)) {
            ref_gene_id = arr[1];
        }
        if (match(fields[i], /gene_name "([^"]+)"/, arr)) {
            gene_name = arr[1];
        }
    }
    if (ref_gene_id != "" && gene_name != "") {
        print ref_gene_id, gene_name
    }
} awk ' $TRANSCRIPT_DIR/stringtie.all.transcripts.gtf | sort -u > $TRANSCRIPT_DIR/ref_gene_id_name.tsv

# END
