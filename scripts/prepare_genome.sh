#!/usr/bin/env bash

while getopts o: option
  do
    case "${option}"
      in
      o) GENOME_DIR=${OPTARG};;
    esac
done

mkdir -p $GENOME_DIR/retentostat/

# get reference genome - only using sequence
NCBI_LOC=ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/223/565/GCA_000223565.1_PicPas_Mar2011/GCA_000223565.1_PicPas_Mar2011_assembly_structure/
wget --accept fna.gz --no-parent --quiet --recursive --no-directories ftp://$NCBI_LOC -P reference_genome


gunzip $GENOME_DIR/*.gz

name=cbs7435

# Change scaffold ID in the fasta header to be compatible with Pichia genome GFF
# Rename file to the same scaffold ID 
# Change files from fna format of NCBI to fa format. Just an extension change


rm -r $GENOME_DIR/fasta; mkdir $GENOME_DIR/fasta

for f in $GENOME_DIR/*fna; do
	varname=${f#*/}
	varname=${varname%.fna}
	varname=$(echo $varname | tr [':upper:'] [':lower:'])
	varname=${varname/chrmt/mt}

	sed --in-place "/^>/ s/.*/>${name}_$varname/" $f
	mv -- "$f" "$GENOME_DIR/fasta/${name}_${varname}.fa"

done

# concatenate seperate fasta files

cat $GENOME_DIR/fasta/*.fa > $GENOME_DIR/${name}_combined.fa

samtools faidx $GENOME_DIR/${name}_combined.fa


## Move GFF3 file from original data directory

#cp data/cbs7435.gff3 $GENOME_DIR/${name}.gff3
cp data/20240417_GFFCleanup/cbs7435_20240417_complete_temp_gffread.gff3 $GENOME_DIR/${name}.gff3



cat $GENOME_DIR/${name}_combined.fa data/additional_sequences/vhh_plasmid.fa > $GENOME_DIR/retentostat/${name}_retentostat_combined.fa


tail -n+2 data/additional_sequences/vhh_plasmid.gff3 | cat $GENOME_DIR/${name}.gff3 - > $GENOME_DIR/retentostat/${name}_retentostat.gff3

#list and gtf of ensembl protein coding genes

grep $'\t''CDS' $GENOME_DIR/retentostat/${name}_retentostat.gff3 > $GENOME_DIR/retentostat/${name}_protein.gff3
grep $'\t''CDS' $GENOME_DIR/retentostat/${name}_retentostat.gff3 | awk {'print $9'} | awk -F '[=;.]' '{print $2}' | uniq > $GENOME_DIR/retentostat/protein.coding.genes.txt


# List of other ncRNAs. To avoid hits for genes with ncRNAs mentioned in description, perform on first columns


awk '$3 ~ /^(ncRNA|rRNA|snoRNA|tRNA|spliceosomalRNA|scpRNA|miRNA|Mt_rRNA|Mt_tRNA|processed_pseudogene|pseudogene|ribozyme|rRNA|scaRNA|snoRNA|snRNA|sRNA)$/' $GENOME_DIR/retentostat/${name}_retentostat.gff3 | awk {'print $9'} | awk -F '[=;.]' '{print $2}' | uniq > $GENOME_DIR/other.noncoding.genes.txt

awk '$3 ~ /^(ncRNA|snoRNA|tRNA|spliceosomalRNA|scpRNA|miRNA|Mt_rRNA|Mt_tRNA|processed_pseudogene|pseudogene|ribozyme|scaRNA|snoRNA|snRNA|sRNA)$/' $GENOME_DIR/retentostat/${name}_retentostat.gff3 | awk {'print $9'} | awk -F '[=;.]' '{print $2}' | uniq > $GENOME_DIR/other.noncoding.no_rrna_genes.txt
 


# No lincRNA in Pichia CBS7435 annotation

gffread -w $GENOME_DIR/retentostat/${name}_retentostat_transcripts.fa \
	-g $GENOME_DIR/retentostat/${name}_retentostat_combined.fa \
	$GENOME_DIR/retentostat/${name}_retentostat.gff3


# END
