### Work notebook for project:

#### Downloading SRA files - need to include a file from NCBI with the SRA data attached to it to then loop through and use SRR number in below format. Can also use it for study design too

```bash
mkdir -p data/retentostat/sra
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do
	prefetch $sample  -O data/retentostat/sra
done
```

```bash
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do
	fasterq-dump data/retentostat/sra/$sample --split-files --outdir data/retentostat/raw
	pigz data/retentostat/raw/${sample}_1.fastq 
	pigz data/retentostat/raw/${sample}_2.fastq 
done
```

# FASTQC and MULTIQC of Raw reads
```bash
mkdir -p QC/retentostat/multiQC
fastqc -q -t 32 -o QC/retentostat/raw data/retentostat/raw/*.gz
multiqc QC/retentostat/raw -n raw -o QC/retentostat/multiQC
```

## Cutadapt and FastQC
```bash
mkdir -p data/retentostat/cutadapt
mkdir -p QC/retentostat/post_cut

awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do
	cutadapt  \
		-j 32 \
		-m 1 \
		-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
		-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
		-o data/retentostat/cutadapt/${sample}_1.fastq.gz \
		-p data/retentostat/cutadapt/${sample}_2.fastq.gz \
		   data/retentostat/raw/${sample}_1.fastq.gz \
		   data/retentostat/raw/${sample}_1.fastq.gz \
		   1>data/retentostat/cutadapt/${sample}_report.txt 
done

fastqc -q -t 32 -o QC/retentostat/post_cut data/retentostat/cutadapt/*.gz
multiqc QC/retentostat/post_cut -n post_cut -o QC/retentostat/multiQC
```

## Trimmomatic and FastQC

```bash
mkdir -p data/retentostat/preprocessed/paired \
	data/retentostat/preprocessed/unpaired
while read -r sample; do 
	java -jar ../bin/trimmomatic-0.36.jar PE \
		 -threads 1 \
		 data/retentostat/cutadapt/${sample}_1.fastq.gz \
		 data/retentostat/cutadapt/${sample}_2.fastq.gz \
		 data/retentostat/preprocessed/paired/${sample}_1.fastq.gz \
		 data/retentostat/preprocessed/unpaired/${sample}_1.fastq.gz \
		 data/retentostat/preprocessed/paired/${sample}_2.fastq.gz \
		 data/retentostat/preprocessed/unpaired/${sample}_2.fastq.gz \
		 SLIDINGWINDOW:4:20 \
		 MINLEN:36 \
		 -trimlog data/retentostat/preprocessed/${sample}.trimmomatic.log &
done < <(awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv)
echo "Waiting"
wait
echo "Trim quality finished" 
```


# Read Mapping

## Download and prepare reference genome.
```bash
./scripts/prepare_genome.sh -o reference_genome
```

## Make STAR Index

```bash
mkdir -p reference_genome/retentostat/star_index
STAR \
     --runThreadN 32 \
     --runMode genomeGenerate \
     --sjdbOverhang 149 \
     --genomeChrBinNbits 16 \
     --genomeSAindexNbases 10 \
     --genomeDir reference_genome/retentostat/star_index \
     --genomeFastaFiles reference_genome/retentostat/cbs7435_retentostat_combined.fa \
     --sjdbGTFfile reference_genome/retentostat/cbs7435_retentostat.gff3 \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFtagExonParentGene locus_tag
```


## Mapt to CBS7435 Genome
```bash
mkdir data/retentostat/mapped
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do 
	STAR \
		--runThreadN 16 \
		--readFilesIn data/retentostat/preprocessed/paired/${sample}_1.fastq.gz \
			data/retentostat/preprocessed/paired/${sample}_2.fastq.gz \
		--genomeDir reference_genome/retentostat/star_index \
		--readFilesCommand pigz -dc \
		--outFileNamePrefix data/retentostat/mapped/${sample} \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--quantMode GeneCounts \
		--outReadsUnmapped Fastx \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignIntronMin 20 \
		--outFilterType BySJout \
		--sjdbScore 1 \
		--limitBAMsortRAM 10000000000
done
```


## Genome guided assembly
### string tie assembly
```bash
mkdir -p transcriptome_assembly/individual_gtfs
while read -r sample; do
	stringtie \
		-p 1 \
		-G reference_genome/retentostat/cbs7435_retentostat.gff3 \
		--rf \
		-j 5 \
		-f 0.03 \
		-o transcriptome_assembly/individual_gtfs/${sample}.gtf \
		data/retentostat/mapped/${sample}Aligned.sortedByCoord.out.bam &
done < <(awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv)
echo "Waiting for stringtie assembly"
wait
echo "Assembly finished"
```


### Merge individual stringtie assemblies and compare to previous annotation. 

In this code the stringtie is separated into transcripts originating from annotated protein coding genes and
non-coding RNAs


```bash
./scripts/stringtie_merge.sh -t transcriptome_assembly \
	-g reference_genome/retentostat/cbs7435_retentostat.gff3 \
	-m transcriptome_assembly/stringtie_original.gtf \
	-r reference_genome/retentostat
```


# Salmon quantification

## Generate transcripts from stringtie transcriptome for use in Salmon. 
## Additionally, append whole genome sequence for use in decoy
```bash
gffread -w reference_genome/retentostat/cbs7435_retentostat_transcripts.fa \
	    -g reference_genome/retentostat/cbs7435_retentostat_combined.fa \
		transcriptome_assembly/stringtie.all.transcripts.gtf 

cat reference_genome/retentostat/cbs7435_retentostat_transcripts.fa \
    reference_genome/retentostat/cbs7435_retentostat_combined.fa \
	> reference_genome/retentostat/cbs7435_retentostat_transcripts_w_decoy.fa
```

### Generate decoy txt file for salmon - contains chromosome names

```bash
grep "^>" reference_genome/retentostat/cbs7435_retentostat_combined.fa \
	| cut -d ">" -f 2 \
	> reference_genome/retentostat/decoys.txt  
```

## lncRNA annotation

```bash
mkdir lncrna_annotation
```

## Salmon
### 1. Build index opt: --keepDuplicates
```bash
salmon index \
	-t reference_genome/retentostat/cbs7435_retentostat_transcripts_w_decoy.fa \
	-i data/retentostat/salmon/retentostat_index \
	-k 31 \
	--decoys reference_genome/retentostat/decoys.txt 
```

### 2. Quant

```bash
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv |  while read sample; do
	salmon quant --numBootstraps 100 \
		--gcBias \
		--seqBias \
		--threads 32 \
		--validateMappings \
		-l ISR  \
		-g transcriptome_assembly/stringtie.all.transcripts.gtf  \
		-i data/retentostat/salmon/retentostat_index \
		-1 data/retentostat/preprocessed/paired/${sample}_1.fastq.gz \
		-2 data/retentostat/preprocessed/paired/${sample}_2.fastq.gz  \
		-o data/retentostat/salmon/quant/${sample} 
done
```

### Merge TPM of quant files for TPM cutoff in lncRNA filtering
```bash
salmon quantmerge --output lncrna_annotation/TPM/transcript_tpm_all_samples.tsv \
	--quants data/retentostat/salmon/quant/* 

```


## FEELnc requires a GTF. Also want to generate fasta sequences of all annotated features in reference genome, separated into coding/non-coding.
-C coding only
-T GTF output
-F Keep all GFF attributes
--keep-genes: in transcript only mode, preserve gene records
--force-exons: make sure lowest level featires are "exon"

```bash
gffread reference_genome/cbs7435.gff3 -C -T -F --force-exons --keep-genes -o reference_genome/cbs7435_codingtranscripts.gtf
gffread reference_genome/cbs7435.gff3 --nc -T -F --keep-genes -o reference_genome/cbs7435_noncodingtranscripts.gtf
mkdir lncrna_annotation/data_fasta
gffread -w lncrna_annotation/data_fasta/non_coding_cdna_test.fa -g reference_genome/cbs7435_combined.fa reference_genome/cbs7435_noncodingtranscripts.gtf
gffread -w lncrna_annotation/data_fasta/coding_cdna_test.fa -g reference_genome/cbs7435_combined.fa reference_genome/cbs7435_codingtranscripts.gtf
```

## Download datasets for coding and non-coding features for use in coding potential predictions
Concatenate the extracted coding and non-coding features identified from Pichia. Use awk 1 to ensure new line is inserted
```bash
curl http://www.rnabinding.com/CPPred/Integrated-testing/Integrated_coding_RNA_test.fa -o lncrna_annotation/data_fasta/Integrated_coding_RNA_test.fa
curl http://www.rnabinding.com/CPPred/Integrated-testing/Integrated_ncrna_test.fa -o lncrna_annotation/data_fasta/Integrated_ncRNA_test.fa

awk 1 lncrna_annotation/data_fasta/Integrated_coding_RNA_test.fa lncrna_annotation/data_fasta/coding_cdna_test.fa > lncrna_annotation/data_fasta/combined_coding_rna_test.fa
awk 1 lncrna_annotation/data_fasta/Integrated_ncRNA_test.fa lncrna_annotation/data_fasta/non_coding_cdna_test.fa > lncrna_annotation/data_fasta/combined_non_coding_rna_test.fa
```


### FEELNc analysis
```bash
export FEELNCPATH=bin/FEELnc
mkdir -p lncrna_annotation/FEELnc
FEELnc_filter.pl \
	-i transcriptome_assembly/non_protein_coding_stringtie.gtf \
	-a reference_genome/cbs7435_codingtranscripts.gtf \
	-o lncrna_annotation/FEELnc/file.log \
	--monoex=1 \
	--size=200 \
	-p 32 \
	> lncrna_annotation/FEELnc/candidate_lncRNA.gtf

FEELnc_codpot.pl \
	-i lncrna_annotation/FEELnc/candidate_lncRNA.gtf \
	-m shuffle \
	-a reference_genome/cbs7435_codingtranscripts.gtf \
	-g reference_genome/cbs7435_combined.fa \
	--verbosity=0 \
	--keeptmp \
	--outdir lncrna_annotation/FEELnc/feelnc_codpot_out/ 

FEELnc_codpot.pl \
	-i lncrna_annotation/FEELnc/candidate_lncRNA.gtf \
	-l lncrna_annotation/data_fasta/combined_non_coding_rna_test.fa \
	-a lncrna_annotation/data_fasta/combined_coding_rna_test.fa \
	-g reference_genome/cbs7435_combined.fa \
	--verbosity=0 \
	--keeptmp \
	--outdir lncrna_annotation/FEELnc/feelnc_codpot_out_integrated/ 
```


### Assess FEELNc output using additional protein potential calculators, PFAM search and BLAST against protein and RNA databases
#### Use transdecoder to create a cDNA fasta file identify the longest ORF for each candidate lncRNA
##### Note: Produces pipeliners file (executeable command for pipeliner library) in main directory.

```bash
mkdir -p lncrna_annotation/TRANSDECODER
~/lncRNA_PpNZ/bin/TransDecoder/util/gtf_genome_to_cdna_fasta.pl \
	lncrna_annotation/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf \
	reference_genome/cbs7435_combined.fa \
	> lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa

~/lncRNA_PpNZ/bin/TransDecoder/TransDecoder.LongOrfs \
	-t lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa \
	-S \
	-O lncrna_annotation/TRANSDECODER
```

## Coding potential predictions

### Testing and building models

### Test RNAsamba, CPAT and CPC2 on K. phaffii data
```bash
./scripts/test_classifiers.sh -g lncrna_annotation/data_fasta/coding_cdna_test.fa \
	-n lncrna_annotation/data_fasta/non_coding_cdna_test.fa \
	-l lncrna_annotation
```



#### Classification with models
```bash
./scripts/run_rnasamba.sh \
	-f lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa  \
	-o lncrna_annotation/RNAsamba
```

### CPAT coding prediction for FEELNc candidate lncRNAs 
Maybe not run cpat predictions as would have to use mouse hexamers
Determine protein coding potential using CPAT
```bash
./scripts/run_cpat.sh \
	-f lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa \
	-o lncrna_annotation/CPAT \
	-c lncrna_annotation/data_fasta/combined_coding_rna_test.fa \
	-n lncrna_annotation/data_fasta/combined_non_coding_rna_test.fa 
Rscript ./scripts/10Fold_CrossValidation.r \
	"lncrna_annotation/CPAT/combined.feature.xls" \
	"lncrna_annotation/CPAT/10Fold_CV"
```



### CPC2 coding prediction for FEELnc candidate lncRNAs
Determine protein coding potential using CPC2

```bash
bin/CPC2_standalone-1.0.1/bin/CPC2.py \
  -i lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa \
  -o lncrna_annotation/CPC2/CPC2.analysis
```


## Database searches
### Protein domains - PFAM
Search for known PFAM domains
```bash
./scripts/run_hmmscan.sh \
	-t 4 \
	-e 1e-5 \
	-p lncrna_annotation/TRANSDECODER/longest_orfs.pep \
	-o lncrna_annotation/PFAM
```

### Blast against SWISSPROT, RFAM and MIRBASE databases
```bash
./scripts/run_blast.sh \
	-t 32 \
	-e 1e-5 \
	-s lncrna_annotation/SWISSPROT \
	-p lncrna_annotation/TRANSDECODER/longest_orfs.pep \
	-r lncrna_annotation/RNACentral \
	-n lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa
```

### RNAcentral
python scripts/search_RNAcentral.py lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa lncrna_annotation/RNACentral 1e-5 90


## Filtering
### First stage - take union of all analysis
lncrna_annotation/FEELnc/feelnc_codpot_out_integrated/candidate_lncRNA.gtf.lncRNA.gtf"
```bash
Rscript ./scripts/filter_lncrna.R \
	lncrna_annotation/FEELnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf \
	lncrna_annotation/CPC2/CPC2.analysis.txt \
	lncrna_annotation/CPAT/CPAT.analysis.ORF_prob.tsv \
	lncrna_annotation/RNAsamba/RNAsamba_classification_ensemble.tsv \
	lncrna_annotation/SWISSPROT/blastp.outfmt6 \
	lncrna_annotation/RNACentral/RNACentral.tsv \
	lncrna_annotation/RNACentral/Rfam_Families.tsv \
	lncrna_annotation/PFAM/pfam_domain_transcripts \
	lncrna_annotation/TPM/transcript_tpm_all_samples.tsv \
	lncrna_annotation/FEELnc/lncRNA_classes.txt \
	transcriptome_assembly/non_protein_coding_stringtie.gtf \
	reference_genome/cbs7435_combined.fa \
	lncrna_annotation
```

### Second stage - classification 
```bash
./scripts/filter_monoexonic_lncrnas.sh \
	-o lncrna_annotation \
	-g reference_genome/cbs7435_codingtranscripts.gtf \
	-l lncrna_annotation/firstpass_filter/all_lncrna.gtf \
	-a transcriptome_assembly/non_protein_coding_stringtie.gtf \
	-r transcriptome_assembly/ref_gene_id_name.tsv
```
#### Classify lncRNA types
```bash
Rscript ./scripts/simplify_class.R \
	lncrna_annotation/classification/first_classification.txt \
	lncrna_annotation/classification/simple_classification.txt

Rscript ./scripts/simplify_class.R \
	lncrna_annotation/classification/conservative_classification.txt \
	lncrna_annotation/classification/conservative_simple_classification.txt
```


### lncRNA characteristics
#### Want to generate fasta files and bed files of intergenic regions etc
```bash
samtools faidx \
	reference_genome/cbs7435_combined.fa

cut -f1,2 reference_genome/cbs7435_combined.fa.fai > reference_genome/cbs7435chromSizes.txt

cat transcriptome_assembly/stringtie.all.transcripts.gtf | \
	awk '$1 ~ /^#/ {print $0;next} {print $0 | \
	"sort -k1,1 -k4,4n -k5,5n"}' \
	> transcriptome_assembly/stringtie.all.transcripts.sorted.gtf

grep -v ^vhh_plasmid transcriptome_assembly/stringtie.all.transcripts.sorted.gtf \
	> transcriptome_assembly/stringtie.all.transcripts.sorted.vhhdropped.gtf

bedtools complement \
	-i transcriptome_assembly/stringtie.all.transcripts.sorted.vhhdropped.gtf \
	-g reference_genome/cbs7435chromSizes.txt \
	> transcriptome_assembly/cbs7435_intergenic.bed

bedtools getfasta \
	-fi reference_genome/cbs7435_combined.fa \
	-bed transcriptome_assembly/cbs7435_intergenic.bed \
	-fo lncrna_annotation/data_fasta/intergenic_dna.fa

gffread transcriptome_assembly/stringtie.all.transcripts.gtf \
	--bed \
	-o transcriptome_assembly/stringtie.all.transcripts.bed
```

## Functional predictions
### Combine original and lncRNA annotation
```bash
gffread reference_genome/cbs7435.gff3 \
	--force-exons \
	--sort-alpha \
	--keep-genes \
	-T \
	-o reference_genome/CBS7435_updated.gtf

awk 1 reference_genome/CBS7435_updated.gtf \
	lncrna_annotation/all_lncrna.gtf \
	> reference_genome/CBS7435_with_lncrna.gtf

gffread reference_genome/CBS7435_with_lncrna.gtf \
	--force-exons \
	--sort-alpha \
	--keep-genes \
	-T \
	-o reference_genome/CBS7435_updated.gtf

gffread lncrna_annotation/all_lncrna.gtf \
	-w lncrna_annotation/lncrna_sequences.fa \
	-g reference_genome/cbs7435_combined.fa
```

#### Plot characteristics and details about lncRNAs, other genes etc.
```bash
mkdir lncrna_annotation/plots
./scripts/gaps_to_GTF.py reference_genome/cbs7435_combined.fa > reference_genome/cbs7435_gaps.gtf
./scripts/plot_chromosomes.py \
	reference_genome/cbs7435_codingtranscripts.gtf \
	lncrna_annotation/all_lncrna.gtf \
	reference_genome/cbs7435_gaps.gtf \
	reference_genome/cbs7435chromSizes.txt \
	lncrna_annotation/plots

./scripts/plot_characteristics.py  \
	lncrna_annotation/lncrna_sequences.fa \
	lncrna_annotation/data_fasta/non_coding_cdna_test.fa \
	lncrna_annotation/data_fasta/intergenic_dna.fa \
	lncrna_annotation/data_fasta/coding_cdna_test.fa \
	lncrna_annotation/monoexonic_filter/first_lncrna_list.txt \
	reference_genome/other.noncoding.genes.txt \
	reference_genome/retentostat/protein.coding.genes.txt \
	transcriptome_assembly/stringtie.all.transcripts.bed \
	lncrna_annotation/TPM/transcript_tpm_all_samples.tsv \
	lncrna_annotation/plots

Rscript ./scripts/plot_genomealignments.R 
```


## Differential expression analysis
### DESeq2 differential expression - Salmon
```bash
Rscript ./scripts/deseq_analysis.R \
	"data/retentostat/salmon/quant/" \
	"data/retentostat_sample_information.tsv" \
	"lncrna_annotation/classification/simple_classification.txt" \
	"lncrna_annotation/monoexonic_filter/first_lncrna_list.txt" \
	"transcriptome_assembly/stringtie.all.transcripts.gtf" \
	"differential_expression/retentostat" \
	"lncrna_annotation/TPM/TPM_ISAER.tsv" \
	"data/kphaffi_annotation_info.tsv" \
	"functional_predictions/triplexes/interacting_pairs_coding.txt" \
	"functional_predictions/nearby/nearby_genes.txt" \
	"functional_predictions"
```

```bash
gffread lncrna_annotation/all_lncrna.gtf \
	--keep-genes \
	-F \
	-o lncrna_annotation/all_lncrna.gff3
	
cat reference_genome/cbs7435.gff3 \
	lncrna_annotation/all_lncrna.gff3 \
	| gffread -F \
		--keep-genes \
		-o reference_genome/cbs7435_with_lncrna.gff3
```


### Determine pairs of lncRNA and genes that interact via triplexes or are within ~1kb of each other
```bash
mkdir -p functional_predictions
./scripts/predict_interactions.sh \
	-l lncrna_annotation/lncrna_sequences.fa \
	-f cbs7435_combined.fa \
	-g reference_genome/CBS7435_updated.gtf \
	-s transcriptome_assembly/stringtie.all.transcripts.gtf \
	-r reference_genome/retentostat \
	-o functional_predictions \
	-k lncrna_annotation/all_lncrna.gtf \
	-c cbs7435chromSizes.txt \
	-p reference_genome/protein.coding.genes.txt
```



### Download and prepare the GO/KEGG term analysis - GS115 to CBS7435
```bash
./scripts/replace_go.sh \
	-o reference_genome/GO_analysis \
	-i data/cbs7435_gs115.tsv \
wget --quiet -nc \
	-O reference_genome/GO_analysis/kp_kegg.json \
	'https://www.genome.jp/kegg-bin/download_htext?htext=ppa00001&format=json&filedir=kegg/brite/ppa'
./scripts/format_kegg.py
```

## Enrichment of triplex pairs
```bash
Rscript ./scripts/triplex_enrichment.R \
	functional_predictions/triplexes/interacting_pairs_coding.txt \
	functional_predictions/triplexes/pair_enrichment.xlsx
```

### Determine correlation coefficients, generate density distributions, GO/KEGG analysis and clustering 
```bash
Rscript ./scripts/retentostat_functional.R \
	"data/retentostat/salmon/quant/" \
	"data/retentostat_sample_information.tsv" \
	"lncrna_annotation/classification/simple_classification.txt" \
	"lncrna_annotation/monoexonic_filter/first_lncrna_list.txt" \
	"transcriptome_assembly/stringtie.all.transcripts.gtf" \
	"functional_predictions/triplexes/interacting_pairs_coding.txt" \
	"functional_predictions/nearby/nearby_genes.txt" \
	"functional_predictions"
```


```bash
Rscript ./scripts/clustering_assesment.R \
	functional_predictions/clustering \
	transcriptome_assembly/stringtie.all.transcripts.gtf \
	lncrna_annotation/monoexonic_filter/first_lncrna_list.txt
```
