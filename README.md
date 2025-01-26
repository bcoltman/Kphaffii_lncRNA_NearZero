# Identification of lncRNAs in _Komagataella phaffii_ and their responses to extremely slow growth rates. 

## This repository contains the code for analysis of:
## Growth-rate associated expression levels of long non-coding RNAs in extremely slow growing _Komagataella phaffii_, Coltman _et_al_, 2025

## DOI
Our most recent Zenodo DOI generated for the journal submission is: 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14742669.svg)](https://doi.org/10.5281/zenodo.14742669)

## Contact
- benjamin.coltman@univie.ac.at

## Downloading and reproducing results

The results presented in Coltman, 2025 can be reproduced by cloning this repository, defining a conda environment and executing the included scripts. Raw data required for running these scripts is either made avilable for download in [data](data) or is downloaded by the scripts. The programmes and the versions used can be found in [requirements.txt](requirements.txt).

## Work notebook for project:

### 1. Download reads and quality control
#### a. Downloading reads from NCBI SRA
 Downloading SRA files - need to include a file from NCBI with the SRA data attached to it to then loop through and use SRR number in below format. Can also use it for study design too

```bash
mkdir -p data/retentostat/sra
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do
prefetch $sample  -O data/retentostat/sra
fasterq-dump data/retentostat/sra/$sample \
	--threads 32 \
	--split-files \
	--outdir data/retentostat/raw
pigz -p 32 data/retentostat/raw/${sample}_1.fastq 
pigz -p 32 data/retentostat/raw/${sample}_2.fastq 
done
```

#### b. Initial QC of raw reads
```bash
mkdir -p QC/retentostat/multiQC QC/retentostat/raw
fastqc -q -t 32 -o QC/retentostat/raw data/retentostat/raw/*.gz
multiqc QC/retentostat/raw -n raw -o QC/retentostat/multiQC
```

#### c. Cutadapt and QC
```bash
mkdir -p data/retentostat/cutadapt \
	QC/retentostat/post_cut
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}; next}{ print $(f["accession"])}' data/retentostat_sample_information.tsv | while read sample; do
cutadapt  \
	-j 32 \
	-m 1 \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o data/retentostat/cutadapt/${sample}_1.fastq.gz \
	-p data/retentostat/cutadapt/${sample}_2.fastq.gz \
	   data/retentostat/raw/${sample}_1.fastq.gz \
	   data/retentostat/raw/${sample}_2.fastq.gz \
	   1>data/retentostat/cutadapt/${sample}_report.txt 
done
fastqc -q -t 32 -o QC/retentostat/post_cut data/retentostat/cutadapt/*.gz
multiqc QC/retentostat/post_cut -n post_cut -o QC/retentostat/multiQC
```

#### d. Trimmomatic and QC
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
mkdir QC/retentostat/post_trim
fastqc -q -t 32 -o QC/retentostat/post_trim data/retentostat/preprocessed/paired/*.gz
multiqc QC/retentostat/post_trim -n post_trim -o QC/retentostat/multiQC
```

### 2. Mapping reads to genome and generating transcriptome
#### a. Download and prepare reference genome for CBS7435
```bash
./scripts/prepare_genome.sh -o reference_genome
```

#### b. Make STAR Index of CBS7435 reference genome
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

#### c. Map reads to STAR index CBS7435 reference enome
```bash
mkdir -p data/retentostat/mapped
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

#### d. Reference-assisted transcriptome assembly using Stringtie
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

#### e. Merge individual stringtie assemblies and compare changes
```bash
./scripts/stringtie_merge.sh -t transcriptome_assembly \
	-g reference_genome/retentostat/cbs7435_retentostat.gff3 \
	-m transcriptome_assembly/stringtie_original.gtf \
	-r reference_genome/retentostat
```

### 3. Quantify transcript abundance using Salmon
#### a. Prepare for generating Salmon index 
##### i. Extract transcript annotations and add intergenic regions for decoy
```bash
gffread -w reference_genome/retentostat/cbs7435_retentostat_transcripts.fa \
	    -g reference_genome/retentostat/cbs7435_retentostat_combined.fa \
		transcriptome_assembly/stringtie.all.transcripts.gtf 

cat reference_genome/retentostat/cbs7435_retentostat_transcripts.fa \
    reference_genome/retentostat/cbs7435_retentostat_combined.fa \
	> reference_genome/retentostat/cbs7435_retentostat_transcripts_w_decoy.fa
```

##### ii. Generate decoy txt file for Salmon containing chromosome names
```bash
grep "^>" reference_genome/retentostat/cbs7435_retentostat_combined.fa \
	| cut -d ">" -f 2 \
	> reference_genome/retentostat/decoys.txt  
```

#### b. Build Salmon index
```bash
salmon index \
	-t reference_genome/retentostat/cbs7435_retentostat_transcripts_w_decoy.fa \
	-i data/retentostat/salmon/retentostat_index \
	-k 31 \
	--decoys reference_genome/retentostat/decoys.txt 
```

#### c. Quantify counts with Salmon
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

#### d. Merge TPM counts fo all samples (only used for filtering lncRNAs)
```bash
mkdir lncrna_annotation
salmon quantmerge --output lncrna_annotation/TPM/transcript_tpm_all_samples.tsv \
	--quants data/retentostat/salmon/quant/* 
```

### 4. Initial predictions of lncRNAs
#### a. Prepare existing annotations download datasets
##### i. Extract sequences of all existing coding and non-coding annotations for use in lncRNA prediction tools 
```bash
mkdir lncrna_annotation/data_fasta
gffread reference_genome/cbs7435.gff3 -C \
	-T \
	-F \
	--force-exons \
	--keep-genes \
	-o reference_genome/cbs7435_codingtranscripts.gtf
gffread reference_genome/cbs7435.gff3 \
	--nc \
	-T \
	-F \
	--keep-genes \
	-o reference_genome/cbs7435_noncodingtranscripts.gtf
gffread -w lncrna_annotation/data_fasta/non_coding_cdna_test.fa \
	-g reference_genome/cbs7435_combined.fa \
	reference_genome/cbs7435_noncodingtranscripts.gtf
gffread -w lncrna_annotation/data_fasta/coding_cdna_test.fa \
	-g reference_genome/cbs7435_combined.fa \
	reference_genome/cbs7435_codingtranscripts.gtf
```

##### ii. Download coding and non-coding features dataset from CPPred for use in coding potential predictions
Concatenate the extracted coding and non-coding features identified from Pichia. Use awk 1 to ensure new line is inserted
```bash
curl http://www.rnabinding.com/CPPred/Integrated-testing/Integrated_coding_RNA_test.fa \
	-o lncrna_annotation/data_fasta/Integrated_coding_RNA_test.fa
curl http://www.rnabinding.com/CPPred/Integrated-testing/Integrated_ncrna_test.fa \
	-o lncrna_annotation/data_fasta/Integrated_ncRNA_test.fa
awk 1 lncrna_annotation/data_fasta/Integrated_coding_RNA_test.fa \
	lncrna_annotation/data_fasta/coding_cdna_test.fa \
	> lncrna_annotation/data_fasta/combined_coding_rna_test.fa
awk 1 lncrna_annotation/data_fasta/Integrated_ncRNA_test.fa \
	lncrna_annotation/data_fasta/non_coding_cdna_test.fa \
	> lncrna_annotation/data_fasta/combined_non_coding_rna_test.fa
```

#### b. Predict lncRNAs with FEELnc
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

#### c. Identify the longest potential ORFs from FEELnc predicted lncRNAs to use for PFAM search and BLAST against protein and RNA databases
 Note: Produces pipeliners file (executeable command for pipeliner library) in main directory.
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

### 5. Filtering FEELnc predicted lncRNAs 
#### a. Assess additional coding potential prediction tools (RNAsamba, CPAT and CPC2) on pre-existing _K. phaffii_ annotations

We observed that RNAsamba predictions vary slightly between runs. Subsequently, we include the pre-trained models in "/data" which are copied during code execution. The model training code is within the scripts and can be uncommented.
```bash
./scripts/test_classifiers.sh -g lncrna_annotation/data_fasta/coding_cdna_test.fa \
	-n lncrna_annotation/data_fasta/non_coding_cdna_test.fa \
	-l lncrna_annotation
```

#### b. Classify lncRNAs with RNAsamba, CPAT and CPC2
##### i. RNAsamba
```bash
./scripts/run_rnasamba.sh \
	-f lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa  \
	-o lncrna_annotation/RNAsamba
```

##### ii. CPAT (not used in final filtering, but included for comparison)
Uses mouse hexamers
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

##### iii. CPC2 (not used in final filtering, but included for comparison)
```bash
bin/CPC2_standalone-1.0.1/bin/CPC2.py \
  -i lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa \
  -o lncrna_annotation/CPC2/CPC2.analysis
```

#### c. Assess the longest ORFs for possible PFAM domains
```bash
./scripts/run_hmmscan.sh \
	-t 4 \
	-e 1e-5 \
	-p lncrna_annotation/TRANSDECODER/longest_orfs.pep \
	-o lncrna_annotation/PFAM
```

#### d. Blast against SWISSPROT databases
```bash
./scripts/run_blast.sh \
	-t 32 \
	-e 1e-5 \
	-s lncrna_annotation/SWISSPROT \
	-p lncrna_annotation/TRANSDECODER/longest_orfs.pep
```

#### e. Query RNAcentral with lncRNA sequences
```bash
mkdir lncrna_annotation/RNACentral
python scripts/search_RNAcentral.py lncrna_annotation/TRANSDECODER/candidate_lncrna_cdna.fa lncrna_annotation/RNACentral 1e-5 90
```

#### f. First filtering stage - take union of all analysis
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

#### g. Second stage
##### i. Filtering and classification with respect to other genes
```bash
./scripts/filter_monoexonic_lncrnas.sh \
	-o lncrna_annotation \
	-g reference_genome/cbs7435_codingtranscripts.gtf \
	-l lncrna_annotation/firstpass_filter/all_lncrna.gtf \
	-a transcriptome_assembly/non_protein_coding_stringtie.gtf \
	-r transcriptome_assembly/ref_gene_id_name.tsv
```

##### ii. Simplification of classifications
```bash
Rscript ./scripts/simplify_class.R \
	lncrna_annotation/classification/first_classification.txt \
	lncrna_annotation/classification/simple_classification.txt
```

### 6. Assessing identified lncRNAs, characteristics and predicting potential functions
#### a. Prepare files for comparing characteristics of lncRNAs to other coding and non-coding genes and intergenic regions
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

#### b. Combine original annotation with filtered lncRNA annotations
```bash
gffread reference_genome/cbs7435.gff3 \
	--force-exons \
	--sort-alpha \
	--keep-genes \
	-T \
	-o reference_genome/cbs7435_updated.gtf
awk 1 reference_genome/cbs7435_updated.gtf \
	lncrna_annotation/firstpass_filter/all_lncrna.gtf \
	| gffread --force-exons \
	--sort-alpha \
	--keep-genes \
	-T \
	-o reference_genome/cbs7435_updated.gtf
gffread lncrna_annotation/firstpass_filter/all_lncrna.gtf \
	-w lncrna_annotation/lncrna_sequences.fa \
	-g reference_genome/cbs7435_combined.fa
```


#### c. Plot characteristics and details about lncRNAs, other genes etc.
```bash
mkdir lncrna_annotation/plots
./scripts/gaps_to_GTF.py reference_genome/cbs7435_combined.fa > reference_genome/cbs7435_gaps.gtf
./scripts/plot_chromosomes.py \
	reference_genome/cbs7435_codingtranscripts.gtf \
	lncrna_annotation/firstpass_filter/all_lncrna.gtf \
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
```

#### d. Determine pairs of lncRNA and genes that interact via triplexes or are within ~1kb of each other
```bash
mkdir -p functional_predictions
./scripts/predict_interactions.sh \
	-l lncrna_annotation/lncrna_sequences.fa \
	-f cbs7435_combined.fa \
	-g reference_genome/cbs7435_updated.gtf \
	-s transcriptome_assembly/stringtie.all.transcripts.gtf \
	-r reference_genome/retentostat \
	-o functional_predictions \
	-k lncrna_annotation/firstpass_filter/all_lncrna.gtf \
	-c cbs7435chromSizes.txt \
	-p reference_genome/retentostat/protein.coding.genes.txt
```

#### e. Investigate enrichment of annotated functions (KEGG or GO) of coding genes predicted to interact with lncRNAs via triplexes
##### i. Download and prepare the GO/KEGG term analysis - GS115 to CBS7435
```bash
./scripts/replace_go.sh \
	-o reference_genome/GO_analysis \
	-i data/cbs7435_gs115.tsv \
wget --quiet -nc \
	-O reference_genome/GO_analysis/kp_kegg.json \
	'https://www.genome.jp/kegg-bin/download_htext?htext=ppa00001&format=json&filedir=kegg/brite/ppa'
./scripts/format_kegg.py
```

##### ii. Perform enrichment analysis
```bash
Rscript ./scripts/triplex_enrichment.R \
	functional_predictions/triplexes/interacting_pairs_coding.txt \
	functional_predictions/triplexes/pair_enrichment.xlsx
```

### 7. Assess growth-rate dependent expression of lncRNAs
#### DESeq2 differential expression using Salmon quantified abundances & functional assessments
##### Assess correlation coefficients of lncRNAs with other coding genes (neighbours, interacting pairs, or otherwise), co-clustering genes and enrichments (GO/KEGG) of functions in these sub compostions
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

### 8. Plot coverage plots of lncRNA
##### Split aligned reads by strand - indexing BAMs too

```bash
list=("SRR25915888" "SRR25915870" "SRR25915879" "SRR25915871" "SRR25915867" "SRR25915875")

for item in "${list[@]}"; do
    scripts/split_bam_file_by_strand.sh data/retentostat/mapped/${item}Aligned.sortedByCoord.out.bam $item data/retentostat/mapped_split_by_strand &
done


for item in "${list[@]}"; do
    samtools index -@ 6 data/retentostat/mapped/${item}Aligned.sortedByCoord.out.bam &
done
```
##### Plot coverage plots 

```bash
Rscript scripts/plot_coverage_lncRNA.R 
```



