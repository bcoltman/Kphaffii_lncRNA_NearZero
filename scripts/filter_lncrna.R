#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: filter_lncrna.R
## Purpose: Filter FEELnc-identified lncRNA transcripts using multiple coding potential metrics
## Author: Adapted script
## Date Created: Dec-2020
##
## -----------------------------------------------------------------------------

suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))

# Command-line arguments
args <- commandArgs(TRUE)
candidate_lncrna_gtf <- args[1]
cpc_output <- args[2]
cpat_output <- args[3]
rnasamba_output <- args[4]
swissprot_blast_output <- args[5]
rna_central_output <- args[6]
rfam_rna_central_output <- args[7]
pfam_domain_hit_output <- args[8]
tpm_expression_matrix <- args[9]
feelnc_lncrna_classification <- args[10]
stringtie_gtf <- args[11]
reference_genome <- args[12]
output_dir <- args[13]

# Create output directory
output_path <- file.path(output_dir, "firstpass_filter")
dir.create(output_path, showWarnings = FALSE)



# Load FEELnc GTF, removing mtDNA transcripts
feelnc_gtf <- rtracklayer::import(candidate_lncrna_gtf)
feelnc_gtf <- feelnc_gtf[seqnames(feelnc_gtf) != "MT"]
feelnc_lncrna_transcripts <- unique(feelnc_gtf$transcript_id)

write("---------------FEELNC-------------------", stdout())
write(paste(length(feelnc_lncrna_transcripts), "candidate lncRNAs identified by FEELNc"), stdout())

# Load coding potential data from CPC2, CPAT, and RNAsamba
cpc2_lncrna_transcripts <- subset(read.table(cpc_output, sep = "\t"), V8 == "noncoding")$V1
write(paste(length(cpc2_lncrna_transcripts), "CPC2 Candidates"),stdout())


cpat <- read.table(cpat_output, sep = "\t", header = TRUE)
cpat_lncrna_transcripts <- setdiff(
  unique(gsub("_.*$", "", subset(cpat, Coding_prob < 0.4)$ID)),
  unique(gsub("_.*$", "", subset(cpat, Coding_prob > 0.4)$ID))
)
write(paste(length(cpat_lncrna_transcripts), "CPAT Candidates"),stdout())

rnasamba <- subset(read.table(rnasamba_output, sep = "\t", header = TRUE), classification == "noncoding")$sequence_name
rnasamba_lncrna_transcripts <- gsub(" .*$", "", rnasamba)
write(paste(length(rnasamba_lncrna_transcripts), "RNAsamba Candidates"),stdout())

# Overlap results of the coding potential calculators
# lncrna_all_noncoding <- intersect(
  # feelnc_lncrna_transcripts,
  # intersect(cpat_lncrna_transcripts, rnasamba_lncrna_transcripts)
# )

lncrna_all_noncoding <- intersect(
  feelnc_lncrna_transcripts, rnasamba_lncrna_transcripts)


# overlap results of 3 coding potential calculators

# write("---------------CPAT, RNAsamba-------------------", stdout())
# write(paste(length(intersect(rnasamba_lncrna_transcripts, feelnc_lncrna_transcripts)), "candidate lncRNAs common to FEELNc, RNAsamba"), stdout())
write(paste(length(intersect(cpc2_lncrna_transcripts, feelnc_lncrna_transcripts)), "candidate lncRNAs common to FEELNc, CPC2"), stdout())
write(paste(length(intersect(cpat_lncrna_transcripts, feelnc_lncrna_transcripts)), "candidate lncRNAs common to FEELNc, CPAT"), stdout())
# write(paste(length(lncrna_all_noncoding), "candidate lncRNAs common to FEELNc, RNAsamba & CPAT"), stdout())
write(paste(length(lncrna_all_noncoding), "candidate lncRNAs common to FEELNc & RNAsamba"), stdout())



# Filter results with SWISSPROT, PFAM, RNAcentral, RFAM, and TPM
write("---------------Filtering with External Databases-------------------", stdout())

# SWISSPROT
if (file.exists(swissprot_blast_output)) {
  swissprot_hits <- unique(substr(read.table(swissprot_blast_output, sep = "\t")$V1, 1, nchar(swissprot_blast_output) - 3))
  lncrna_all_noncoding <- lncrna_all_noncoding[!(lncrna_all_noncoding %in% swissprot_hits)]
}
write(paste(length(swissprot_hits), "SWISSPROT hits"), stdout())
write(paste(length(lncrna_all_noncoding), "remaining after SWISSPROT filtering"), stdout())

# PFAM
if (file.exists(pfam_domain_hit_output)) {
  pfam_hits <- unique(read.table(pfam_domain_hit_output)$V1)
  lncrna_all_noncoding <- lncrna_all_noncoding[!(lncrna_all_noncoding %in% pfam_hits)]
}
write(paste(length(pfam_hits), "PFAM hits"), stdout())
write(paste(length(lncrna_all_noncoding), "remaining after PFAM domain filtering"), stdout())

# RNAcentral and RFAM filtering
if (file.exists(rna_central_output)) {
  rna_central_hits <- unique(read.table(rna_central_output, sep = "\t", header = TRUE)$Transcript)
  #lncrna_all_noncoding <- lncrna_all_noncoding[!(lncrna_all_noncoding %in% rna_central_hits)]
}
write(paste(length(rna_central_hits), "RNAcentral hits"), stdout())
write(paste(length(lncrna_all_noncoding), "remaining after RNA central filtering"), stdout())

if (file.exists(rfam_rna_central_output)) {
  rfam_hits <- unique(read.table(rfam_rna_central_output, sep = "\t", header = TRUE)$Transcript)
  lncrna_all_noncoding <- lncrna_all_noncoding[!(lncrna_all_noncoding %in% rfam_hits)]
}
write(paste(length(rfam_hits), "RFAM hits"), stdout())
write(paste(length(lncrna_all_noncoding), "remaining after RFAM-RNA central filtering"), stdout())


write("---------------lncRNA - RNA Central-------------------", stdout())


# TPM filter
tpm_matrix <- read.table(tpm_expression_matrix, sep = "\t", header = TRUE, row.names = 1)
lncrna_tpm_matrix <- tpm_matrix[rownames(tpm_matrix) %in% lncrna_all_noncoding, ]
lncrna_transcripts <- rownames(lncrna_tpm_matrix[apply(lncrna_tpm_matrix, 1, max) > 1, ])
write(paste(length(lncrna_transcripts), "candidates with TPM > 1 in at least one sample"), stdout())

# Save filtered transcripts and other outputs
stringtie_txdb <- suppressMessages(makeTxDbFromGFF(stringtie_gtf))
lncrna_genes <- suppressMessages(select(stringtie_txdb, keys = lncrna_transcripts, columns = "GENEID", keytype = "TXNAME"))

save(lncrna_transcripts, lncrna_genes, file = file.path(output_path, "lncRNA_filtering.rData"))

write("Writing final lists of transcripts", stdout())
write(lncrna_transcripts, file.path(output_path, "all_lncrna_transcripts.txt"))

# Export final GTF and FASTA
system(paste("grep -wFf", file.path(output_path, "all_lncrna_transcripts.txt"), stringtie_gtf, ">", file.path(output_path, "all_lncrna.gtf")))
system(paste("~/lncRNA_PpNZ/bin/TransDecoder/util/gtf_genome_to_cdna_fasta.pl", file.path(output_path, "all_lncrna.gtf"), reference_genome, ">", file.path(output_path, "lncRNA.fasta")))

write("Completed: Data saved", stdout())
write(paste("Novel lncRNA transcripts identified:", length(lncrna_transcripts)), stdout())
