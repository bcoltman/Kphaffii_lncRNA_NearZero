#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: retentostat_deseq.r
##
## Purpose of script: To carry out count based differential expression analyis 
## for the full StringTie transcriptome assembly. Uses IsoformSwitchAnalyzeR to correct stringtie issues (merging etc).
## The differential expression results for  lncRNAs are parsed and the classififcaiton and orthology information is added 
##
## -----------------------------------------------------------------------------
##
## Info: Genes differentially expressed with +/- 1.5 fold difference, along with a
## BH adjusted p-value and baseMean >=100 are considered significant
##
## -----------------------------------------------------------------------------

# load libraries
# load libraries



#library(data.table)
#suppressMessages(library("biomaRt"))
#suppressMessages(library("GenomicRanges"))
#suppressMessages(library("magrittr"))
#suppressMessages(library("stringr"))
#suppressMessages(library("ComplexHeatmap"))
#suppressMessages(library("gprofiler2"))


suppressMessages(library("IsoformSwitchAnalyzeR"))
suppressMessages(library("DESeq2"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("DEFormats"))
suppressMessages(library("edgeR"))
suppressMessages(library("pheatmap"))

options(bitmapType='cairo')

# Input arguments
args <- commandArgs(TRUE)
quant_dir <- args[1]
sample_info_file <- args[2]
lncrna_classification <- args[3]
lncrna_list <- args[4]
# "lncrna_annotation/monoexonic_filter/conservative_lncrna_list.txt", 
gtf <- args[5]
results_dir <- args[6]
tpm <- args[7]

dir.create(file.path(results_dir, "/volcano"),recursive=TRUE)
dir.create(file.path(results_dir, "/heatmaps"),recursive=TRUE)
dir.create(file.path(results_dir, "/PCA"),recursive=TRUE)

# Set default directories (for testing purposes)
#quant_dir <- "data/retentostat/salmon/quant/"
#sample_info_file <- "data/retentostat_sample_information.tsv"
#lncrna_classification <- "lncrna_annotation/classification/simple_classification.txt" 
#lncrna_list <- "lncrna_annotation/monoexonic_filter/first_lncrna_list.txt"
# "lncrna_annotation/monoexonic_filter/conservative_lncrna_list.txt", 
#gtf <- "transcriptome_assembly/stringtie.all.transcripts.gtf"
#results_dir <- "differential_expression/retentostat"
#tpm <- "lncrna_annotation/TPM/TPM_ISAER.tsv"


volcano_ggplot <- function(volcano_data, comparison, scm_values, results_dir) {
			vgp <- ggplot(volcano_data) +
			geom_point(aes(x = log2FoldChange, 
				       y = -log10(padj), 
				       color = pointclass, 
				       fill = pointclass, 
				       size = pointsize)) + 
			scale_size_identity() + # w/o ggplot scales marker size range between 1-6
			scale_color_manual(name="Colour", 
				values=scm_values, 
				limits=unique(volcano_data$pointclass), 
				labels=unique(volcano_data$pointclass))  +
			xlab(bquote(~ Log[2] ~ "fold change")) + 
			ylab("-log10 BH adjusted p-value") +
			geom_hline(
				yintercept = -log10(0.05),
				linetype = "dashed", color = "gray") + 
			#geom_vline(
			#    	xintercept = c(-0.5849625, 0.5849625),
			#    	linetype = "dashed",
			#    	color = "light gray") +  
			theme_bw() +
			theme(
			    	legend.position = "none",
			    	plot.title = element_text(size = rel(1.5), hjust = 0.5),
			    	axis.title = element_text(size = rel(1.25)),
			    	panel.grid.major = element_blank(),
			    	panel.grid.minor = element_blank(),
			    	strip.background = element_blank(),
			    	strip.text = element_text(face = "bold", size = 7)) 

			# save
			ggsave(paste0(results_dir,"/volcano/", comparison, "_lncrna_volcano.png"), 
			       plot = vgp, 
			       height = 5, width = 5, units = "in", dpi = 600)
}


classification <- function(sig_results, lncRNA_classification, results_dir) {

		# classification match columns
		match_classification <- lncRNA_classification[match(rownames(sig_results), lncRNA_classification$lncRNA_gene), "summary_class"]
		match_classification[is.na(match_classification)] <- "intergenic"

		match_partner <- lncRNA_classification[match(rownames(sig_results), lncRNA_classification$lncRNA_gene), "partnerRNA_gene"]
		match_partner[is.na(match_partner)] <- "intergenic"

		match_partner_name <- lncRNA_classification[match(rownames(sig_results), lncRNA_classification$lncRNA_gene), "partnerRNA_gene_name"]
		match_partner_name[is.na(match_partner_name)] <- "intergenic"
		
		match_distance <- lncRNA_classification[match(rownames(sig_results), lncRNA_classification$lncRNA_gene), "distance"]
		match_distance[is.na(match_distance)] <- "intergenic"

		# create a summary table

		lncrna_de_matrix <- data.frame(cbind(
		  sig_results[, c("baseMean","log2FoldChange","lfcSE","pvalue","padj")],
		  match_classification, 
		  match_distance,
		  match_partner,
		  match_partner_name
		))


		colnames(lncrna_de_matrix) <- c("baseMean",
		  "log2FoldChange",
		  "lfcSE",
		  "pvalue",
		  "padj",
		  "lncRNA type",
		  "Match distance",
		  "Closest gene",
		  "Closest gene name")

		return(lncrna_de_matrix)
		
		}



# volcano plot function
volcano_plot <- function(de_results) {
  #de_results <- de_results[de_results$baseMean >= 100, ]

  # identify significant
  de_results$threshold_DE <- de_results$padj < 0.05 
#  de_results$threshold_DE <- de_results$padj < 0.05 & 
#                             abs(de_results$log2FoldChange) >= 0.5849625
  volcano_data <- cbind(de_results$log2FoldChange, 
                    de_results$padj, 
                    as.logical(de_results$threshold_DE))
  
  
  # data frame for all required information
  volcano_data <- data.frame(
    ensemblid = rownames(de_results),
    log2FoldChange = de_results$log2FoldChange,
    padj = de_results$padj,
    threshold = de_results$threshold_DE)
  volcano_data <- na.omit(volcano_data)

  volcano_data <- volcano_data %>% dplyr::mutate(pointcolor = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ "#F9B288",
    log2FoldChange < 0 & threshold == 1 ~ "#A2D9F9",
    threshold == 0 ~ "gray"
  ))
  volcano_data <- volcano_data %>% dplyr::mutate(pointclass = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ "Upregualted",
    log2FoldChange < 0 & threshold == 1 ~ "Downregulated",
    threshold == 0 ~ "Not Significant"
  ))
  volcano_data <- volcano_data %>% dplyr::mutate(pointsize = dplyr::case_when(
    log2FoldChange > 0 & threshold == 1 ~ 2.5,
    log2FoldChange < 0 & threshold == 1 ~ 2.5,
    threshold == 0 ~ 1
  ))
  return(volcano_data)
}

process_comparison <- function(desired_comparisons, comparison_interpretation, DESeqDataSet, lncrna_genes, lncRNA_classification, file_name, results_dir) {
    # Create a blank workbook
    OUT <- createWorkbook()

    for (i in 1:length(desired_comparisons)) {
        comparison <- desired_comparisons[i]
        interpretation <- comparison_interpretation[i]

	res <- results(DESeqDataSet, alpha=0.05, name=comparison)
	#lfcThreshold=log2(1.5), altHypothesis="greaterAbs"
	shrunk_results <- lfcShrink(DESeqDataSet, coef=comparison, type = "apeglm", res=res)

	sig_de_results <- subset(
	    shrunk_results,
	    padj < 0.05
	)

#        suppressMessages(shrunk_results <- lfcShrink(DESeqDataSet, coef = comparison, type = "apeglm"))

        # set differential expression criteria
#        sig_de_results <- subset(
#            shrunk_results,
#            abs(log2FoldChange) >= 0.5849625 & padj < 0.05
#        )

        sig_de_results <- sig_de_results[
            order(sig_de_results$log2FoldChange, decreasing = TRUE),]

        #write.table(as.data.frame(sig_de_results),
        #            file = paste0(results_dir,"/", comparison, ".tsv"),
        #            quote = FALSE, col.names = NA, sep = "\t")

        # Export all results - all transcripts and both DE and Non DE
        addWorksheet(OUT, paste0(interpretation, "_all"))
        writeData(OUT, rowNames = TRUE, sheet = paste0(interpretation, "_all"),
                  x = as.data.frame(shrunk_results))

	# extract only lncRNA genes
        all_lncrna_results <- shrunk_results[rownames(shrunk_results) %in% lncrna_genes, ]
	
	# Plot volcano of lncRNA
	volcano_data <- volcano_plot(all_lncrna_results)
	scm_values <- unique(volcano_data$pointcolor)
	names(scm_values) <- unique(volcano_data$pointclass)
	volcano_ggplot(volcano_data, interpretation, scm_values, results_dir)

	# Add sheet of just lncRNA results
	lncrna_matrix <- classification(all_lncrna_results, lncRNA_classification, results_dir)
        addWorksheet(OUT, paste0(interpretation, "_all_lncRNA"))
        writeData(OUT, rowNames = TRUE, sheet = paste0(interpretation, "_all_lncRNA"),
                  x = as.data.frame(lncrna_matrix))
	
	# Now, want to add sheet of just DE lncRNA 
	sig_de_results_lncrna <- sig_de_results[rownames(sig_de_results) %in% lncrna_genes, ]

        if (nrow(sig_de_results_lncrna) == 0) {
            print(paste0("No DE lncRNA in ", interpretation))
        } else {  
            lncrna_de_matrix <- classification(sig_de_results_lncrna, lncRNA_classification, results_dir)

            addWorksheet(OUT, paste0(interpretation, "_DE_lncRNAs"))
            writeData(OUT, rowNames = TRUE, sheet = paste0(interpretation, "_DE_lncRNAs"),
                      x = as.data.frame(lncrna_de_matrix))
        }
    }

    # Export the file
    saveWorkbook(OUT, file_name)
}


plot_pca_comps <- function(pca_data, comp1, comp2, results_dir, percentVar, norm_method, title) { 
	
	pc_var1 <- paste0("PC", comp1)
	pc_var2 <- paste0("PC", comp2)
	
	pcaPlot <- ggplot2::ggplot(pca_data, aes(!!sym(pc_var1), !!sym(pc_var2), color=samplePoint, shape=condition)) +
  		ggplot2::geom_point(size=3) +
  	ggplot2::xlab(paste0(pc_var1, ": ", percentVar[comp1], "% variance")) +
  	ggplot2::ylab(paste0(pc_var2, ": ", percentVar[comp2], "% variance")) + 
  	ggplot2::coord_fixed()
	
	ggplot2::ggsave(paste0(results_dir, "/PCA/pca_plot_", title, norm_method, "_", comp1, comp2, ".png"), 
       plot = pcaPlot, 
       height = 5, width = 5, units = "in", dpi = 600)
}

# Read the GTF file
#gtf_data <- read.table(gtf, sep="\t", header=FALSE, stringsAsFactors=FALSE,quote="")

# Filter rows where the third column is "transcript"
#transcript_rows <- gtf_data %>% filter(V3 == "transcript")

# Extract gene_id and gene_name from the attributes column (9th column)
#extract_attributes <- function(attribute_string, key) {
#  match <- str_match(attribute_string, paste0(key, ' "([^"]+)"'))
#  return(match[,2])
#}

#transcript_rows$gene_id <- sapply(transcript_rows$V9, extract_attributes, key="gene_id")
#transcript_rows$ref_gene_id <- sapply(transcript_rows$V9, extract_attributes, key="ref_gene_id")
#transcript_rows$gene_name <- sapply(transcript_rows$V9, extract_attributes, key="gene_name")

# Select relevant columns: gene_id and gene_name
#transcript_info <- transcript_rows %>% select(gene_id, gene_name)

# Remove duplicates (if any)
#transcript_info <- distinct(transcript_info)


# lncrna transcript list
lncrna_transcripts <- read.table(lncrna_list,  
                                 header = F, stringsAsFactors = F)$V1

#lncrna_transcripts_conservative <- read.table(
 #                                header = F, stringsAsFactors = F)$V1

# load gene GTF
suppressMessages(stringtie_txdb <- makeTxDbFromGFF(gtf))
#  "transcriptome_assembly/non_protein_coding_stringtie.gtf"
#)

# collapse lncRNA transcripts to genes 

# will return a warning: 'select()' returned 1:1 mapping between keys and columns. Therefore, suppress


suppressMessages(lncrna_genes_all <- AnnotationDbi::select(stringtie_txdb,
                                            keys = lncrna_transcripts,
                                            columns = c("TXCHROM", 
                                                        "GENEID", 
                                                        "TXSTRAND", 
                                                        "TXSTART", 
                                                        "TXEND"), 
                                            keytype = "TXNAME"))

# lncrna genes
lncrna_genes <- lncrna_genes_all[!is.na(lncrna_genes_all$GENEID), ]$GENEID

# add the FEELnc classifiation for lncRNA type
lncRNA_classification <- read.csv(lncrna_classification, 
                                  header = T, stringsAsFactors = F, row.names = 1)[, 2:10]


#lncRNA_classification$partnerRNA_gene_name <- transcript_rows[match(lncRNA_classification$partnerRNA_gene, transcript_rows$ref_gene_id), "gene_name"]


sample_information <- read.table(sample_info_file, header=TRUE, stringsAsFactors=FALSE)

sample_information$samplePoint <- as.factor(sample_information$samplePoint)
sample_information$condition <- as.factor(sample_information$condition)
sample_information$samplePoint <- relevel(sample_information$samplePoint, "0.025")
sample_information$condition <- relevel(sample_information$condition, "SC")

sample_information$sampleID <- sample_information$accession

#rownames(sample_information) <- sample_information$sampleID


# Generate a simple matrix for importing via IsoformSwitchAnalyzeR

isoform_design <- sample_information[,c("sampleID", "sample")]
#isoform_design <- sample_information[,c("accession", "sample")]
colnames(isoform_design) <- c("sampleID","condition")


# Import salmon quantifications
salmonQuant <- importIsoformExpression(
    parentDir = quant_dir, addIsofomIdAsColumn = FALSE, quiet=TRUE)



# Create an isoform switch list. As the package isn't being used for isoform switch analysis, but mainly to deal with stringtie issues, then won't try and correct for unwanted effects as leads to issue during the removal (matrix independency issues, limma:remove_batch_effects issues)

aSwitchList <- importRdata(
	isoformCountMatrix   = salmonQuant$counts,
	isoformRepExpression = salmonQuant$abundance,
	designMatrix         = isoform_design,
	isoformExonAnnoation = gtf, 
	detectUnwantedEffects = FALSE, quiet=TRUE)



#"transcriptome_assembly/stringtie.all.transcripts.gtf",
#export gene counts. Summarises transcript counts to gene level, where appropriate
df_gene <- extractGeneExpression(aSwitchList, addGeneNames = FALSE, addIdsAsColumns=FALSE)
df_gene <- extractGeneExpression(aSwitchList, addGeneNames = FALSE, addIdsAsColumns=FALSE,extractCounts=FALSE)
write.table(df_gene, file=tpm, sep="\t", quote=F, col.names=NA)
#write.table(df_gene, file="lncrna_annotation/TPM/TPM_ISAER.tsv", sep="\t", quote=F, col.names=NA)
#paste(results_dir, "/edgeR_normalized_counts.tsv",  sep="")


##################################
## Quality control - PCA plots and heatmaps 

all_dds <- DESeqDataSetFromMatrix(countData = round(df_gene),
                                 colData = sample_information,
                                 design = ~ samplePoint)

counts_matrix <- assay(all_dds)
dge <- DEFormats::as.DGEList(all_dds)

# Calculagte between-sample normalization factors (TMM) using edge R. 
dge <- edgeR::calcNormFactors(dge)
dge_cpm <- edgeR::cpm(dge, log=F, normalized.lib.sizes=T)
write.table(dge_cpm, file=paste(results_dir, "/edgeR_normalized_counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)

# Convert back to a DESeqDataSet to keep the same workflow. However, need to rename normalization column for vst to work with edgeR normalsiation. 
dge_dds <- DEFormats::as.DESeqDataSet(dge)
names(colData(dge_dds))[stringr::str_detect(names(colData(dge_dds)), "norm.factors")] <- "sizeFactor"

# vst transformation (log transformation)
dge_rld <- DESeq2::vst(dge_dds, blind=TRUE)

dge_rld_mat <- SummarizedExperiment::assay(dge_rld)

colnames(dge_rld_mat) <- sample_information$sample_replicate[match(colnames(dge_rld_mat),sample_information$sampleID)]

#dge_lncrna_rld_mat <- dge_rld_mat[lncrna_genes,]
dge_lncrna_rld_mat <- dge_rld_mat[rownames(dge_rld_mat) %in% lncrna_genes, ]




#Plot Heatmap for sample correlation

dge_rld_cor <- cor(dge_rld_mat)
#heatMap <- pheatmap::pheatmap(dge_rld_cor)
#png(paste(results_dir,"/heatmaps/edgeR_heatmap.png",sep=""), )
heatMap <- pheatmap::pheatmap(dge_rld_cor,filename=paste(results_dir,"/heatmaps/edgeR_heatmap.png",sep=""))
#dev.off()

# save the order for the other plots
row_order <- heatMap$tree_row$order
col_order <- heatMap$tree_col$order
# PCAs

dge_pca <- prcomp(t(dge_rld_mat))

# Create data frame with metadata and all PC values for input to ggplot
dge_pcaData <- cbind(sample_information, dge_pca$x)
dge_percentVar <- round(100 * summary(dge_pca)$importance[2,])

plot_pca_comps(dge_pcaData, 1,2,results_dir, dge_percentVar, "edgeR", "all")
plot_pca_comps(dge_pcaData, 2,3,results_dir, dge_percentVar, "edgeR", "all")
plot_pca_comps(dge_pcaData, 3,4,results_dir, dge_percentVar, "edgeR", "all")

### For lncRNAs now

# Plot lncRNA based heatmap

dge_lncrna_rld_cor <- cor(dge_lncrna_rld_mat)

# to get it as same order as first one but currently does not adjust the dendograms
heatMap <- pheatmap::pheatmap(dge_lncrna_rld_cor[row_order, col_order], 
		cluster_rows=F, cluster_cols=F,
		filename=paste(results_dir,"/heatmaps/edgeR_heatmap_lncrna.png",sep=""))

## PCA plots
dge_lncrna_pca <- prcomp(t(dge_lncrna_rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
dge_lncrna_pcaData <- cbind(sample_information, dge_lncrna_pca$x)
dge_lncrna_percentVar <- round(100 * summary(dge_lncrna_pca)$importance[2,])

plot_pca_comps(dge_lncrna_pcaData, 1,2,results_dir, dge_lncrna_percentVar, "edgeR", "lncrna")
plot_pca_comps(dge_lncrna_pcaData, 2,3,results_dir, dge_lncrna_percentVar, "edgeR", "lncrna")
plot_pca_comps(dge_lncrna_pcaData, 3,4,results_dir, dge_lncrna_percentVar, "edgeR", "lncrna")

######


# DESeq2 normalized results 

# calculate size factors (DESeq2 uses median of ratios method) and apply to DESeqDataSet
# Median of ratios accounts for sequencing depth and RNA
# Median of ratios = counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene

all_dds <- DESeq2::estimateSizeFactors(all_dds)
dds_cpm <- DESeq2::counts(all_dds, normalized=TRUE)
write.table(dds_cpm, file=paste(results_dir, "/DESeq_normalized_counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)


## Some QC
# blind =TRUE results in transformation unbiased to sample condition information
dds_rld <- DESeq2::vst(all_dds, blind=TRUE)

dds_rld_mat <- SummarizedExperiment::assay(dds_rld)

# Check if all values in each row are the same
row_not_all_same <- apply(dds_rld_mat, 1, function(row) !all(row == row[1]))

# Subset the matrix to keep only rows where not all values are the same
dds_rld_mat <- dds_rld_mat[row_not_all_same, , drop = FALSE]

colnames(dds_rld_mat) <- sample_information$sample_replicate[match(colnames(dds_rld_mat),sample_information$sampleID)]

dds_lncrna_rld_mat <- dds_rld_mat[rownames(dds_rld_mat) %in% lncrna_genes,]

#Plot Heatmap for sample correlation

dds_rld_cor <- cor(dds_rld_mat)
dds_heatMap <- pheatmap::pheatmap(dds_rld_cor, filename=paste(results_dir,"/heatmaps/DESeq2_heatmap.png",sep=""))

# save the order for the other plots
row_order <- dds_heatMap$tree_row$order
col_order <- dds_heatMap$tree_col$order

# PCAs

dds_pca <- prcomp(t(dds_rld_mat))

# Create data frame with metadata and all PC values for input to ggplot
dds_pcaData <- cbind(sample_information, dds_pca$x)
dds_percentVar <- round(100 * summary(dds_pca)$importance[2,])

plot_pca_comps(dds_pcaData, 1,2,results_dir, dds_percentVar, "DESeq2", "all")
plot_pca_comps(dds_pcaData, 2,3,results_dir, dds_percentVar, "DESeq2", "all")
plot_pca_comps(dds_pcaData, 3,4,results_dir, dds_percentVar, "DESeq2", "all")

# Plot lncRNA based heatmap
dds_lncrna_rld_cor <- cor(dds_lncrna_rld_mat)
# to get it as same order as first one but currently does not adjust the dendograms
dds_heatMap <- pheatmap::pheatmap(dds_lncrna_rld_cor[row_order, col_order], 
		cluster_rows=F, cluster_cols=F,
		filename=paste(results_dir,"/heatmaps/DESeq2_heatmap_lncrna.png",sep=""))

dds_lncrna_pca <- prcomp(t(dds_lncrna_rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
dds_lncrna_pcaData <- cbind(sample_information, dds_lncrna_pca$x)
dds_lncrna_percentVar <- round(100 * summary(dds_lncrna_pca)$importance[2,])

plot_pca_comps(dds_lncrna_pcaData, 1,2,results_dir, dds_lncrna_percentVar, "DESeq2", "lncrna")
plot_pca_comps(dds_lncrna_pcaData, 2,3,results_dir, dds_lncrna_percentVar, "DESeq2", "lncrna")
plot_pca_comps(dds_lncrna_pcaData, 3,4,results_dir, dds_lncrna_percentVar, "DESeq2", "lncrna")


#############################################
## Now want to look at a number of scenarios. 
# 1) identify genes/lncRNAs with a growth-associated expression behaviour.
# 2) Identify genes differentially expressed (DE) between different growth rates
# 3) Identify DE between the same sample point but between SC and LC

# 1 & 2: Will only use the samples from the SC retentostat (no 0.1 and LC). Thus, can avoid having to account for any batch related issues and also avoids any full rank issues as no R6 samplePoint for LC, and no 0.1 samplePoint for SC or LC condition

# Subset sample_information and gene counts

sc_information <- sample_information[grep("SC",sample_information$condition),]
sc_information <- droplevels(sc_information)
sc_gene <- df_gene[,sc_information$sampleID]
sc_gene <- df_gene[,sc_information$sampleID]

# Round it as counts from salmon not itnegrs. Done like this in tximport, so not a problem
sc_dds <- DESeqDataSetFromMatrix(countData = round(sc_gene),
                                 colData = sc_information,
                                 design = ~ samplePoint)

## 1) Identify subset of genes that are better explained by a full model accounting for the growth rate, than the intercept alone (reduced model)
# calculate differential expression using the DESeq wrapper function

suppressMessages(sc_dds_lrt <- DESeq(sc_dds, test="LRT", reduced=~1))

# set differential expression criteria
de_results <- results(sc_dds_lrt,
                      lfcThreshold = 0,
                      independentFiltering = T, alpha=0.01)

# 2) Now, want to compare between samplePoints using the "standard" wald test of DESeq
suppressMessages(sc_dds_wt <- DESeq(sc_dds))

desired_comparisons <- c("samplePoint_R3_vs_0.025", 
			 "samplePoint_R6_vs_0.025", 
			 "samplePoint_R10_vs_0.025")
comparison_interpretation <- c("R3_vs_0.025","R6_vs_0.025","R10_vs_0.025")

file_name <- paste(results_dir, "/SCRetentostat_DESeq2_results.xlsx",sep="")
suppressMessages(if (file.exists(file_name)) {file.remove(file_name)})

# generate excel files and volcano plots for desired comparisons
process_comparison(desired_comparisons, comparison_interpretation, sc_dds_wt, lncrna_genes, lncRNA_classification, file_name, results_dir) 


# 3). Include all SC and LC samples. Removed R6 comparison as there was no R6 samplePoint for the LC condition. Removed SS_01 as not in SC or LC


sc_lc_information <- sample_information[-grep("C0.1|R6", sample_information$samplePoint),]
sc_lc_information <- droplevels(sc_lc_information)
sc_lc_gene <- df_gene[,sc_lc_information$sampleID]

# Round it as counts from salmon not itnegrs. Done like this in tximport, so not a problem
sc_lc_dds <- DESeqDataSetFromMatrix(countData = round(sc_lc_gene),
                                 colData = sc_lc_information,
                                 design = ~ condition + samplePoint + condition:samplePoint)


suppressMessages(sc_lc_dds_wt <- DESeq(sc_lc_dds))


desired_comparisons <- c("conditionLC.samplePointR10" , "conditionLC.samplePointR3", "condition_LC_vs_SC")
comparison_interpretation <- c("R10_LC_vs_R10_SC", "R3_LC_vs_R3_SC", "LC_vs_SC")

file_name <- paste(results_dir, "/SCLC_Retentostat_DESeq2_results.xlsx",sep="")
suppressMessages(if (file.exists(file_name)) {file.remove(file_name)})

process_comparison(desired_comparisons, comparison_interpretation, sc_lc_dds_wt, lncrna_genes, lncRNA_classification, file_name, results_dir) 






