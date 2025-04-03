#!/usr/bin/Rscript
## -----------------------------------------------------------------------------
##
## Script name: deseq_analysis.R
##
## Purpose of script: To carry out count based differential expression analyis 
## for the full StringTie transcriptome assembly. Uses IsoformSwitchAnalyzeR to correct stringtie issues (merging etc).
## The differential expression results for  lncRNAs are parsed and the classififcaiton and orthology information is added 
##


suppressMessages(library("IsoformSwitchAnalyzeR"))
suppressMessages(library("DESeq2"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("DEFormats"))
suppressMessages(library("edgeR"))
suppressMessages(library("pheatmap"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("gprofiler2"))
suppressMessages(library("reshape2"))
suppressMessages(library("MPLNClust"))
suppressMessages(library("Mfuzz"))
suppressMessages(library("viridis"))
suppressMessages(library("tidyr"))
suppressMessages(library("cowplot"))
suppressMessages(library("GenomicFeatures"))

options(bitmapType='cairo')

# Input arguments
args <- commandArgs(TRUE)
quant_dir <- args[1]
sample_info_file <- args[2]
lncrna_classification <- args[3]
lncrna_list <- args[4]
gtf <- args[5]
de_dir <- args[6]
tpm <- args[7]
annotation_info <- args[8]
interacting_pairs_promoter <- args[9]
interacting_pairs_exon <- args[10]
interacting_pairs_all <- args[11]
neighbouring_pairs <- args[12]
func_dir <- args[13]
#
#quant_dir <- "data/retentostat/salmon/quant/"
#sample_info_file <- "data/retentostat_sample_information.tsv"
#lncrna_classification <-"lncrna_annotation/classification/simple_classification.txt"
#lncrna_list <-"lncrna_annotation/monoexonic_filter/first_lncrna_list.txt"
#gtf <-"transcriptome_assembly/stringtie.all.transcripts.gtf"
#de_dir <-"differential_expression/retentostat"
#tpm <-"lncrna_annotation/TPM/TPM_ISAER.tsv"
#annotation_info <-"data/kphaffi_annotation_info.tsv"
#interacting_pairs_promoter <-"functional_predictions/triplexes/interacting_pairs_promoter.txt"
#interacting_pairs_exon <-"functional_predictions/triplexes/interacting_pairs_exon.txt"
#interacting_pairs_all <-"functional_predictions/triplexes/interacting_pairs_all.txt"
#neighbouring_pairs <-"functional_predictions/nearby/nearby_genes.txt"
#func_dir <-"functional_predictions"

dir.create(file.path(de_dir, "/heatmaps"),recursive=TRUE)


dir.create(file.path(func_dir, "/clustering/mfuzz"),recursive=TRUE)
dir.create(file.path(func_dir, "/clustering/mpln"),recursive=TRUE)


volcano_plot <- function(de_results, svalue_threshold = 0.005) {
  # Identify significant results based on thresholds
  #lfc_threshold = log2(1.5), 
  de_results <- as.data.frame(de_results) %>%
    dplyr::mutate(
	  threshold_DE = svalue < svalue_threshold,
      # threshold_DE = padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold,
      pointcolor = dplyr::case_when(
        log2FoldChange > 0 & threshold_DE ~ "#f99988", ##F9B288
        log2FoldChange < 0 & threshold_DE ~ "#A2D9F9",
        TRUE ~ "gray"
      ),
      pointclass = dplyr::case_when(
        log2FoldChange > 0 & threshold_DE ~ "Upregulated",
        log2FoldChange < 0 & threshold_DE ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      pointsize = if_else(threshold_DE, 3, 1)
    )
	
  # Organize data and ensure factor levels for consistent plotting
  de_results$pointclass <- factor(de_results$pointclass, 
                                  levels = c("Downregulated", "Not Significant", "Upregulated"))
  
  # Return processed data
  # return(de_results %>% select(log2FoldChange, padj, pointcolor, pointclass, pointsize))
  return(de_results %>% select(log2FoldChange, svalue, pointcolor, pointclass, pointsize))
}

plot_pca_comps <- function(pca_data, comp1, comp2, percentVar) { 
  
    pc_var1 <- paste0("PC", comp1)
    pc_var2 <- paste0("PC", comp2)

    pcaPlot <- ggplot2::ggplot(pca_data, aes_string(x = pc_var1, y = pc_var2, color = "samplePoint", shape = "condition")) + #, shape = "condition")) +
        ggplot2::geom_point(size = 3) +
        ggplot2::xlab(paste0(pc_var1, ": ", percentVar[comp1], "% variance")) +
        ggplot2::ylab(paste0(pc_var2, ": ", percentVar[comp2], "% variance")) + 
        ggplot2::coord_fixed() +
        theme_bw() + 
        scale_color_discrete(name = "Sample Point") +
        scale_shape_discrete(name = "Condition") + #, labels=levels(pca_data$samplePoint)) +  
        theme(legend.key.size = unit(2, "lines"),   # Increase legend marker size
              legend.text = element_text(size = 14),  # Increase legend text size
              legend.title = element_blank(), #element_text(element_blank(), size = 14), # Increase legend title size
              axis.text.x = element_text(size = 14),  # Increase x-axis tick label size
              axis.text.y = element_text(size = 14),
              legend.spacing = unit(0.4, 'cm'),       # Increase spacing between title and markers
              legend.spacing.x = unit(0.4, 'cm'),
              panel.background = element_blank(),     # Remove grey background
              panel.grid.major = element_line(linetype = "dashed", color = "light gray"),  # Keep major gridlines
              panel.grid.minor = element_blank())   # Keep minor gridlines) 
 
    return(pcaPlot)
}


volcano_ggplot <- function(volcano_data, comparison, scm_values, x_limits, y_limits) {
  vgp <- ggplot(volcano_data, aes(x = log2FoldChange, 
                                   y = -log10(svalue), #padj
                                   color = pointclass, 
                                   size = pointsize)) +
    geom_point() + 
    scale_size_identity() + 
    scale_color_manual(name = "Colour", 
                       values = scm_values, 
                       limits = levels(droplevels(volcano_data$pointclass)), 
                       labels = levels(droplevels(volcano_data$pointclass))) +
    xlab(bquote(~ Log[2] ~ "fold change")) + 
    ylab("-log10 s-value") + #BH adjusted p-value
    geom_hline(yintercept = -log10(0.005),
               linetype = "dashed", color = "gray") + 
    # geom_vline(xintercept = c(-0.5849625, 0.5849625),
               # linetype = "dashed",
               # color = "light gray") +  
    xlim(x_limits) +
    ylim(y_limits) +
    theme_bw() + guides(colour=guide_legend(ncol=3, override.aes = list(size = 3))) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(1, "lines"),   # Increase legend marker size
      legend.text = element_text(size = 14),  # Increase legend text size
      legend.spacing = unit(0.4, 'cm'),       # Increase spacing between title and markers
      legend.spacing.x = unit(0.6, 'cm'),
      axis.text.x = element_text(size = 14),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7)
    )
  return(vgp)
}


classification <- function(sig_results, lncRNA_classification, id_name_description) {

  # Helper function to get matched values with NA replaced by "intergenic"
  get_class_value <- function(col_name) {
    match_value <- lncRNA_classification[match(rownames(sig_results), lncRNA_classification$lncRNA_gene), col_name]
    replace(match_value, is.na(match_value), "intergenic")
  }
  
  # Helper function for retrieving partner description
  get_partner_description <- function() {
    partner_ids <- get_class_value("partnerRNA_gene")
    partner_desc <- id_name_description$description[match(partner_ids, id_name_description$id)]
    replace(partner_desc, is.na(partner_desc), "intergenic")
  }
  
  # Create classification data
  lncrna_de_matrix <- data.frame(
    baseMean = sig_results$baseMean,
    log2FoldChange = sig_results$log2FoldChange,
    lfcSE = sig_results$lfcSE,
    svalue = sig_results$svalue,
    # pvalue = sig_results$pvalue,
    # padj = sig_results$padj,
    `lncRNA type` = get_class_value("summary_class"),
    `Match distance` = get_class_value("distance"),
    `Closest gene` = get_class_value("partnerRNA_gene"),
    `Closest gene name` = get_class_value("partnerRNA_gene_name"),
    `Closest gene description` = get_partner_description()
  )
  
  return(lncrna_de_matrix)
}


process_comparison <- function(desired_comparisons, comparison_interpretation, DESeqDataSet, lncrna_genes, lncRNA_classification, file_name, de_dir) {
    volcano_data_list <- list()
    all_lncrna_results_list <- list()  # List to store data frames
    OUT <- createWorkbook()
	
    for (i in 1:length(desired_comparisons)) {
        comparison <- desired_comparisons[i]
        interpretation <- comparison_interpretation[i]
		
		res <- results(DESeqDataSet, alpha=0.005, name=comparison)
		suppressMessages(res <- lfcShrink(DESeqDataSet, coef=comparison, type = "apeglm", res=res, svalue=TRUE))
		
		# res <- results(DESeqDataSet, alpha=0.1, name=comparison, lfcThreshold=log2(1.5), altHypothesis="greaterAbs")
		# suppressMessages(res <- lfcShrink(DESeqDataSet, coef=comparison, type = "apeglm", res=res))
		print(summary(res))
		
		sig_de_results <- subset(
			res,
			svalue < 0.005
			)
			
		sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing = TRUE),]
		
		# extract only lncRNA genes
		all_lncrna_results <- res[rownames(res) %in% lncrna_genes, ]
		lncrna_matrix <- classification(all_lncrna_results, lncRNA_classification, id_name_description)
		sig_de_results_lncrna <- sig_de_results[rownames(sig_de_results) %in% lncrna_genes, ]
		
		# Plot volcano of lncRNA
		volcano_data <- volcano_plot(all_lncrna_results)
		volcano_data_list[[interpretation]] <- volcano_data
		
		# Export all results - all transcripts and both DE and Non DE
		addWorksheet(OUT, paste0(interpretation, "_all"))
		sr_df <- as.data.frame(res)
		sr_df$ID <- rownames(sr_df)
		rownames(sr_df) <- NULL
		sr_df <- sr_df %>% select(ID, everything())
		sr_df$name <- id_name_description$name[match(rownames(res), id_name_description$id)]
		sr_df$description <- id_name_description$description[match(rownames(res), id_name_description$id)]
		writeData(OUT, rowNames = FALSE, sheet = paste0(interpretation, "_all"),
					  x = sr_df)			  
		# Separate the changing columns
        dynamic_cols <- lncrna_matrix[, c("log2FoldChange", "lfcSE", "svalue")]
        colnames(dynamic_cols) <- paste(interpretation, colnames(dynamic_cols), sep="_")
        constant_cols <- lncrna_matrix[, !(colnames(lncrna_matrix) %in% c("log2FoldChange", "lfcSE", "svalue"))]
        combined <- cbind(constant_cols, dynamic_cols)
        combined$Gene <- rownames(combined)
        rownames(combined) <- NULL
        all_lncrna_results_list[[i]] <- combined
		
        addWorksheet(OUT, paste0(interpretation, "_all_lncRNA"))
        writeData(OUT, rowNames = TRUE, sheet = paste0(interpretation, "_all_lncRNA"),
                  x = as.data.frame(lncrna_matrix))
	
        if (nrow(sig_de_results_lncrna) == 0) {
            print(paste0("No DE lncRNA in ", interpretation))
			} else { 
				lncrna_de_matrix <- classification(sig_de_results_lncrna, lncRNA_classification, id_name_description)
				addWorksheet(OUT, paste0(interpretation, "_DE_lncRNAs"))
				writeData(OUT, rowNames = TRUE, sheet = paste0(interpretation, "_DE_lncRNAs"), x = as.data.frame(lncrna_de_matrix))
			}
	}
    merge_columns <- c("Gene", "baseMean", "lncRNA.type", "Match.distance", "Closest.gene", "Closest.gene.name", "Closest.gene.description")
# Combine all data frames into one with interpretation-based column levels
    final_lncrna_results <- Reduce(function(x, y) merge(x, y, by =  merge_columns, all = TRUE), all_lncrna_results_list)
    # List of dynamic columns (log2FoldChange, lfcSE, svalue) from all iterations
    dynamic_cols_order <- setdiff(colnames(final_lncrna_results), merge_columns)
    final_col_order <- c("Gene", "baseMean", dynamic_cols_order, "lncRNA.type", "Match.distance", "Closest.gene", "Closest.gene.name", "Closest.gene.description")
    final_lncrna_results <- final_lncrna_results[, final_col_order]
    rownames(final_lncrna_results) <- final_lncrna_results$Gene
    final_lncrna_results$Gene <- NULL
    addWorksheet(OUT, "lncRNA_All_comparisons")
    writeData(OUT, rowNames = TRUE, sheet = "lncRNA_All_comparisons", x = as.data.frame(final_lncrna_results))
    saveWorkbook(OUT, file_name, overwrite=TRUE)
    return(volcano_data_list)
}



# Define a function to process the counts matrix and produce a median matrix
create_median_matrix <- function(counts_matrix, sig_de_results, non_lc_info, sample_info, sample_order) {
  # Subset matrix to include only significant DE genes and order by sample ID
  sig_matrix <- counts_matrix[rownames(sig_de_results), as.character(non_lc_info$sampleID)]
  
  # Rename columns to use sample replicate names for consistency
  colnames(sig_matrix) <- sample_info$sample_replicate[match(colnames(sig_matrix), sample_info$sampleID)]
  print(paste0("Significant dimensions ", paste0(dim(sig_matrix), collapse = ' x ')))
  
  # Transpose and add sample point information
  transposed <- as.data.frame(t(sig_matrix))
  transposed$samplePoint <- sample_info$samplePoint[match(rownames(transposed), sample_info$sample_replicate)]
  
  # Calculate median counts by sample point
  median_df <- aggregate(transposed[, 1:(ncol(transposed) - 1)], list(transposed$samplePoint), median)
  
  # Set row names and reorder by sample point
  rownames(median_df) <- median_df$Group.1
  median_df <- median_df[, -1]
  median_df <- median_df[match(sample_order, rownames(median_df)), ]
  
  # Convert to matrix and transpose for clustering or analysis
  median_matrix <- data.matrix(t(median_df))
  print(paste0("Final dimensions ", paste0(dim(median_matrix), collapse = ' x ')))
  
  return(median_matrix)
}

VisualizeLine <- function(dataset,
                              clusterMembershipVector = NA,
                              fileName = paste0('Plot_', Sys.Date()),
                              output_dir = "") {

  coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
  
  # Log transform the dataset
  dataset <- log(dataset + 1)

  # Combine log-transformed dataset with cluster membership vector
  DataPlusLabs <- cbind(dataset, cluster=clusterMembershipVector)

  # Convert DataPlusLabs to a dataframe
  df <- as.data.frame(DataPlusLabs)
  df$genes <- rownames(df)
  
  df_long <- reshape2::melt(df, id.vars = c("genes","cluster"))
  cluster_means <- df_long %>% group_by(variable, cluster) %>% summarise(mean_value = mean(value, na.rm = TRUE))
  cluster_means <- as.data.frame(cluster_means)


  # Plot using ggplot2 with facet_grid

  png(file = paste0(func_dir, "/clustering/mpln/",fileName, ".png"),width=800, height=600)
  p <- ggplot(data = df_long, aes(x = variable, y = value, color = as.factor(cluster), group = genes)) +
  geom_line() +
  geom_line(data = df_long %>% group_by(variable, cluster) %>% summarise(mean_value = mean(value, na.rm = TRUE)),
            aes(group = cluster, x = variable, y = mean_value),color = "black") +
  scale_color_manual(values = coloursBarPlot) +
  labs(x = "Samples", y = "Expression (log counts)")+
  theme_minimal() + facet_wrap(~ cluster, ncol = 4) + 
  theme(legend.position="none",axis.text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text = element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
   
  print(p)
  dev.off()
}

#############################################################################################################
### END OF FUNCTIONS

# lncrna transcript list
lncrna_transcripts <- read.table(lncrna_list,  
                                 header = F, stringsAsFactors = F)$V1

# load gene GTF
suppressMessages(stringtie_txdb <- makeTxDbFromGFF(gtf))

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
lncrna_genes <- unique(lncrna_genes)


# add the FEELnc classifiation for lncRNA type
lncRNA_classification <- read.csv(lncrna_classification, 
                                  header = T, stringsAsFactors = F, row.names = 1)[, 2:10]

id_name_description <- read.table(annotation_info, header=TRUE, quote="", sep="\t")


sample_information <- read.table(sample_info_file, header=TRUE, stringsAsFactors=FALSE)
sample_information$samplePoint <- as.factor(sample_information$samplePoint)
sample_information$condition <- as.factor(sample_information$condition)
sample_information$samplePoint <- relevel(sample_information$samplePoint, "0.025")
sample_information$condition <- relevel(sample_information$condition, "SC")
sample_information$sampleID <- sample_information$accession


# Generate a simple matrix for importing via IsoformSwitchAnalyzeR

isoform_design <- sample_information[,c("sampleID", "sample")]
colnames(isoform_design) <- c("sampleID","condition")

# Import salmon quantifications
salmonQuant <- importIsoformExpression(
    parentDir = quant_dir, addIsofomIdAsColumn = FALSE, quiet=TRUE)


# Create an isoform switch list. As the package isn't being used for isoform switch 
# analysis, but mainly to deal with stringtie issues, then won't try and correct 
# for unwanted effects as leads to issue during the removal (matrix independency 
# issues, limma:remove_batch_effects issues)

aSwitchList <- importRdata(
	isoformCountMatrix   = salmonQuant$counts,
	isoformRepExpression = salmonQuant$abundance,
	designMatrix         = isoform_design,
	isoformExonAnnoation = gtf, 
	detectUnwantedEffects = FALSE, quiet=TRUE)


#export gene counts. Summarises transcript counts to gene level, where appropriate
df_gene <- extractGeneExpression(aSwitchList, addGeneNames = FALSE, addIdsAsColumns=FALSE,extractCounts=FALSE)
write.table(df_gene, file=tpm, sep="\t", quote=F, col.names=NA)

##################################
## Quality control - PCA plots and heatmaps 
non_lc_information <- sample_information[-grep("LC",sample_information$condition),]
non_lc_information <- droplevels(non_lc_information)
non_lc_gene <- df_gene[,non_lc_information$sampleID]

# Round it as counts from salmon not itnegrs. Done like this in tximport, so not a problem
non_lc_dds <- DESeqDataSetFromMatrix(countData = round(non_lc_gene),
                                 colData = non_lc_information,
                                 design = ~ samplePoint)
								 

######

###############################	
# Normalisation and transformations
###############################

# DESeq2 normalized results 

# calculate size factors (DESeq2 uses median of ratios method) and apply to DESeqDataSet
# Median of ratios accounts for sequencing depth and RNA
# Median of ratios = counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene

non_lc_dds <- DESeq2::estimateSizeFactors(non_lc_dds)
dds_cpm <- DESeq2::counts(non_lc_dds, normalized=TRUE)
write.table(dds_cpm, file=paste(de_dir, "/DESeq_normalized_counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)


## Some QC
# blind =TRUE results in transformation unbiased to sample condition information
dds_rld <- DESeq2::vst(non_lc_dds, blind=TRUE)
dds_rld_mat <- SummarizedExperiment::assay(dds_rld)

# Check if all values in each row are the same
row_not_all_same <- apply(dds_rld_mat, 1, function(row) !all(row == row[1]))

# Subset the matrix to keep only rows where not all values are the same
dds_rld_mat <- dds_rld_mat[row_not_all_same, , drop = FALSE]
norm_counts_matrix <- dds_rld_mat

colnames(dds_rld_mat) <- non_lc_information$sample_replicate[match(colnames(dds_rld_mat), non_lc_information$sampleID)]
dds_lncrna_rld_mat <- dds_rld_mat[rownames(dds_rld_mat) %in% lncrna_genes,]

###############################	
# Heatmap
###############################

# correlation for all transcripts
dds_rld_cor <- cor(dds_rld_mat)
# for lncrna
dds_lncrna_rld_cor <- cor(dds_lncrna_rld_mat)

# Plot heatmap for all and extract the row and col order
dds_heatMap <- pheatmap::pheatmap(dds_rld_cor, filename=paste(de_dir,"/heatmaps/DESeq2_heatmap.png",sep=""))
row_order <- dds_heatMap$tree_row$order
col_order <- dds_heatMap$tree_col$order

# Plot lncRNA based heatmap
# to get it as same order as first one but currently does not adjust the dendograms
dds_heatMap <- pheatmap::pheatmap(dds_lncrna_rld_cor[row_order, col_order], 
		cluster_rows=F, cluster_cols=F,
		filename=paste(de_dir,"/heatmaps/DESeq2_heatmap_lncrna.png",sep=""))

###############################	
# PCAs
###############################
dds_pca <- prcomp(t(dds_rld_mat))
# Create data frame with metadata and all PC values for input to ggplot
dds_pcaData <- cbind(non_lc_information, dds_pca$x)
dds_pcaData$samplePoint <- factor(dds_pcaData$samplePoint, levels = c("C0.1", "0.025", "R3", "R6", "R10"))
dds_percentVar <- round(100 * summary(dds_pca)$importance[2,])


dds_lncrna_pca <- prcomp(t(dds_lncrna_rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
dds_lncrna_pcaData <- cbind(non_lc_information, dds_lncrna_pca$x)
dds_lncrna_pcaData$samplePoint <- factor(dds_lncrna_pcaData$samplePoint, levels = c("C0.1", "0.025", "R3", "R6", "R10"))
dds_lncrna_percentVar <- round(100 * summary(dds_lncrna_pca)$importance[2,])

## Plots are generated later after DE analysis

###############################	
# Pairwise wald test for DE
###############################


sc_information <- sample_information[grep("SC",sample_information$condition),]
sc_information <- droplevels(sc_information)
sc_gene <- df_gene[,sc_information$sampleID]
sc_gene <- df_gene[,sc_information$sampleID]

# Round it as counts from salmon not itnegrs. Done like this in tximport, so not a problem
sc_dds <- DESeqDataSetFromMatrix(countData = round(sc_gene),
                                 colData = sc_information,
                                 design = ~ samplePoint)



# 2) Now, want to compare between samplePoints using the "standard" wald test of DESeq
suppressMessages(sc_dds_wt <- DESeq(sc_dds))

desired_comparisons <- c("samplePoint_R3_vs_0.025", 
			 "samplePoint_R6_vs_0.025", 
			 "samplePoint_R10_vs_0.025")
comparison_interpretation <- c("R3_vs_0.025","R6_vs_0.025","R10_vs_0.025")

file_name <- paste(de_dir, "/SCRetentostat_DESeq2_results.xlsx",sep="")
suppressMessages(if (file.exists(file_name)) {file.remove(file_name)})

# generate excel files and volcano plots for desired comparisons
# Call the process_comparison function
volcano_data_list <- process_comparison(desired_comparisons, comparison_interpretation, sc_dds_wt, lncrna_genes, lncRNA_classification, file_name, de_dir)


###############################	
# Plotting volcanos and PCAs
###############################


# Combine volcano data for global axis limits
all_volcano_data <- do.call(rbind, volcano_data_list)
global_x_limits <- c(-max(abs(all_volcano_data$log2FoldChange), na.rm=TRUE), max(abs(all_volcano_data$log2FoldChange), na.rm=TRUE))

# global_y_limits <- c(0, max(-log10(all_volcano_data$padj), na.rm = TRUE))
global_y_limits <- c(0, max(-log10(all_volcano_data$svalue), na.rm = TRUE))


# Generate ggplot objects with consistent axes
volcano_plots <- lapply(names(volcano_data_list), function(interpretation) {
    volcano_data <- volcano_data_list[[interpretation]]
    scm_values <- unique(volcano_data$pointcolor)
    names(scm_values) <- unique(volcano_data$pointclass)
    volcano_ggplot(volcano_data, interpretation, scm_values, global_x_limits, global_y_limits)
})

# Extract individual plots
volcano_plot_1 <- volcano_plots[[1]]
volcano_plot_2 <- volcano_plots[[2]]
volcano_plot_3 <- volcano_plots[[3]]


# Generate PCA plots
pca_plot_1 <- plot_pca_comps(dds_pcaData, 1, 2, dds_percentVar)
pca_plot_2 <- plot_pca_comps(dds_lncrna_pcaData, 1,2, dds_lncrna_percentVar)


# Combine PCA plots with a shared legend
pca_combined <- plot_grid(pca_plot_1 + theme(legend.position = "none") + theme(aspect.ratio = 2/3),
                          pca_plot_2 + theme(legend.position = "none") + theme(aspect.ratio = 2/3),
                          labels = c("A", "B"),
                          ncol = 2, align = 'hv')

# Extract legend from one of the PCA plots
pca_legend <- get_legend(pca_plot_1 + theme(legend.position = "top",
                                            legend.title = element_text(size = 14),
                                            legend.text = element_text(size = 14)))

# Combine the PCA combined plot with its legend
pca_combined_with_legend <- plot_grid(pca_combined, pca_legend, ncol = 1, rel_heights = c(1, 0.08))

# Combine Volcano plots with a shared legend
volcano_combined <- plot_grid(volcano_plot_1 + theme(legend.position = "none"),
                              volcano_plot_2 + theme(legend.position = "none"),
                              volcano_plot_3 + theme(legend.position = "none"),
                              labels = c("C", "D", "E"),
                              ncol = 3, align = 'hv')

# Extract legend from one of the Volcano plots
volcano_legend <- get_legend(volcano_plot_2 + theme(legend.position = "bottom"))

# Combine the Volcano combined plot with its legend
volcano_combined_with_legend <- plot_grid(volcano_combined, volcano_legend, ncol = 1, rel_heights = c(1, 0.08))

# Arrange the final figure with PCA on top and Volcano plots at the bottom
five_panel_figure <- plot_grid(pca_combined_with_legend,
                               volcano_combined_with_legend,
                               ncol = 1, rel_heights = c(1, 1))

# Save the combined figure
ggsave(paste0(de_dir, "/FivePanel_PCA_Volcano.png"), plot = five_panel_figure, width = 12, height = 10)

#################################################################################
# Functional analysis
#################################################################################

###############################
# Process count matrices - median etc for clustering and correlation
###############################


counts_matrix <- SummarizedExperiment::assay(non_lc_dds)
print(paste0("Original dimensions ", paste0(dim(counts_matrix), collapse=' x ')))

## Identify subset of genes that are better explained by a full model accounting for the growth rate
# , than the intercept alone (reduced model)
# calculate differential expression using the DESeq wrapper function

suppressMessages(non_lc_dds_lrt <- DESeq(non_lc_dds, test="LRT", reduced=~1))
de_results <- results(non_lc_dds_lrt,
                      lfcThreshold = 0,
                      independentFiltering = T, alpha=0.0005)
print(summary(de_results))

sig_de_results <- subset(de_results,  padj < 0.0005) #, baseMean >= 100)
DE_genes <- rownames(sig_de_results)

print("There are x lncRNAs")
print(length(lncrna_genes))

print("x of them are DE")
print(length(lncrna_genes[lncrna_genes %in% DE_genes]))


# Use the function with the original counts matrix
counts_median_matrix <- create_median_matrix(
  counts_matrix = counts_matrix, 
  sig_de_results = sig_de_results, 
  non_lc_info = non_lc_information, 
  sample_info = sample_information, 
  sample_order = c("C0.1", "0.025", "R3", "R6", "R10")
)

# Use the function with the transformed matrix
transformed_median_matrix <- create_median_matrix(
  counts_matrix = norm_counts_matrix, 
  sig_de_results = sig_de_results, 
  non_lc_info = non_lc_information, 
  sample_info = sample_information, 
  sample_order = c("C0.1", "0.025", "R3", "R6", "R10")
)



###############################
# Zscore heatmaps of neighbouring and interacting 
###############################

nonlc_zscore_matrix <- t(scale(t(dds_rld_mat)))
nonlc_zscore_matrix <- na.omit(nonlc_zscore_matrix)


triplex_pairs_prom <- read.table(interacting_pairs_promoter)
colnames(triplex_pairs_prom) <- c("Source", "Target")

triplex_pairs_exon <- read.table(interacting_pairs_exon)
colnames(triplex_pairs_exon) <- c("Source", "Target")

triplex_pairs_all <- read.table(interacting_pairs_all)
colnames(triplex_pairs_all) <- c("Source", "Target")

neighbour_pairs <- read.table(neighbouring_pairs)
colnames(neighbour_pairs) <- c("Source", "Target")


# Not all Source and Target present due to overlap issues (~250 out of 7500)
triplex_pairs_all <- triplex_pairs_all[(triplex_pairs_all$Source %in% rownames(nonlc_zscore_matrix))&(triplex_pairs_all$Target %in% rownames(nonlc_zscore_matrix)),]

# Not all Source and Target present due to overlap issues (~250 out of 7500)
triplex_pairs_prom <- triplex_pairs_prom[(triplex_pairs_prom$Source %in% rownames(nonlc_zscore_matrix))&(triplex_pairs_prom$Target %in% rownames(nonlc_zscore_matrix)),]

# Not all Source and Target present due to overlap issues (~250 out of 7500)
triplex_pairs_exon <- triplex_pairs_exon[(triplex_pairs_exon$Source %in% rownames(nonlc_zscore_matrix))&(triplex_pairs_exon$Target %in% rownames(nonlc_zscore_matrix)),]

# Not all Source and Target present due to overlap issues (~250 out of 7500)
neighbour_pairs <- neighbour_pairs[(neighbour_pairs$Source %in% rownames(nonlc_zscore_matrix))&(neighbour_pairs$Target %in% rownames(nonlc_zscore_matrix)),]


## Create heatmap annotations - top and bottom

df <- data.frame(gene = c(rep("lncRNA", 15), rep("Coding", 15)))
ann_labels=list(gene = list(title = "", at = c("lncRNA", "Coding"), labels=c("lncRNA", "Coding")))
ann_col = list(gene = c("lncRNA" = "CornflowerBlue", "Coding" = "firebrick2"))
t_ha = HeatmapAnnotation(df = df, col = ann_col, annotation_legend_param = ann_labels, show_annotation_name = FALSE)

bottom_df <- data.frame(sample = rep(c(rep("C0.1", 3), rep("0.025", 3),rep("R3", 3),rep("R6", 3), rep("R10", 3)),2))
ann_labels=list(sample = list(title = "", at = c("C0.1", "0.025", "R3", "R6", "R10"), labels = c("C0.1", "0.025", "R3", "R6", "R10")))
ann_col = list(sample = c("C0.1" = "grey1", "0.025" = "grey25","R3" = "grey50","R6" = "grey75","R10" = "grey90"))
b_ha = HeatmapAnnotation(df = bottom_df, col = ann_col, annotation_legend_param = ann_labels, show_annotation_name = FALSE)


## Plot heatmaps
paired_znorm <- cbind(nonlc_zscore_matrix[triplex_pairs_all$Source,],nonlc_zscore_matrix[triplex_pairs_all$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/triplexes/NonLC_TriplexAllZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()


## Plot heatmaps
paired_znorm <- cbind(nonlc_zscore_matrix[triplex_pairs_prom$Source,],nonlc_zscore_matrix[triplex_pairs_prom$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/triplexes/NonLC_TriplexPromZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()


## Plot heatmaps
paired_znorm <- cbind(nonlc_zscore_matrix[triplex_pairs_exon$Source,],nonlc_zscore_matrix[triplex_pairs_exon$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/triplexes/NonLC_TriplexExonZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()






paired_znorm <- cbind(nonlc_zscore_matrix[neighbour_pairs$Source,],nonlc_zscore_matrix[neighbour_pairs$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/nearby/NonLC_NeighbourZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()

###############################
# Correlation coeffiecient and density plots
###############################

# Calculate correlation matrices
cor_matrix <- cor(t(nonlc_zscore_matrix), method = "pearson")
lncrna_cor_matrix <- cor_matrix[rownames(cor_matrix) %in% unique(lncrna_genes),]
lncrna_cor_matrix <- lncrna_cor_matrix[rowSums(is.na(lncrna_cor_matrix)) != ncol(lncrna_cor_matrix) - 1,]

# Save lncRNA correlation matrix
write.table(lncrna_cor_matrix, file = paste0(func_dir, "/lncrna_correlation_matrix.csv"), quote = FALSE, row.names = TRUE, sep = ",")

# Extract significant correlations
correlated <- lapply(rownames(lncrna_cor_matrix), function(lncrna) {
  names(lncrna_cor_matrix[lncrna, abs(lncrna_cor_matrix[lncrna,]) > 0.8])
})
names(correlated) <- rownames(lncrna_cor_matrix)

# Save correlation results
cor_export <- data.frame(x = names(correlated), correlated_partners = I(unlist(lapply(correlated, paste, collapse = ","))))
write.table(cor_export, file = paste0(func_dir, "/correlation_results.csv"), quote = FALSE, row.names = FALSE, sep = ",")
write.table(cor_matrix, file = paste0(func_dir, "/correlation_matrix.csv"), quote = FALSE, row.names = TRUE, sep = ",")

# Melt the lncRNA correlation matrix
all <- melt(lncrna_cor_matrix)
colnames(all) <- c("Source", "Target", "cor")

# Expand and filter the 'all' dataframe
all_expanded <- all %>%
  separate_rows(Source, sep = "\\|") %>%
  separate_rows(Target, sep = "\\|") %>%
  filter(Source != Target)

# Generate interacting and neighboring data
interacting_all <- merge(all_expanded, triplex_pairs_all, by = c("Source", "Target")) %>% filter(Source != Target)

interacting_prom <- merge(all_expanded, triplex_pairs_prom, by = c("Source", "Target")) %>% filter(Source != Target)

interacting_exon <- merge(all_expanded, triplex_pairs_exon, by = c("Source", "Target")) %>% filter(Source != Target)


neighbours <- merge(all_expanded, neighbour_pairs, by = c("Source", "Target")) %>% filter(Source != Target)

# Identify independent entries in 'all_expanded' that are not in 'neighbours' or 'interacting'
independent <- all_expanded %>%
  anti_join(interacting_prom, by = c("Source", "Target")) %>%
  anti_join(neighbours, by = c("Source", "Target"))

# anti_join(interacting_all, by = c("Source", "Target")) %>%


  
# Combine all data types, including the new 'Independent' category
all_neighbour_interacting <- rbind(
  #data.frame(all_expanded, Type = "All"),
  data.frame(neighbours, Type = "Neighbouring"),
  data.frame(interacting_prom, Type = "Interacting - Promoter"),
  data.frame(independent, Type = "Independent")
)


# Plot density distributions for All, Neighbouring, and Interacting types
png(file = paste0(func_dir, "/All_Neighbour_Interacting_DC.png"), units = 'cm', res = 450, width = 5 * 450 / 72, height = 5 * 450 / 72)
ggplot(all_neighbour_interacting, aes(x=cor, group = Type, col = Type)) +
  theme_classic(base_size = 25) +
  geom_density(linewidth = 2, bounds = c(-1, 1)) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 0.7), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 12))) +
  labs(
    #title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
    y = "Density", x = "Correlation Coefficient"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("Independent" = "#F8766D", "Neighbouring" = "#619CFF", "Interacting - Promoter" = "#00BA38"))
dev.off()

#y=after_stat(density*n/nrow(all_neighbour_interacting))

# Density scaling and plotting for DE and NDE categories
plot_scaled_density <- function(data, type_filter, file_name, title, colors) {
  # Calculate scaling factor to align densities
  group1_density <- density(data$cor[data$Type == type_filter])
  group2_density <- density(data$cor[data$Type == "lncRNA DE"])
  group3_density <- density(data$cor[data$Type == "lncRNA NDE"])
  scaling_factor <- sum(group1_density$y) / (sum(group2_density$y) + sum(group3_density$y))
  
  # Generate plot
  png(file = file_name, units = 'cm', res = 450, width = 5 * 450 / 72, height = 5 * 450 / 72)
  i <- ggplot(data, aes(cor, group = Type, col = Type)) +
    theme_classic(base_size = 25) +
    geom_density(data = subset(data, Type == type_filter), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(data, Type == "lncRNA DE"), aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(data, Type == "lncRNA NDE"), aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 0.7), expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 12))) +
    labs(
      #title = title,
      y = "Density", x = "Correlation Coefficient"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors)
  print(i)
  dev.off()
}



neighbour_de <- rbind(
    data.frame(neighbours, Type = "Neighbours"),
    data.frame(neighbours[neighbours$Source %in% DE_genes, ], Type = "lncRNA DE"),
    data.frame(neighbours[!neighbours$Source %in% DE_genes, ], Type = "lncRNA NDE"))

# Plot Neighbor and Interacting DE/NDE densities
plot_scaled_density(
  neighbour_de,
  "Neighbours",
  file_name = paste0(func_dir, "/nearby/Neighbours_Neig_All_lncRNA_DE_NDE.png"),
  title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
  colors = c("Neighbours" = "#619CFF", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange")
)


interacting_gp_de <- rbind(
    data.frame(interacting_all, Type = "Interacting - Gene + Promoter"),
    data.frame(interacting_all[interacting_all$Source %in% DE_genes, ], Type = "lncRNA DE"),
    data.frame(interacting_all[!interacting_all$Source %in% DE_genes, ], Type = "lncRNA NDE"))

plot_scaled_density(
  interacting_gp_de,
  "Interacting - Gene + Promoter",
  file_name = paste0(func_dir, "/triplexes/Triplexes_Int_All_lncRNA_DE_NDE.png"),
  title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
  colors = c("Interacting - Gene + Promoter" = "#848384", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange")
)

interacting_prom_de <- rbind(
    data.frame(interacting_prom, Type = "Interacting - Promoter"),
    data.frame(interacting_prom[interacting_prom$Source %in% DE_genes, ], Type = "lncRNA DE"),
    data.frame(interacting_prom[!interacting_prom$Source %in% DE_genes, ], Type = "lncRNA NDE"))

plot_scaled_density(
  interacting_prom_de,
  "Interacting - Promoter",
  file_name = paste0(func_dir, "/triplexes/Triplexes_Int_Prom_lncRNA_DE_NDE.png"),
  title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
  colors = c("Interacting - Promoter" = "#00BA38", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange")
)


interacting_exon_de <- rbind(
    data.frame(interacting_exon, Type = "Interacting - Exon"),
    data.frame(interacting_exon[interacting_exon$Source %in% DE_genes, ], Type = "lncRNA DE"),
    data.frame(interacting_exon[!interacting_exon$Source %in% DE_genes, ], Type = "lncRNA NDE"))

plot_scaled_density(
  interacting_exon_de,
  "Interacting - Exon",
  file_name = paste0(func_dir, "/triplexes/Triplexes_Int_Exon_lncRNA_DE_NDE.png"),
  title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
  colors = c("Interacting - Exon" = "#41c5bb", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange")
)


independent_de <- rbind(
    data.frame(independent, Type = "Independent"),
    data.frame(independent[independent$Source %in% DE_genes, ], Type = "lncRNA DE"),
    data.frame(independent[!independent$Source %in% DE_genes, ], Type = "lncRNA NDE"))

plot_scaled_density(
  independent_de,
  "Independent",
  file_name = paste0(func_dir, "/Independent_All_lncRNA_DE_NDE.png"),
  title = "Density distribution of correlation coefficients of \nlncRNAs and other genes",
  colors = c("Independent" = "#F8766D", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange")
)


######



interacting_overview <- rbind(
	data.frame(interacting_all[interacting_all$Source %in% DE_genes, ], Type = "Interacting - Gene + Promoter"),
	data.frame(interacting_exon[interacting_exon$Source %in% DE_genes, ], Type = "Interacting - Exon"),
	data.frame(interacting_prom[interacting_prom$Source %in% DE_genes, ], Type = "Interacting - Promoter"))



# Calculate scaling factor to align densities
group1_density <- density(interacting_overview$cor[interacting_overview$Type == "Interacting - Gene + Promoter"])
group2_density <- density(interacting_overview$cor[interacting_overview$Type == "Interacting - Promoter"])
group3_density <- density(interacting_overview$cor[interacting_overview$Type == "Interacting - Exon"])
scaling_factor <- sum(group1_density$y) / (sum(group2_density$y) + sum(group3_density$y))

colors = c("Interacting - Gene + Promoter" = "#848384", "Interacting - Promoter" = "#00BA38", "Interacting - Exon" = "#41c5bb")

# Generate plot
png(file = paste0(func_dir, "/triplexes/Triplexes_Int_Overview_lncRNA_DE_NDE.png"), units = 'cm', res = 450, width = 5 * 450 / 72, height = 5 * 450 / 72)
ggplot(interacting_overview, aes(cor, group = Type, col = Type)) +
theme_classic(base_size = 25) +
geom_density(data = subset(interacting_overview, Type == "Interacting - Gene + Promoter"), linewidth = 2, bounds = c(-1, 1)) +
geom_density(data = subset(interacting_overview, Type == "Interacting - Promoter"), aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
geom_density(data = subset(interacting_overview, Type == "Interacting - Exon"), aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
coord_cartesian(xlim = c(-1, 1), expand = FALSE) + #ylim = c(0, 0.7), 
guides(color = guide_legend(override.aes = list(size = 12))) +
labs(y = "Density", x = "Correlation Coefficient") +
theme(plot.title = element_text(hjust = 0.5)) +
scale_color_manual(values = colors)
dev.off()


###############################	
# MFuzz clustering
###############################

# eset <- ExpressionSet(assayData=median_matrix)
eset <- ExpressionSet(assayData=transformed_median_matrix)
eset.r <- filter.NA(eset, thres=0.25)
eset.f <- fill.NA(eset.r, mode='knn')
eset.s <- standardise(eset.f)

m1 <- mestimate(eset.s)

# 1.17 for all samples
# 2.01 for median

# Sometimes fails
eset.s <- filter.NA(eset.s)

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/clustering/mfuzz/cselection.png"))
tmp <- cselection(eset.s,m=m1,crange=seq(2,20,1),repeats=5,visu=TRUE)
dev.off()

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/clustering/mfuzz/Dmin.png"))
tmp <- Dmin(eset.s,m=m1,crange=seq(2,20,1),repeats=10,visu=TRUE)#
dev.off()


opt_c <- 4
cl <- mfuzz(eset.s,c=opt_c,m=m1)

colo <- viridis_pal()(10)

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/clustering/mfuzz/clusters_opt.png"))
mfuzz.plot(eset.s,cl=cl,mfrow=c(2,ceiling(opt_c/2)),colo=colo,new.window=FALSE, time.labels=colnames(exprs(eset.s)))
dev.off()

png(width=20,height=20,units="cm",res=400, filename=paste0(func_dir, "/clustering/mfuzz/clusters_opt_overlap.png"))
O <- overlap(cl)
Ptmp <- overlap.plot(cl,over=O,thres=0.1)
dev.off()

cores <- acore(eset.s, cl, min.acore=0.5)

# Create a new Excel workbook
output_wb <- createWorkbook()

# Loop through the list of data frames (cores)
for (i in 1:length(cores)) {
  cluster_data <- cores[[i]]  # Get the current cluster's data frame
  cluster_label <- paste("Cluster", i)  # Create a unique sheet name for the cluster

# Extract gene IDs as a vector and split them
  cluster_genes <- unlist(lapply(as.character(cluster_data$NAME), function(gene_id) {
  parts <- unlist(strsplit(gene_id, "\\|"))
  if (length(parts) > 1) {
    return(parts[-1])  # Remove the first part if it contains "|"
  } else {
    return(gene_id)  # Keep the original gene ID if it doesn't contain "|"
  }
   }))

 # Perform gene ontology analysis using gprofiler2's gost function for GO enrichment
# upload_GMT_file(gmtfile="reference_genome/GO_analysis/gprofiler_converted.gmt")
  GO_name <- "gp__Uku7_mqlG_9fI"  
  suppressMessages(GO_result <- gost(cluster_genes, organism = GO_name))
  
  # Perform gene ontology analysis using gprofiler2's gost function for the kegg ontology
# upload_GMT_file(gmtfile="reference_genome/GO_analysis/Kegg_Ontology.gmt")
  kegg_name <- "gp__PEAO_RSQ8_tq0" 
  suppressMessages(kegg_result <- gost(cluster_genes, organism = kegg_name))

# Create a list to store the results for this cluster
  cluster_result <- list(
    cluster_label = cluster_label,
    genes = cluster_genes,
    GO_enrichment = NA,
    kegg_enrichment = NA)

 # Check if both results are empty, if not, fill in the data frame
  if (length(GO_result$result) > 0) {
    cluster_result$GO_enrichment <- GO_result$result
    num_empty_rows <- nrow(cluster_result$GO_enrichment)
  } else { num_empty_rows <- 0 }

  
  if (length(kegg_result$result) > 0) {
    cluster_result$kegg_enrichment <- kegg_result$result
  }

# Create a new sheet with the cluster label as the sheet name
  addWorksheet(output_wb, sheetName = cluster_label)

# Write the genes to the first column
  writeData(output_wb, sheet = cluster_label, cluster_result$genes, startCol = 1, startRow = 1)

  writeData(output_wb, sheet = cluster_label, "GO enrichment", 
              startCol = 5, startRow = 1)

  # Write the GO_enrichment data to the sheet
  writeData(output_wb, sheet = cluster_label, cluster_result$GO_enrichment, startCol = 5, startRow = 2)
  writeData(output_wb, sheet = cluster_label, "Kegg enrichment", 
              startCol = 5, startRow = num_empty_rows + 5)

  # Write the kegg_enrichment data after the empty rows
  writeData(output_wb, sheet = cluster_label, cluster_result$kegg_enrichment, 
            startCol = 5, startRow = num_empty_rows + 6)
  
}


addWorksheet(output_wb, sheetName = "lncRNA membership")
membership <- cbind(cluster_label=cl$cluster, cl$membership)
lncrna_membership <- membership[rownames(membership) %in% unique(lncrna_genes),]
writeData(output_wb, sheet = "lncRNA membership", lncrna_membership, 
              startCol = 1, startRow = 1, rowNames=TRUE)

# Save the workbook to the specified file
saveWorkbook(output_wb, file =paste0(func_dir, "/clustering/mfuzz/fuzzy_cluster_data.xlsx"), overwrite=TRUE)


###############################	
# MPLN clustering
###############################

mplnResults <- mplnVariational(dataset=counts_median_matrix, 
				membership="none", 
				gmin=1, gmax=20, 
				initMethod="kmeans", 
				nInitIterations=3, 
				normalize="Yes")
save(mplnResults, file=paste0(func_dir, "/clustering/mpln/sig_mplnResults.RData"))




########################
# Process MPLN Clust results
#########################

# Can load the earlier run!
#load(file=paste0(func_dir, "/clustering/mpln/sig_mplnResults.RData"))


png(file=paste0(func_dir, "/clustering/mpln/significant_loglikelihood.png"), width=800, height=600)
par(mfrow = c(1, 2))
graphics::matplot(mplnResults$logLikelihood, xlab = "Run",
                  ylab = "logL", type = c("b"), pch = 1, lty = 2) 
ICvalues <- matrix(c(mplnResults$BICresults$allBICvalues,
              mplnResults$ICLresults$allICLvalues,
              mplnResults$AICresults$allAICvalues,
              mplnResults$AIC3results$allAIC3values),
              ncol=4) 
graphics::matplot(ICvalues, xlab = "Run", ylab = "Information criteria value", 
                  type = c("b"), pch = 1, col = 1:4) 
legend("top", inset = c(- 0.4, 0), legend = c("BIC", "ICL", "AIC", "AIC3"), 
        col = 1:4, pch = 1, horiz = TRUE, bty = "n")
dev.off()

print(summary(warnings()))

opt_clusters <- 7


clusters <- data.frame(ID=rownames(mplnResults$dataset),Cluster=mplnResults$allResults[[opt_clusters]]$clusterlabels)

# Initialize an empty data frame to store the processed data
processed_clusters <- data.frame()

# Process each row in the input data
for (i in 1:nrow(clusters)) {
  # Split the first column by "|"
  ids <- strsplit(as.character(clusters[i, 1]), "\\|")[[1]]
  n <- length(ids)
  
  if (n == 1) {
    # If only one element, keep it as is
    processed_clusters <- rbind(processed_clusters, clusters[i, ])
  } else {
    # If multiple elements, create multiple rows with the same value in the second column
    for (j in 2:n) {
      new_row <- data.frame(ids[j], clusters[i, 2])
      colnames(new_row) <- colnames(clusters)
      processed_clusters <- rbind(processed_clusters, new_row)
    }
  }
}


# Get unique cluster labels
unique_clusters <- sort(unique(processed_clusters[, 2]))

# Create an Excel workbook to store the data
output_wb <- createWorkbook()

# Loop through unique cluster labels
for (cluster_id in unique_clusters) {

  cluster_label <- paste("Cluster", cluster_id)  # Create a unique sheet name
  
  # Filter the data for the current cluster label
  cluster_data <- subset(processed_clusters, processed_clusters[, 2] == cluster_id)
  
  # Extract gene IDs as a vector
  cluster_genes <- cluster_data[, 1]
  
  # Perform gene ontology analysis using gprofiler2's gost function for GO enrichment
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/gprofiler_converted.gmt")
  GO_name <- "gp__Uku7_mqlG_9fI"  
  GO_result <- gost(cluster_genes, organism = GO_name)
  
  # Perform gene ontology analysis using gprofiler2's gost function for the kegg ontology
# upload_GMT_file(gmtfile="functional_predictions/GO_analysis/Kegg_Ontology.gmt")
  kegg_name <- "gp__PEAO_RSQ8_tq0" 
  kegg_result <- gost(cluster_genes, organism = kegg_name)

  # Create a list to store the results for this cluster
  cluster_result <- list(
    cluster_label = cluster_label,
    genes = cluster_genes,
    GO_enrichment = NA,
    kegg_enrichment = NA
  )
 # Check if both results are empty, if not, fill in the data frame
  if (length(GO_result$result) > 0) {
    cluster_result$GO_enrichment <- GO_result$result
    num_empty_rows <- nrow(cluster_result$GO_enrichment)
  } else { num_empty_rows <- 0 }

  
  if (length(kegg_result$result) > 0) {
    cluster_result$kegg_enrichment <- kegg_result$result
  }

# Create a new sheet with the cluster label as the sheet name
  addWorksheet(output_wb, sheetName = cluster_label)

# Write the genes to the first column
  writeData(output_wb, sheet = cluster_label, cluster_result$genes, startCol = 1, startRow = 1)
  writeData(output_wb, sheet = cluster_label, "GO enrichment", 
              startCol = 5, startRow = 1)

  # Write the GO_enrichment data to the sheet
  writeData(output_wb, sheet = cluster_label, cluster_result$GO_enrichment, startCol = 5, startRow = 2)

  # Add two empty rows after the GO_enrichment data
  for (i in 1:2) {
    writeData(output_wb, sheet = cluster_label, "", 
              startCol = 5, startRow = num_empty_rows + i + 2)
  }

  writeData(output_wb, sheet = cluster_label, "Kegg enrichment", 
              startCol = 5, startRow = num_empty_rows + 5)

  # Write the kegg_enrichment data after the empty rows
  writeData(output_wb, sheet = cluster_label, cluster_result$kegg_enrichment, 
            startCol = 5, startRow = num_empty_rows + 6)


}

optimal <- data.frame(mplnResults$allResults[[opt_clusters]]$clusterlabels)
colnames(optimal) <- "Optimal"
rownames(optimal) <- rownames(mplnResults$dataset)

membership <- as.data.frame(mplnResults$allResults[[opt_clusters]]$probaPost)
colnames(membership) <- paste("Cluster", 1:ncol(membership), sep = "")
rownames(membership) <- rownames(mplnResults$dataset)

membership <- cbind(optimal, membership)

addWorksheet(output_wb, sheetName = "lncRNA membership")
lncrna_membership <- membership[rownames(membership) %in% unique(lncrna_genes),]
writeData(output_wb, sheet = "lncRNA membership", lncrna_membership, 
              startCol = 1, startRow = 1, rowNames=TRUE)

# Save the workbook to the specified file
saveWorkbook(output_wb, file =paste0(func_dir, "/clustering/mpln/mpln_cluster_data.xlsx"), overwrite=TRUE)


##################
#### Plot MPLN Clust results
##################


MPLNVisuals <- VisualizeLine(dataset = mplnResults$dataset,
	clusterMembershipVector =mplnResults$allResults[[opt_clusters]]$clusterlabels,
	fileName = '/Line',
	output_dir=paste0(func_dir, "/clustering/mpln"))

setwd(paste0(func_dir,"/clustering/mpln/"))

mplnViz <- mplnVisualizeBar(vectorObservations=1:nrow(mplnResults$dataset), 
			probabilities=mplnResults$allResults[[opt_clusters]]$probaPost,
			clusterMembershipVector=mplnResults$allResults[[opt_clusters]]$clusterlabels,
			format = 'png',fileName = 'PosteriorProbabilities',
			printPlot=TRUE)

mplnViz <- mplnVisualizeHeatmap(dataset=mplnResults$dataset, 
			clusterMembershipVector=mplnResults$allResults[[opt_clusters]]$clusterlabels,
			printPlot=TRUE,fileName = 'Heatmap',
			format = 'png')

png(file="Alluvial_AlternativeClustering.png", width=800, height=600)
MPLNVisualsAlluvial <- mplnVisualizeAlluvial(nObservations = nrow(mplnResults$dataset),
			firstGrouping = mplnResults$BICresults$BICmodelSelectedLabels,
			secondGrouping = mplnResults$ICLresults$ICLmodelSelectedLabels,
			thirdGrouping =  mplnResults$AICresults$AICmodelSelectedLabels,
			fourthGrouping = mplnResults$AIC3results$AIC3modelSelectedLabels, 
			fileName = "Alluvial",
			printPlot = FALSE) 
dev.off()


print(summary(warnings()))



