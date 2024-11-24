
suppressMessages(library("IsoformSwitchAnalyzeR"))
suppressMessages(library("DESeq2"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("gprofiler2"))
suppressMessages(library("reshape2"))

#options(bitmapType='cairo')

# Input arguments
args <- commandArgs(TRUE)
quant_dir <- args[1]
sample_info_file <- args[2]
lncrna_classification <- args[3]
lncrna_list <- args[4]
# "lncrna_annotation/monoexonic_filter/conservative_lncrna_list.txt", 
gtf <- args[5]
interacting_pairs <- args[6]
neighbouring_pairs <- args[7]
outdir <- args[8]


# Set default directories (for testing purposes)
#quant_dir <- "data/retentostat/salmon/quant/"
#sample_info_file <- "data/retentostat_sample_information.tsv"
#lncrna_classification <- "lncrna_annotation/classification/simple_classification.txt" 
#lncrna_list <- "lncrna_annotation/monoexonic_filter/first_lncrna_list.txt"
## "lncrna_annotation/monoexonic_filter/conservative_lncrna_list.txt", 
#gtf <- "transcriptome_assembly/stringtie.all.transcripts.gtf"
#interacting_pairs <- "functional_predictions/triplexes/interacting_pairs.txt"
#neighbouring_pairs <- "functional_predictions/nearby/nearby_genes.txt"
#outdir <- "functional_predictions"

system(paste("mkdir -p ", outdir, "/clustering/mfuzz", sep = ""))
system(paste("mkdir -p ", outdir, "/clustering/mpln", sep = ""))


# lncrna transcript list
lncrna_transcripts <- read.table(lncrna_list,  
                                 header = F, stringsAsFactors = F)$V1

#lncrna_transcripts_conservative <- read.table(
 #                                header = F, stringsAsFactors = F)$V1

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

# add the FEELnc classifiation for lncRNA type
lncRNA_classification <- read.csv(lncrna_classification, 
                                  header = T, stringsAsFactors = F, row.names = 1)[, 2:10]



sample_information <- read.table(sample_info_file, header=TRUE, stringsAsFactors=FALSE)

sample_information$samplePoint <- as.factor(sample_information$samplePoint)
sample_information$condition <- as.factor(sample_information$condition)
sample_information$samplePoint <- relevel(sample_information$samplePoint, "0.025")
sample_information$condition <- relevel(sample_information$condition, "SC")

sample_information$sampleID <- sample_information$accession

#rownames(sample_information) <- sample_information$sampleID


# Generate a simple matrix for importing via IsoformSwitchAnalyzeR

isoform_design <- sample_information[,c("sampleID", "sample")]
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

##################################

all_dds <- DESeqDataSetFromMatrix(countData = round(df_gene), colData = sample_information, design = ~ samplePoint)

#####


## Some QC
# blind =TRUE results in transformation unbiased to sample condition information
dds_rld <- DESeq2::vst(all_dds, blind=TRUE)

dds_rld_mat <- SummarizedExperiment::assay(dds_rld)

colnames(dds_rld_mat) <- sample_information$sample_replicate[match(colnames(dds_rld_mat),sample_information$sampleID)]
condition <- sample_information$condition[match(colnames(dds_rld_mat),sample_information$sample_replicate)]


# Check if all values in each row are the same
row_not_all_same <- apply(dds_rld_mat, 1, function(row) !all(row == row[1]))

# Subset the matrix to keep only rows where not all values are the same
dds_rld_mat <- dds_rld_mat[row_not_all_same, , drop = FALSE]


print("Clear")


dds_lncrna_rld_mat <- dds_rld_mat[lncrna_genes,]

# Remove rows of zeros
#dds_rld_mat <- dds_rld_mat[rowMax(dds_rld_mat) > 0,]



### Find DE
non_lc_information <- sample_information[-grep("LC", sample_information$condition),]
non_lc_information <- droplevels(non_lc_information)
non_lc_gene <- df_gene[,non_lc_information$sampleID]

# Round it as counts from salmon not itnegrs. Done like this in tximport, so not a problem
non_lc_dds <- DESeqDataSetFromMatrix(countData = round(non_lc_gene),
                                 colData = non_lc_information,
                                 design = ~ samplePoint)

counts_matrix <- SummarizedExperiment::assay(non_lc_dds)

## Identify subset of genes that are better explained by a full model accounting for the growth rate, than the intercept alone (reduced model)
# calculate differential expression using the DESeq wrapper function

suppressMessages(non_lc_dds_lrt <- DESeq(non_lc_dds, test="LRT", reduced=~1))

# set differential expression criteria
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

sig_nonlc_matrix <- counts_matrix[rownames(sig_de_results), c(as.character(non_lc_information$sampleID))]


## Calculate median counts between samples
transposed <- as.data.frame(t(sig_nonlc_matrix))

colnames(sig_nonlc_matrix) <- sample_information$sample_replicate[match(colnames(sig_nonlc_matrix), sample_information$sampleID)]

transposed$samplePoint <- sample_information$samplePoint[match(rownames(transposed), sample_information$sampleID)]
median_df <- aggregate(transposed[,1:length(colnames(transposed))-1], list(transposed$samplePoint), median)
# rename rows accoridng to smapling point and drop smapling poitn columns after
rownames(median_df) <- median_df$Group.1
median_df <- median_df[,-1]
# reorder on decreasing gr
median_df <- median_df[match(c("C0.1", "0.025", "R3", "R6", "R10"),rownames(median_df)),]
# finally convert abck to matrix
median_matrix <- data.matrix(t(median_df))



print(paste0("Original dimensions ", paste0(dim(counts_matrix), collapse=' x ')))
#print(paste0("Reduced by var and mean dimensions ", paste0(dim(reduced_matrix), collapse=' x ')))
print(paste0("Significant dimensions ", paste0(dim(sig_nonlc_matrix), collapse=' x ')))
print(paste0("Final dimensions ", paste0(dim(median_matrix), collapse=' x ')))

#DE <- read.table("differential_expression/results/DE_R10_vs_SS_0.025.tsv")
#DE_genes <- rownames(DE)
#####

nonlc_dds_rld_mat <- dds_rld_mat[,condition %in% c("C","SC")]
nonlc_zscore_matrix <- t(scale(t(nonlc_dds_rld_mat)))
nonlc_zscore_matrix <- na.omit(nonlc_zscore_matrix)


all_zscore_matrix <- t(scale(t(dds_rld_mat)))
all_zscore_matrix <- na.omit(all_zscore_matrix)

# Some rows contain NaNs, remove
#zscore_matrix <- zscore_matrix[!rowSums(is.na(zscore_matrix)),] 

triplex_pairs <- read.table(interacting_pairs)
colnames(triplex_pairs) <- c("Source", "Target")
# Not all Source and Target present due to overlap issues (~250 out of 7500)
#triplex_pairs <- triplex_pairs[(triplex_pairs$Source %in% rownames(all_zscore_matrix))&(triplex_pairs$Target %in% rownames(all_zscore_matrix)),]

neighbour_pairs <- read.table(neighbouring_pairs)
colnames(neighbour_pairs) <- c("Source", "Target")
# Not all Source and Target present due to overlap issues (~250 out of 7500)
#neighbour_pairs <- neighbour_pairs[(neighbour_pairs$Source %in% rownames(all_zscore_matrix))&(neighbour_pairs$Target %in% rownames(all_zscore_matrix)),]

#df <- data.frame(type = c(rep("lncRNA", 24), rep("Coding", 24)))

#paired_znorm <- cbind(all_zscore_matrix[triplex_pairs$Source,],all_zscore_matrix[triplex_pairs$Target,])

#png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/triplexes/All_TriplexZScoreCorrelation.png"))
#ha = HeatmapAnnotation(df = df, col = list(type = c("lncRNA" = "cyan3", "Coding" = "cyan4")), show_annotation_name = TRUE)
#ht <- Heatmap(as.matrix(paired_znorm), top_annotation = ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
#ComplexHeatmap::draw(ht)
#dev.off()

#paired_znorm <- cbind(all_zscore_matrix[neighbour_pairs$Source,],all_zscore_matrix[neighbour_pairs$Target,])

#png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/nearby/All_NeighbourZScoreCorrelation.png"))
#ha = HeatmapAnnotation(df = df, col = list(type = c("lncRNA" = "cyan3", "Coding" = "cyan4")), show_annotation_name = TRUE)
#ht <- Heatmap(as.matrix(paired_znorm), top_annotation = ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
#draw(ht)
#dev.off()


# Not all Source and Target present due to overlap issues (~250 out of 7500)
triplex_pairs <- triplex_pairs[(triplex_pairs$Source %in% rownames(nonlc_zscore_matrix))&(triplex_pairs$Target %in% rownames(nonlc_zscore_matrix)),]

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
paired_znorm <- cbind(nonlc_zscore_matrix[triplex_pairs$Source,],nonlc_zscore_matrix[triplex_pairs$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/triplexes/NonLC_TriplexZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()

paired_znorm <- cbind(nonlc_zscore_matrix[neighbour_pairs$Source,],nonlc_zscore_matrix[neighbour_pairs$Target,])

png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/nearby/NonLC_NeighbourZScoreCorrelation.png"))
ht <- Heatmap(as.matrix(paired_znorm), top_annotation = t_ha, bottom_annotation = b_ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
ComplexHeatmap::draw(ht)
dev.off()

#df <- data.frame(gene = c(rep("lncRNA", 15), rep("Coding", 15)), sample = rep(c(rep("C0.1", 3), rep("0.025", 3),rep("R3", 3),rep("R6", 3), rep("R10", 3)),2))
#ann_labels=list(sample = list(title = "sample",
#			    at = c("C0.1", "0.025", "R3", "R6", "R10"),
#			    labels = c("C0.1", "0.025", "R3", "R6", "R10")), 
#		gene = list(title = "gene",
#			    at = c("lncRNA", "Coding"), 
#			    labels=c("lncRNA", "Coding")))
#ann_colors = list(sample = c("C0.1" = "grey1", "0.025" = "grey20","R3" = "grey40","R6" = "grey60","R10" = "grey80"),
#		gene=c("lncRNA" = "DarkGreen", "Coding" = "DarkKhaki"))
#ha = HeatmapAnnotation(df = df, col = ann_colors, annotation_legend_param = ann_labels)

#png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/triplexes/NonLC_TriplexZScoreCorrelation.png"))
#ht <- Heatmap(as.matrix(paired_znorm), top_annotation = ha, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, col = cm.colors(256), show_row_names = FALSE, name = "z-score", column_title="Correlation of expression levels for lncRNA and gene pairs")
#ComplexHeatmap::draw(ht)
#dev.off()




#CPM.Tr.M <- rm_0_rows

#### Correlation density

cor_matrix <- cor(t(nonlc_zscore_matrix), method = "spearman")
lncrna_cor_matrix <- cor_matrix[rownames(cor_matrix) %in% unique(lncrna_genes),]

correlated <- vector(mode="list", length=length(rownames(lncrna_cor_matrix)))
names(correlated) <- rownames(lncrna_cor_matrix)

lncrna_cor_matrix <- lncrna_cor_matrix[rowSums(is.na(lncrna_cor_matrix)) != ncol(lncrna_cor_matrix) -1,]

write.table(lncrna_cor_matrix,file=paste0(outdir, "/lncrna_correlation_matrix.csv"),quote=FALSE,row.names=TRUE, sep=",")


for(lncrna in rownames(lncrna_cor_matrix)){
	correlated[[lncrna]] <- names(lncrna_cor_matrix[lncrna,(abs(lncrna_cor_matrix[lncrna,]) >0.8) == TRUE])}


cor_export <- data.frame(x=names(correlated),correlated_partners=I(unlist(lapply(correlated, paste, collapse=","))))

write.table(cor_export,file=paste0(outdir, "/correlation_results.csv"),quote=FALSE,row.names=F, sep=",")

write.table(cor_matrix,file=paste0(outdir, "/correlation_matrix.csv"),quote=FALSE,row.names=TRUE, sep=",")

all <- melt(lncrna_cor_matrix)
colnames(all)<-c("Source", "Target", "cor")


suppressMessages(library("tidyr"))
# Expand the 'all' dataframe
all_expanded <- all %>%
  separate_rows(Source, sep = "\\|") %>%
  separate_rows(Target, sep = "\\|")

all_expanded <- all_expanded %>%
  dplyr::filter(Source != Target)

interacting <- merge(all_expanded,triplex_pairs, by=c("Source", "Target"))
interacting <- interacting %>%
  dplyr::filter(Source != Target)


neighbours <- merge(all_expanded,neighbour_pairs, by=c("Source", "Target"))
neighbours <- neighbours %>%
  dplyr::filter(Source != Target)



######

all_neighbour <-rbind(data.frame(all_expanded, Type="All"), data.frame(neighbours, Type="Neighbouring"))
all_interacting <-rbind(data.frame(all_expanded, Type="All"), data.frame(interacting, Type="Interacting"))
all_neighbour_interacting <-rbind(data.frame(all_expanded, Type="All"), data.frame(neighbours, Type="Neighbouring"), data.frame(interacting, Type="Interacting"))

# Plot neighbour, interacting and all 

png(file = paste0(outdir, "/All_Neighbour_Interacting_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_neighbour_interacting, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + coord_cartesian(xlim=c(-1, 1), ylim=c(0, 0.75), expand=FALSE) + #xlim(-1,1)  + ylim(0,0.75) + 
guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Plot neighbour and all 
png(file = paste0(outdir, "/nearby/Neighbours_All_All_Neighbour_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_neighbour, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Plot interacting and all 
png(file = paste0(outdir, "/triplexes/Triplexes_All_All_Int_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_interacting, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1) + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()


####################
neighbour_de <- rbind(data.frame(neighbours, Type="Neighbours"), data.frame(neighbours[neighbours$Source %in% DE_genes,], Type="lncRNA DE"), data.frame(neighbours[!neighbours$Source %in% DE_genes,], Type="lncRNA NDE"))
png(file = paste0(outdir, "/nearby/Neighbours_NeigAll_lncRNA_DE.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(neighbour_de, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1) + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5)) +
scale_color_manual(values = c("Neighbours" = "#619CFF", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange"))
dev.off()

# Calculate density scaling factors
group1_density <- density(neighbour_de$cor[neighbour_de$Type == "Neighbours"])
group2_density <- density(neighbour_de$cor[neighbour_de$Type == "lncRNA DE"])
group3_density <- density(neighbour_de$cor[neighbour_de$Type == "lncRNA NDE"])

# Sum of densities for Group2 and Group3 should equal Group1
scaling_factor <- sum(group1_density$y) / (sum(group2_density$y) + sum(group3_density$y))
png(file = paste0(outdir, "/nearby/Neighbours_Neig_All_lncRNA_DE_NDE.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
# Plot with geom_density and scaling applied
ggplot(neighbour_de, aes(cor, group = Type, col = Type)) + 
    theme_classic(base_size = 25) + 
    geom_density(data = subset(neighbour_de, Type == "Neighbours"), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(neighbour_de, Type == "lncRNA DE"),
                 aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(neighbour_de, Type == "lncRNA NDE"),
                 aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    coord_cartesian(xlim=c(-1, 1), ylim=c(0, 0.75), expand=FALSE) +
#xlim(-1, 1) + 
 #   ylim(0, 0.75) +
    guides(color = guide_legend(override.aes = list(size = 12))) +
    labs(title = paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"),
         y = "Density", x = "Correlation Coefficient") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("Neighbours" = "#619CFF", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange"))
dev.off()

interacting_de <- rbind(data.frame(interacting, Type="Interacting"), data.frame(interacting[interacting$Source %in% DE_genes,], Type="lncRNA DE"), data.frame(interacting[!interacting$Source %in% DE_genes,], Type="lncRNA NDE"))
png(file = paste0(outdir, "/triplexes/Triplexes_Int_All_lncRNA_DE.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(interacting_de, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1) + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5)) +
scale_color_manual(values = c("Interacting" = "#00BA38", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange"))
dev.off()


# Calculate density scaling factors
group1_density <- density(interacting_de$cor[interacting_de$Type == "Interacting"])
group2_density <- density(interacting_de$cor[interacting_de$Type == "lncRNA DE"])
group3_density <- density(interacting_de$cor[interacting_de$Type == "lncRNA NDE"])

# Sum of densities for Group2 and Group3 should equal Group1
scaling_factor <- sum(group1_density$y) / (sum(group2_density$y) + sum(group3_density$y))
png(file = paste0(outdir, "/triplexes/Triplexes_Int_All_lncRNA_DE_NDE.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
# Plot with geom_density and scaling applied
ggplot(interacting_de, aes(cor, group = Type, col = Type)) + 
    theme_classic(base_size = 25) + 
    geom_density(data = subset(interacting_de, Type == "Interacting"), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(interacting_de, Type == "lncRNA DE"),
                 aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    geom_density(data = subset(interacting_de, Type == "lncRNA NDE"),
                 aes(y = after_stat(density * scaling_factor)), linewidth = 2, bounds = c(-1, 1)) +
    coord_cartesian(xlim=c(-1, 1), ylim=c(0, 0.75), expand=FALSE) +
#xlim(-1, 1) + 
    #ylim(0, 0.75) +
    guides(color = guide_legend(override.aes = list(size = 12))) +
    labs(title = paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"),
         y = "Density", x = "Correlation Coefficient") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("Interacting" = "#00BA38", "lncRNA DE" = "darkorchid", "lncRNA NDE" = "darkorange"))
dev.off()


##################
## Interacting specific plots

all_int_lncRNA_DE <-rbind(data.frame(all_expanded[all_expanded$Source %in% DE_genes,], Type="All"), data.frame(interacting[interacting$Source %in% DE_genes,], Type="Interacting"))

png(file = paste0(outdir, "/triplexes/Triplexes_DElncRNA_All_Int_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_int_lncRNA_DE, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \ndifferentially expressed lncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

all_int_Target_DE <-rbind(data.frame(all_expanded[all_expanded$Target %in% DE_genes,], Type="All"), data.frame(interacting[interacting$Target %in% DE_genes,], Type="Interacting"))

png(file = paste0(outdir, "/triplexes/Triplexes_DETarget_All_Int_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_int_Target_DE, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other differentially expressed genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

interacting <- data.frame(interacting, lncRNA="NDE",Partner="NDE",stringsAsFactors = FALSE)
interacting$Partner[interacting$Target %in% DE_genes] <- "DE"
interacting$lncRNA[interacting$Source %in% DE_genes] <- "DE"

# Create a new column for class
interacting <- interacting %>%
  dplyr::mutate(Class = paste(lncRNA, Partner, sep = "_"))

write.table(interacting,file=paste0(outdir, "/triplexes/triplex_pairs_correlation.csv"),quote=FALSE,row.names=TRUE, sep=",")


PartnerDE_interacting <- interacting[interacting$Partner == "DE",]

png(file=paste0(outdir, "/triplexes/Triplexes_DETarget_DE_NDE_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(data.frame(PartnerDE_interacting), aes(cor, group=lncRNA, col=lncRNA)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and differentially expressed interacting genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

lncRNADE_interacting <- interacting[interacting$lncRNA == "DE",]

png(file=paste0(outdir, "/triplexes/Triplexes_DElncRNA_DE_NDE_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(data.frame(lncRNADE_interacting), aes(cor, group=Partner, col=Partner)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\ndifferentially expressed lncRNAs and interacting genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()


#png(file = paste0(outdir, "/triplexes/TriplexesAllHHistogramCorrelation.png"))
#ggplot(data.frame(interacting, Type="Interacting"), aes(cor, group=Type, col=Type)) + 
#geom_histogram() + xlim(-1,1) +
#labs(title=paste("Density distribution of correlation coefficients of \nlncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5)) + theme_classic()
#dev.off()


#png(file=paste0(outdir, "/triplexes/TriplexesStackedDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(DE_interacting), aes(cor, group=lncRNA_DE, col=lncRNA_DE, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="stack", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and interacting differentially expressed genes"), y="Density", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#png(file=paste0(outdir, "/triplexes/TriplexesMarginalDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(DE_interacting), aes(cor, after_stat(count) ,group=lncRNA_DE, col=lncRNA_DE, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="stack", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and interacting differentially expressed genes"), y="Count", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#png(file=paste0(outdir, "/triplexes/TriplexesStackedConditionalDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(DE_interacting), aes(cor, after_stat(count) ,group=lncRNA_DE, col=lncRNA_DE, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="fill", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and interacting differentially expressed genes"), y="Count", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()

###### Neighbouring specific plots

all_lncRNA_DE <-rbind(data.frame(all_expanded[all_expanded$Source %in% DE_genes,], Type="All"), data.frame(neighbours[neighbours$Source %in% DE_genes,], Type="Neighbouring"))

png(file = paste0(outdir, "/nearby/Neighbours_DElncRNA_All_Neighbour_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_lncRNA_DE, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\ndifferentially expressed lncRNAs and other genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

all_Target_DE <-rbind(data.frame(all_expanded[all_expanded$Target %in% DE_genes,], Type="All"), data.frame(neighbours[neighbours$Target %in% DE_genes,], Type="Neighbouring"))

png(file = paste0(outdir, "/nearby/Neighbours_DETarget_All_Neighbour_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(all_Target_DE, aes(cor, group=Type, col=Type)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and other differentially expressed genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

neighbours <- data.frame(neighbours, lncRNA="NDE",Neighbour="NDE",stringsAsFactors = FALSE)

neighbours$Neighbour[neighbours$Target %in% DE_genes] <- "DE"
neighbours$lncRNA[neighbours$Source %in% DE_genes] <- "DE"

# Create a new column for class
neighbours <- neighbours %>%
  dplyr::mutate(Class = paste(lncRNA, Neighbour, sep = "_"))

write.table(neighbours,file=paste0(outdir, "/nearby/neighbouring_pairs_correlation.csv"),quote=FALSE,row.names=TRUE, sep=",")

NeighbourDE_neighbours <- neighbours[neighbours$Neighbour == "DE",]

png(file=paste0(outdir, "/nearby/Neighbours_DETarget_DE_NDE_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(data.frame(NeighbourDE_neighbours), aes(cor, group=lncRNA, col=lncRNA)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and differentially expressed neighbouring genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()

lncRNADE_neighbours <- neighbours[neighbours$lncRNA == "DE",]

png(file=paste0(outdir, "/nearby/Neighbours_DElncRNA_DE_NDE_DC.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
ggplot(data.frame(lncRNADE_neighbours), aes(cor, group=Neighbour, col=Neighbour)) + theme_classic(base_size = 25) + 
geom_density(linewidth=2, bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,0.75) + guides(color = guide_legend(override.aes = list(size = 12))) +
labs(title=paste("Density distribution of correlation coefficients of\ndifferentially expressed lncRNAs and neighbouring genes"), y="Density", x="Correlation Coefficient") + 
theme(plot.title = element_text(hjust = 0.5))
dev.off()


#png(file=paste0(outdir, "/nearby/NeighboursStackedDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(lncRNADE_neighbours), aes(cor, group=Neighbour, col=Neighbour, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="stack", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\ndifferentially expressed lncRNAs and neighbouring genes"), y="Density", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()



#png(file=paste0(outdir, "/nearby/NeighboursMarginalDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(DE_neighbours), aes(cor, after_stat(count) ,group=lncRNA_DE, col=lncRNA_DE, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="stack", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and neighbouring differentially expressed genes"), y="Count", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#png(file=paste0(outdir, "/nearby/NeighboursStackedConditionalDensityCorrelation.png"), units = 'cm',res = 450,width = 5*450/72, height = 5*450/72)
#ggplot(data.frame(DE_neighbours), aes(cor, after_stat(count) ,group=lncRNA_DE, col=lncRNA_DE, fill=lncRNA_DE)) + theme_classic(base_size = 25) + 
#geom_density(linewidth=2, position="fill", bounds=c(-1,1)) + xlim(-1,1)  + ylim(0,NA) + guides(color = guide_legend(override.aes = list(size = 12))) +
#labs(title=paste("Density distribution of correlation coefficients of\nlncRNAs and neighbouring differentially expressed genes"), y="Count", x="Correlation Coefficient") + 
#theme(plot.title = element_text(hjust = 0.5))
#dev.off()
####################################################
### MFuzz clustering ## Uses only DGE genes. 

suppressMessages(library("MPLNClust"))
suppressMessages(library("Mfuzz"))

#eset <- ExpressionSet(assayData=CPM.Tr.M)
#eset <- ExpressionSet(assayData=zscore_matrix)
eset <- ExpressionSet(assayData=median_matrix)
eset.r <- filter.NA(eset, thres=0.25)

#eset.f <- fill.NA(eset.r, mode="mean")
#eset.f <- fill.NA(eset.r, mode='wknn')
eset.f <- fill.NA(eset.r, mode='knn')

eset.s <- standardise(eset.f)

m1 <- mestimate(eset.s)

# 1.17 for all samples
# 2.01 for median

# Sometimes fails
eset.s <- filter.NA(eset.s)


png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/clustering/mfuzz/cselection.png"))

tmp <- cselection(eset.s,m=m1,crange=seq(2,20,1),repeats=5,visu=TRUE)
dev.off()

png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/clustering/mfuzz/Dmin.png"))
tmp <- Dmin(eset.s,m=m1,crange=seq(2,20,1),repeats=10,visu=TRUE)#
dev.off()


opt_c <- 5
# Define a custom color palette from blue to red
#blue_to_red_palette <- colorRampPalette(c("blue", "red"))

# Generate the color palette
#colo <- blue_to_red_palette(10)

suppressMessages(library(viridis))
colo <- viridis_pal()(10)


cl <- mfuzz(eset.s,c=opt_c,m=m1)
png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/clustering/mfuzz/clusters_opt.png"))
mfuzz.plot(eset.s,cl=cl,mfrow=c(2,ceiling(opt_c/2)),colo=colo,new.window=FALSE, time.labels=colnames(exprs(eset.s)))
dev.off()

png(width=20,height=20,units="cm",res=400, filename=paste0(outdir, "/clustering/mfuzz/clusters_opt_overlap.png"))
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


# Extract gene IDs as a vector and split them
  #cluster_genes <- unlist(strsplit(as.character(cluster_data$NAME), "\\|"))

  # Extract gene IDs as a vector
  #cluster_genes <- cluster_data$NAME

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
saveWorkbook(output_wb, file =paste0(outdir, "/clustering/mfuzz/fuzzy_cluster_data.xlsx"), overwrite=TRUE)

###### For all samples

#c5 <- mfuzz(eset.s,c=5,m=m1)
#png(width=20,height=20,units="cm",res=400, filename=paste0(results_dir,"/clustering/mfuzz/#clusters_5.png"))
#mfuzz.plot(eset.s,cl=c5,mfrow=c(3,2),new.window=FALSE, time.labels=colnames(exprs(eset.s)))
#dev.off()

#png(width=20,height=20,units="cm",res=400, filename=paste0(results_dir,"/clustering/mfuzz/clusters_5_overlap.png"))
#O <- overlap(c5)
#Ptmp <- overlap.plot(c5,over=O,thres=0.05)
#dev.off()

#c10 <- mfuzz(eset.s,c=10,m=m1)
#png(width=20,height=20,units="cm",res=400, filename=paste0(results_dir,"/clustering/mfuzz/clusters_10.png"))
#mfuzz.plot(eset.s,cl=c10,mfrow=c(5,5),new.window=FALSE, time.labels=colnames(exprs(eset.s)))
#dev.off()

#png(width=20,height=20,units="cm",res=400, filename=paste0(results_dir,"/clustering/mfuzz/clusters_10_overlap.png"))
#O <- overlap(c10)
#Ptmp <- overlap.plot(c10,over=O,thres=0.05)
#dev.off()


###### 


#### MPLN clustering

mplnResults <- mplnVariational(dataset=median_matrix, 
				membership="none", 
				gmin=1, gmax=20, 
				initMethod="kmeans", #random", 
				nInitIterations=3, 
				normalize="Yes")
save(mplnResults, file=paste0(outdir, "/clustering/mpln/sig_mplnResults.RData"))

png(file=paste0(outdir, "/clustering/mpln/significant_loglikelihood.png"), width=800, height=600)
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


